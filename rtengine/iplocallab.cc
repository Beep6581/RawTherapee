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
 *  2016 Jacques Desmis <jdesmis@gmail.com>
 *  2016 Ingo Weyrich <heckflosse@i-weyrich.de>

 */
#include <cmath>
#include <glib.h>
#include <glibmm.h>

#include "rtengine.h"
#include "improcfun.h"
#include "curves.h"
#include "gauss.h"
#include "iccstore.h"
#include "iccmatrices.h"
#include "color.h"
#include "rt_math.h"
#include "jaggedarray.h"
#ifdef _DEBUG
#include "mytime.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "cplx_wavelet_dec.h"

#define BENCHMARK
#include "StopWatch.h"

#define cliploc( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )

#define CLIPC(a) ((a)>-42000?((a)<42000?(a):42000):-42000)  // limit a and b  to 130 probably enough ?
#define CLIPL(x) LIM(x,0.f,40000.f) // limit L to about L=120 probably enough ?
#define CLIPLOC(x) LIM(x,0.f,32767.f)
#define CLIPLIG(x) LIM(x,0.f, 99.5f)
#define CLIPCHRO(x) LIM(x,0.f, 140.f)
#define CLIPRET(x) LIM(x,-99.5f, 99.5f)
#define CLIP1(x) LIM(x, 0.f, 1.f)

namespace
{

float calcLocalFactor (const float lox, const float loy, const float lcx, const float dx, const float lcy, const float dy, const float ach)
{
//elipse x2/a2 + y2/b2=1
//transition elipsoidal
//x==>lox y==>loy
// a==> dx  b==>dy

    float kelip = dx / dy;
    float belip = sqrt ((rtengine::SQR ((lox - lcx) / kelip) + rtengine::SQR (loy - lcy))); //determine position ellipse ==> a and b
    float aelip = belip * kelip;
    float degrad = aelip / dx;
    float ap = rtengine::RT_PI / (1.f - ach);
    float bp = rtengine::RT_PI - ap;
    return 0.5f * (1.f + xcosf (degrad * ap + bp)); //trigo cos transition

}

}

namespace rtengine
{
using namespace procparams;

extern const Settings* settings;

struct local_params {
    float yc, xc;
    float lx, ly;
    float lxL, lyT;
    float dxx, dyy;
    float iterat;
    int cir;
    float thr;
    int prox;
    int chro, cont, sens, sensh, senscb, sensbn, senstm;
    float ligh;
    int shamo, shdamp, shiter, senssha;
    double shrad;
    double rad;
    double stren;
    int trans;
    bool inv;
    bool curvact;
    bool invrad;
    bool invret;
    bool invshar;
    bool actsp;
    float str;
    int qualmet;
    int qualcurvemet;
    float noiself;
    float noiselc;
    float noisecf;
    float noisecc;
    float mulloc[5];
    float threshol;
    float strengt;
    float gamm;
    float esto;
    float scalt;
    float rewe;
    bool colorena;
    bool blurena;
    bool tonemapena;
    bool retiena;
    bool sharpena;
    bool cbdlena;
    bool denoiena;
};

static void calcLocalParams (int oW, int oH, const LocallabParams& locallab, struct local_params& lp)
{
    int w = oW;
    int h = oH;
    int circr = locallab.circrad;
    float streng = ((float)locallab.stren) / 100.f;
    float gam = ((float)locallab.gamma) / 100.f;
    float est = ((float)locallab.estop) / 100.f;
    float scal_tm = ((float)locallab.scaltm) / 10.f;
    float rewe = ((float)locallab.rewei);

    float thre = locallab.thres / 100.f;
    double local_x = locallab.locX / 2000.0;
    double local_y = locallab.locY / 2000.0;
    double local_xL = locallab.locXL / 2000.0;
    double local_yT = locallab.locYT / 2000.0;
    double local_center_x = locallab.centerX / 2000.0 + 0.5;
    double local_center_y = locallab.centerY / 2000.0 + 0.5;
    double local_dxx = locallab.proxi / 8000.0;//for proxi = 2==> # 1 pixel
    double local_dyy = locallab.proxi / 8000.0;
    float iterati = (float) locallab.proxi;
//    double local_dyy = locallab.proxi;

    if (locallab.qualityMethod == "std") {
        lp.qualmet = 0;
    } else if (locallab.qualityMethod == "enh") {
        lp.qualmet = 1;
    } else if (locallab.qualityMethod == "enhden") {
        lp.qualmet = 2;
    }

    if (locallab.qualitycurveMethod == "none") {
        lp.qualcurvemet = 0;
    } else if (locallab.qualitycurveMethod == "std") {
        lp.qualcurvemet = 1;
    } else if (locallab.qualitycurveMethod == "enh") {
        lp.qualcurvemet = 2;
    }

    float local_noiself = locallab.noiselumf;
    float local_noiselc = locallab.noiselumc;
    float local_noisecf = locallab.noisechrof;
    float local_noisecc = locallab.noisechroc;
    float multi[5];

    for (int y = 0; y < 5; y++) {
        multi[y] = ((float) locallab.mult[y]) / 100.f;
    }

    float thresho = ((float)locallab.threshold ) / 100.f;
    int local_chroma = locallab.chroma;
    int local_sensi = locallab.sensi;
    int local_sensibn = locallab.sensibn;
    int local_sensitm = locallab.sensitm;
    int local_sensih = locallab.sensih;
    int local_sensicb = locallab.sensicb;
    int local_contrast = locallab.contrast;
    float local_lightness = (float) locallab.lightness;
    int local_transit = locallab.transit;
    double radius = (double) locallab.radius;
    double sharradius = ((double) locallab.sharradius) / 100. ;
    int local_sensisha = locallab.sensisha;
    int local_sharamount = locallab.sharamount;
    int local_shardamping = locallab.shardamping;
    int local_shariter = locallab.shariter;
    bool inverse = locallab.invers;
    bool curvacti = locallab.curvactiv;
    bool acti = locallab.activlum;

    bool inverserad = locallab.inversrad;
    bool inverseret = locallab.inversret;
    bool inversesha = locallab.inverssha;
    double strength = (double) locallab.strength;
    float str = (float)locallab.str;
    lp.cir = circr;
    lp.actsp = acti;
    lp.xc = w * local_center_x;
    lp.yc = h * local_center_y;
    lp.lx = w * local_x;
    lp.ly = h * local_y;
    lp.lxL = w * local_xL;
    lp.lyT = h * local_yT;
    lp.chro = local_chroma;
    lp.sens = local_sensi;
    lp.sensh = local_sensih;
    lp.senscb = local_sensicb;
    lp.cont = local_contrast;
    lp.ligh = local_lightness;

    if (lp.ligh >= -2.f && lp.ligh <= 2.f) {
        lp.ligh /= 5.f;
    }

    lp.trans = local_transit;
    lp.rad = radius;
    lp.stren = strength;
    lp.sensbn = local_sensibn;
    lp.inv = inverse;
    lp.curvact = curvacti;
    lp.invrad = inverserad;
    lp.invret = inverseret;
    lp.invshar = inversesha;
    lp.str = str;
    lp.shrad = sharradius;
    lp.senssha = local_sensisha;
    lp.shamo = local_sharamount;
    lp.shdamp = local_shardamping;
    lp.shiter = local_shariter;
    lp.iterat = iterati;
    lp.dxx = w * local_dxx;
    lp.dyy = h * local_dyy;
    lp.thr = thre;
    lp.noiself = local_noiself;
    lp.noiself = local_noiself;
    lp.noiselc = local_noiselc;
    lp.noisecf = local_noisecf;
    lp.noisecc = local_noisecc;


    lp.strengt = streng;
    lp.gamm = gam;
    lp.esto = est;
    lp.scalt = scal_tm;
    lp.rewe = rewe;
    lp.senstm = local_sensitm;

    for (int y = 0; y < 5; y++) {
        lp.mulloc[y] = multi[y];
    }

    lp.threshol = thresho;
    lp.colorena = locallab.expcolor;
    lp.blurena = locallab.expblur;
    lp.tonemapena = locallab.exptonemap;
    lp.retiena = locallab.expreti;
    lp.sharpena = locallab.expsharp;
    lp.cbdlena = locallab.expcbdl;
    lp.denoiena = locallab.expdenoi;

}

static void calcTransition (const float lox, const float loy, const float ach, const local_params& lp, int &zone, float &localFactor)
{
    // returns the zone (0 = outside selection, 1 = transition zone between outside and inside selection, 2 = inside selection)
    // and a factor to calculate the transition in case zone == 1

    zone = 0;

    if (lox >= lp.xc && lox < (lp.xc + lp.lx) && loy >= lp.yc && loy < lp.yc + lp.ly) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lx)) + SQR ((loy - lp.yc) / (ach * lp.ly));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lx)) + SQR ((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor (lox, loy, lp.xc, lp.lx, lp.yc, lp.ly, ach);
            }
        }
    } else if (lox >= lp.xc && lox < lp.xc + lp.lx && loy < lp.yc && loy > lp.yc - lp.lyT) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lx)) + SQR ((loy - lp.yc) / (ach * lp.lyT));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lx)) + SQR ((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor (lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
            }
        }
    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy <= lp.yc && loy > lp.yc - lp.lyT) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lxL)) + SQR ((loy - lp.yc) / (ach * lp.lyT));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lxL)) + SQR ((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor (lox, loy, lp.xc, lp.lxL, lp.yc, lp.lyT, ach);
            }
        }
    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy > lp.yc && loy < lp.yc + lp.ly) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lxL)) + SQR ((loy - lp.yc) / (ach * lp.ly));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lxL)) + SQR ((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor (lox, loy, lp.xc, lp.lxL, lp.yc, lp.ly, ach);
            }
        }
    }
}

void ImProcFunctions::strcurv_data (std::string retistr, int *s_datc, int &siz)
{
    std::string delim[69] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
                             "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
                             "&", "#", "{", "[", "]", "}", "$", "*", "?", ">", "!", ";", "<", "(", ")", "+", "-"
                            };

    int s_size;
    std::size_t posend = retistr.find ("@");

    std::string strend = retistr.substr (posend - 1, 1);
    int longe = 0;

    for (int sl = 0; sl < 69; sl++) {
        if (delim[sl] == strend) {
            longe = sl + 1;
        }
    }

    s_size = longe;

    int s_datcu[s_size + 1];

    std::size_t pose[s_size + 1];
    pose[0] = -1;

    for (int z = 1; z < s_size + 1; z++) {
        pose[z] = retistr.find (delim[z - 1]);
    }


    for (int z = 1; z < s_size + 1; z++) {
        std::string sval = retistr.substr (pose[z - 1] + 1, (pose[z] - pose[z - 1]));
        s_datc[z - 1] = s_datcu[z - 1] = std::stoi (sval.c_str());

    }

    /*
    //here to verify process is good
        std::string cur_str = "";

        for(int j = 0; j < s_size; j++) {
            cur_str = cur_str + std::to_string(s_datcu[j]) + delim[j];
        }
        printf("calc str=%s\n", cur_str.c_str());
    */
    siz = longe;

}



void ImProcFunctions::addGaNoise (LabImage *lab, LabImage *dst, const float mean, const float variance, const int sk)
{
//   BENCHFUN
//Box-Muller method.
// add luma noise to image

    srand (1);

    const float variaFactor = SQR (variance) / sk;
    const float randFactor = 1.f / RAND_MAX;
#ifdef _OPENMP
    #pragma omp parallel
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

                float varia = SQR (kvar) * variaFactor;

                if (!generate) {
                    dst->L[y][x] = LIM (lab->L[y][x] + mean + varia * z1, 0.f, 32768.f);
                    continue;
                }

                int u1 = 0;
                int u2;

                while (u1 == 0) {
                    u1 = rand();
                    u2 = rand();
                }

                float u1f = u1 * randFactor;
                float u2f = u2 * randFactor;

                float2 sincosval = xsincosf (2.f * rtengine::RT_PI * u2f);
                float factor = sqrtf (-2.f * xlogf (u1f));
                z0 = factor * sincosval.y;
                z1 = factor * sincosval.x;

                dst->L[y][x] = LIM (lab->L[y][x] + mean + varia * z0, 0.f, 32768.f);

            }
        }
    }
}

void ImProcFunctions::DeNoise_Local (int call, const struct local_params& lp, LabImage* original, LabImage* transformed, const LabImage &tmp1, int cx, int cy)
{
    // local denoise
    // BENCHFUN
    const float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        const int loy = cy + y;
        const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

        if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
            for (int x = 0; x < transformed->W; x++) {
                transformed->L[y][x] = original->L[y][x];
                transformed->a[y][x] = original->a[y][x];
                transformed->b[y][x] = original->b[y][x];
            }

            continue;
        }

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);
            int begx = int (lp.xc - lp.lxL);
            int begy = int (lp.yc - lp.lyT);

            switch (zone) {
                case 0: { // outside selection and outside transition zone => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                    break;
                }

                case 1: { // inside transition zone
                    float factorx = localFactor;
                    float difL, difa, difb;

                    if (call == 2) { //simpleprocess
                        difL = tmp1.L[loy - begy][lox - begx] - original->L[y][x];
                        difa = tmp1.a[loy - begy][lox - begx] - original->a[y][x];
                        difb = tmp1.b[loy - begy][lox - begx] - original->b[y][x];
                    } else  { //dcrop
                        difL = tmp1.L[y][x] - original->L[y][x];
                        difa = tmp1.a[y][x] - original->a[y][x];
                        difb = tmp1.b[y][x] - original->b[y][x];

                    }

                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;
                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                    break;
                }

                case 2: { // inside selection => full effect, no transition
                    float difL, difa, difb;

                    if (call == 2) { //simpleprocess
                        difL = tmp1.L[loy - begy][lox - begx] - original->L[y][x];
                        difa = tmp1.a[loy - begy][lox - begx] - original->a[y][x];
                        difb = tmp1.b[loy - begy][lox - begx] - original->b[y][x];
                    } else  { //dcrop
                        difL = tmp1.L[y][x] - original->L[y][x];
                        difa = tmp1.a[y][x] - original->a[y][x];
                        difb = tmp1.b[y][x] - original->b[y][x];

                    }

                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                }
            }
        }
    }



}


void ImProcFunctions::cbdl_Local (int call, int sp, float ** buflight, float **loctemp, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
//local CBDL
    BENCHFUN
    const float ach = (float)lp.trans / 100.f;
    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.senscb - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.senscb;

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V (327.68f);
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
                vfloat av = LVFU (original->a[y][i]);
                vfloat bv = LVFU (original->b[y][i]);
                STVF (atan2Buffer[i], xatan2f (bv, av));
                STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
            }

#endif

            for (int x = 0; x < transformed->W; x++) {
                const int lox = cx + x;
                int begx = int (lp.xc - lp.lxL);
                int begy = int (lp.yc - lp.lyT);
                int zone = 0;

                float localFactor = 1.f;
                calcTransition (lox, loy, ach, lp, zone, localFactor);

                if (zone == 0) {
                    continue;
                }


#ifdef __SSE2__
                float rhue = atan2Buffer[x];
                float rchro = sqrtBuffer[x];
#else
                float rhue = xatan2f (original->b[y][x], original->a[y][x]);
                float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif
                //     int zone;

                //retrieve data
                float cli = 1.f;

                //    if (lp.curvact == true) {

                cli = (buflight[loy - begy][lox - begx]);
                //    }

                //parameters for linear interpolation in function of real hue
                float apluscligh = (1.f - cli) / delhu;
                float bpluscligh = 1.f - apluscligh * hueplus;
                float amoinscligh = (cli - 1.f) / delhu;
                float bmoinscligh = 1.f - amoinscligh * huemoins;

                float realcligh = 1.f;


                //   float localFactor = 1.f;
                //   calcTransition (lox, loy, ach, lp, zone, localFactor);
                //prepare shape detection
                float khu = 0.f;
                float kch = 1.f;
                float fach = 1.f;
                float deltachro = fabs (rchro - chromaref);
                float deltahue = fabs (rhue - hueref);

                if (deltahue > rtengine::RT_PI) {
                    deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                //kch to modulate action with chroma
                if (deltachro < 160.f * SQR (lp.senscb / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR (lp.senscb / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;
                }

                if (lp.senscb < 40.f ) {
                    kch = pow (kch, pa * lp.senscb + pb);   //increase under 40
                }


                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if (lp.senscb < 100.f) { //to try...
                    //hue detection
                    if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu )  {
                            realcligh = apluscligh * rhue + bpluscligh;

                            khu  = apl * rhue + bpl;
                        } else if (rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                            realcligh = amoinscligh * rhue + bmoinscligh;

                        } else {
                            khu = 1.f;
                            realcligh = cli;

                        }


//                            kzon = true;
                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realcligh = apluscligh * rhue + bpluscligh;

                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                            realcligh = amoinscligh * rhue + bmoinscligh;

                        } else {
                            khu = 1.f;
                            realcligh = cli;

                        }

//                            kzon = true;
                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins ) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realcligh = apluscligh * rhue + bpluscligh;

                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realcligh = amoinscligh * rhue + bmoinscligh;

                            khu = amo * rhue + bmo;
                        } else {
                            realcligh = cli;

                            khu = 1.f;

                        }

//                            kzon = true;
                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realcligh = apluscligh * rhue + bpluscligh;

                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                            realcligh = amoinscligh * rhue + bmoinscligh;

                        } else {
                            khu = 1.f;
                            realcligh = cli;

                        }

//                            kzon = true;
                    }

                    if (deltaE <  2.8f * lp.senscb) {
                        fach = khu;
                    } else {
                        fach = khu * (ahu * deltaE + bhu);
                    }


                    float kcr = 10.f;

                    if (rchro < kcr) {
                        fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                    }

                    if (lp.qualmet == 1) {
                    } else {
                        fach = 1.f;
                    }

                    //fach = khu ;

                } else {
                    /*
                        float kcr = 8.f;
                        if(lp.senssha > 30.f){
                        if (rchro < kcr) {
                            fach *= (1.f / (kcr)) * rchro;

                        }
                        }
                        */
                }

                float fli = ((100.f + realcligh) / 100.f);//luma transition
                float kcr = 100.f * lp.thr;
                float falL = 1.f;

                if (rchro < kcr && chromaref > kcr) { // reduce artifacts in grey tones near hue spot and improve algorithm
                    falL *= pow (rchro / kcr, lp.iterat / 10.f);
                }

                switch (zone) {
                    case 0: { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        break;
                    }

                    case 1: { // inside transition zone
                        float factorx = localFactor;
                        float difL = 0.f;

                        difL = loctemp[loy - begy][lox - begx] * fli * falL - original->L[y][x];

                        //float difL = loctemp[y][x] - original->L[y][x];
                        difL *= factorx;
                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;

                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        float difL = 0.f;
                        difL = loctemp[loy - begy][lox - begx] * fli * falL - original->L[y][x];

                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;
                    }
                }

                // }
            }
        }
    }
}


void ImProcFunctions::TM_Local (int call, int sp, LabImage * tmp1, float **buflight, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
//local TM
    BENCHFUN
    const float ach = (float)lp.trans / 100.f;
    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.senstm - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.senstm;

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V (327.68f);
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
                vfloat av = LVFU (original->a[y][i]);
                vfloat bv = LVFU (original->b[y][i]);
                STVF (atan2Buffer[i], xatan2f (bv, av));
                STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
            }

#endif

            for (int x = 0; x < transformed->W; x++) {
                const int lox = cx + x;
                const int begx = lp.xc - lp.lxL;
                const int begy = lp.yc - lp.lyT;

                float rL;

                if (lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && (rL = original->L[y][x]) > 3.2768f) {
                    // rL > 3.2768f to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                    int zone = 0;

                    float localFactor = 1.f;
                    calcTransition (lox, loy, ach, lp, zone, localFactor);

                    if (zone == 0) {
                        continue;
                    }

                    //     if (lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {

#ifdef __SSE2__
                    float rhue = atan2Buffer[x];
                    float rchro = sqrtBuffer[x];
#else
                    float rhue = xatan2f (original->b[y][x], original->a[y][x]);
                    float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif
                    // int zone;

                    //retrieve data
                    float cli = 1.f;

                    //    if (lp.curvact == true) {

                    cli = (buflight[loy - begy][lox - begx]);
                    //    }

                    //parameters for linear interpolation in function of real hue
                    float apluscligh = (1.f - cli) / delhu;
                    float bpluscligh = 1.f - apluscligh * hueplus;
                    float amoinscligh = (cli - 1.f) / delhu;
                    float bmoinscligh = 1.f - amoinscligh * huemoins;

                    float realcligh = 1.f;


                    //  float localFactor = 1.f;
                    //  calcTransition (lox, loy, ach, lp, zone, localFactor);
                    //prepare shape detection
                    float khu = 0.f;
                    float kch = 1.f;
                    float fach = 1.f;
                    float deltachro = fabs (rchro - chromaref);
                    float deltahue = fabs (rhue - hueref);

                    if (deltahue > rtengine::RT_PI) {
                        deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                    }

                    float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                    //kch to modulate action with chroma
                    if (deltachro < 160.f * SQR (lp.senstm / 100.f)) {
                        kch = 1.f;
                    } else {
                        float ck = 160.f * SQR (lp.senstm / 100.f);
                        float ak = 1.f / (ck - 160.f);
                        float bk = -160.f * ak;
                        kch = ak * deltachro + bk;
                    }

                    if (lp.senstm < 40.f ) {
                        kch = pow (kch, pa * lp.senstm + pb);   //increase under 40
                    }


                    // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                    if (lp.senstm < 100.f) { //to try...
                        //hue detection
                        if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                            if (rhue >= hueplus - delhu )  {
                                realcligh = apluscligh * rhue + bpluscligh;

                                khu  = apl * rhue + bpl;
                            } else if (rhue < huemoins + delhu)  {
                                realcligh = amoinscligh * rhue + bmoinscligh;

                                khu = amo * rhue + bmo;
                            } else {
                                khu = 1.f;
                                realcligh = cli;

                            }


//                            kzon = true;
                        } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                realcligh = apluscligh * rhue + bpluscligh;

                                khu  = apl * rhue + bpl;
                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                khu = amo * rhue + bmo;
                                realcligh = amoinscligh * rhue + bmoinscligh;

                            } else {
                                khu = 1.f;
                                realcligh = cli;

                            }

//                            kzon = true;
                        }

                        if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins ) {
                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                realcligh = apluscligh * rhue + bpluscligh;

                                khu  = apl * rhue + bpl;
                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                realcligh = amoinscligh * rhue + bmoinscligh;

                                khu = amo * rhue + bmo;
                            } else {
                                khu = 1.f;
                                realcligh = cli;

                            }

//                            kzon = true;
                        } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                realcligh = apluscligh * rhue + bpluscligh;

                                khu  = apl * rhue + bpl;
                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                realcligh = amoinscligh * rhue + bmoinscligh;

                                khu = amo * rhue + bmo;
                            } else {
                                khu = 1.f;
                                realcligh = cli;

                            }

//                            kzon = true;
                        }

                        if (deltaE <  2.8f * lp.senstm) {
                            fach = khu;
                        } else {
                            fach = khu * (ahu * deltaE + bhu);
                        }


                        float kcr = 10.f;

                        if (rchro < kcr) {
                            fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                        }

                        if (lp.qualmet == 1) {
                        } else {
                            fach = 1.f;
                        }

                        //fach = khu ;

                    } else {
                    }

                    float fli = ((100.f + realcligh) / 100.f);//luma transition
                    float kcr = 100.f * lp.thr;
                    float falL = 1.f;

                    if (rchro < kcr && chromaref > kcr) { // reduce artifacts in grey tones near hue spot and improve algorithm
                        falL *= pow (rchro / kcr, lp.iterat / 10.f);
                    }


                    switch (zone) {
                        case 0: { // outside selection and outside transition zone => no effect, keep original values
                            transformed->L[y][x] = original->L[y][x];

                            transformed->a[y][x] = original->a[y][x];
                            transformed->b[y][x] = original->b[y][x];

                            break;
                        }

                        case 1: { // inside transition zone
                            float factorx = localFactor;
                            float difL, difa, difb;

                            difL = tmp1->L[loy - begy][lox - begx] * fli * falL - original->L[y][x];
                            difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];

                            difL *= factorx;
                            difa *= factorx;
                            difb *= factorx;

                            transformed->L[y][x] = original->L[y][x] + difL * kch * fach;


                            transformed->a[y][x] = original->a[y][x] + difa * kch * fach;//same as Luma
                            transformed->b[y][x] = original->b[y][x] + difb * kch * fach;//same as Luma

                            break;
                        }

                        case 2: { // inside selection => full effect, no transition
                            float difL, difa, difb;

                            difL = tmp1->L[loy - begy][lox - begx] * fli * falL - original->L[y][x];
                            difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];

                            transformed->L[y][x] = original->L[y][x] + difL * kch * fach;


                            transformed->a[y][x] = original->a[y][x] + difa * kch * fach;//same as Luma
                            transformed->b[y][x] = original->b[y][x] + difb * kch * fach;//same as Luma
                        }
                    }

                }
            }
        }
    }
}




void ImProcFunctions::BlurNoise_Local (int call, int sp, LabImage * tmp1, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
//local BLUR
    BENCHFUN

    const float ach = (float)lp.trans / 100.f;
    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    constexpr float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    constexpr float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    constexpr float pb = 4.f;
    constexpr float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.sensbn - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.sensbn;

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V (327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < transformed->H; y++) {
            const int loy = cy + y;

            const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

            if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                for (int x = 0; x < transformed->W; x++) {
                    transformed->L[y][x] = original->L[y][x];
                }

                if (!lp.actsp) {
                    for (int x = 0; x < transformed->W; x++) {
                        transformed->a[y][x] = original->a[y][x];
                        transformed->b[y][x] = original->b[y][x];
                    }
                }

                continue;
            }

#ifdef __SSE2__
            int i = 0;

            for (; i < transformed->W - 3; i += 4) {
                vfloat av = LVFU (original->a[y][i]);
                vfloat bv = LVFU (original->b[y][i]);
                STVF (atan2Buffer[i], xatan2f (bv, av));
                STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
            }

#endif



            for (int x = 0, lox = cx + x; x < transformed->W; x++, lox++) {
                int zone = 0;
                float localFactor = 1.f;
                calcTransition (lox, loy, ach, lp, zone, localFactor);

                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];

                    if (!lp.actsp) {
                        transformed->a[y][x] = original->a[y][x];
                        transformed->b[y][x] = original->b[y][x];
                    }

                    continue;
                }

#ifdef __SSE2__
                const float rhue = atan2Buffer[x];
                const float rchro = sqrtBuffer[x];
#else
                const float rhue = xatan2f (original->b[y][x], original->a[y][x]);
                const float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif

                //prepare shape detection
                float kch = 1.f;
                float fach = 1.f;
                float deltachro = fabs (rchro - chromaref);
                float deltahue = fabs (rhue - hueref);

                if (deltahue > rtengine::RT_PI) {
                    deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                //kch to modulate action with chroma
                if (deltachro < 160.f * SQR (lp.sensbn / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR (lp.sensbn / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;
                }

                if (lp.sensbn < 40.f ) {
                    float khu = 0.f;
                    kch = pow (kch, pa * lp.sensbn + pb);   //increase under 40


                    // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                    if (lp.qualmet == 1 && lp.sensbn < 20.f) { //to try...
                        //hue detection
                        if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                            if (rhue >= hueplus - delhu )  {
                                khu  = apl * rhue + bpl;
                            } else if (rhue < huemoins + delhu)  {
                                khu = amo * rhue + bmo;
                            } else {
                                khu = 1.f;
                            }
                        } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                khu  = apl * rhue + bpl;
                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                khu = amo * rhue + bmo;
                            } else {
                                khu = 1.f;
                            }
                        }

                        if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins ) {
                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                khu  = apl * rhue + bpl;
                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                khu = amo * rhue + bmo;
                            } else {
                                khu = 1.f;
                            }
                        } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                khu  = apl * rhue + bpl;
                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                khu = amo * rhue + bmo;
                            } else {
                                khu = 1.f;
                            }
                        }

                        if (deltaE <  2.8f * lp.sensbn) {
                            fach = khu;
                        } else {
                            fach = khu * (ahu * deltaE + bhu);
                        }


                        float kcr = 10.f;

                        if (rchro < kcr) {
                            fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                        }
                    }
                }

                int begx = lp.xc - lp.lxL;
                int begy = lp.yc - lp.lyT;

                switch (zone) {

                    case 1: { // inside transition zone
                        float difL, difa, difb;

                        if (call == 2) {
                            difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                            difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                        } else {
                            difL = tmp1->L[y][x] - original->L[y][x];
                            difa = tmp1->a[y][x] - original->a[y][x];
                            difb = tmp1->b[y][x] - original->b[y][x];


                        }

                        difL *= localFactor;

                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;

                        if (!lp.actsp) {
                            difa *= localFactor;
                            difb *= localFactor;
                            transformed->a[y][x] = original->a[y][x] + difa * kch * fach;//same as Luma
                            transformed->b[y][x] = original->b[y][x] + difb * kch * fach;//same as Luma
                        }

                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        float difL, difa, difb;

                        if (call == 2) {
                            difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                            difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                        } else {
                            difL = tmp1->L[y][x] - original->L[y][x];
                            difa = tmp1->a[y][x] - original->a[y][x];
                            difb = tmp1->b[y][x] - original->b[y][x];

                        }

                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;

                        if (!lp.actsp) {
                            transformed->a[y][x] = original->a[y][x] + difa * kch * fach;//same as Luma
                            transformed->b[y][x] = original->b[y][x] + difb * kch * fach;//same as Luma
                        }
                    }
                }
            }
        }
    }
}

void ImProcFunctions::InverseReti_Local (const struct local_params & lp, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy, int chro)
{
    // BENCHFUN
//inverse local retinex
    float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);

            switch (zone) {
                case 0: { // outside selection and outside transition zone => full effect, no transition
                    if (chro == 0) {
                        transformed->L[y][x] = tmp1->L[y][x];
                    }

                    if (chro == 1) {

                        transformed->a[y][x] = tmp1->a[y][x];
                        transformed->b[y][x] = tmp1->b[y][x];
                    }

                    break;
                }

                case 1: { // inside transition zone
                    float factorx = 1.f - localFactor;

                    if (chro == 0) {
                        float difL = tmp1->L[y][x] - original->L[y][x];
                        difL *= factorx;
                        transformed->L[y][x] = original->L[y][x] + difL;
                    }

                    if (chro == 1) {
                        float difa = tmp1->a[y][x] - original->a[y][x];
                        float difb = tmp1->b[y][x] - original->b[y][x];

                        difa *= factorx;
                        difb *= factorx;

                        transformed->a[y][x] = original->a[y][x] + difa;
                        transformed->b[y][x] = original->b[y][x] + difb;
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



void ImProcFunctions::Reti_Local (int call, float **buflight, float **bufchro, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const struct local_params & lp, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy, int chro)
{

//local retinex
    BENCHFUN {
        const float ach = (float)lp.trans / 100.f;

        //chroma
        constexpr float amplchsens = 2.5f;
        constexpr float achsens = (amplchsens - 1.f) / (100.f - 20.f); //20. default locallab.sensih
        constexpr float bchsens = 1.f - 20.f * achsens;
        const float multchro = lp.sensh * achsens + bchsens;

        //luma

        //skin
        constexpr float amplchsensskin = 1.6f;
        constexpr float achsensskin = (amplchsensskin - 1.f) / (100.f - 20.f); //20. default locallab.sensih
        constexpr float bchsensskin = 1.f - 20.f * achsensskin;
        const float multchroskin = lp.sensh * achsensskin + bchsensskin;

        //transition = difficult to avoid artifact with scope on flat area (sky...)

        constexpr float delhu = 0.1f; //between 0.05 and 0.2

        const float apl = (-1.f) / delhu;
        const float bpl = - apl * hueplus;
        const float amo = 1.f / delhu;
        const float bmo = - amo * huemoins;


        const float pb = 4.f;
        const float pa = (1.f - pb) / 40.f;

        const float ahu = 1.f / (2.8f * lp.sensh - 280.f);
        const float bhu = 1.f - ahu * 2.8f * lp.sensh;

        const float alum = 1.f / (lp.sensh - 100.f);
        const float blum = 1.f - alum * lp.sensh;


#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
#ifdef __SSE2__
            float atan2Buffer[transformed->W] ALIGNED16;
            float sqrtBuffer[transformed->W] ALIGNED16;
            vfloat c327d68v = F2V (327.68f);
#endif

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int y = 0; y < transformed->H; y++)
            {

                const int loy = cy + y;
                const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

                if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                    for (int x = 0; x < transformed->W; x++) {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    continue;
                }

#ifdef __SSE2__
                int i = 0;

                for (; i < transformed->W - 3; i += 4) {
                    vfloat av = LVFU (original->a[y][i]);
                    vfloat bv = LVFU (original->b[y][i]);
                    STVF (atan2Buffer[i], xatan2f (bv, av));
                    STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                    sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
                }

#endif

                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int begx = int (lp.xc - lp.lxL);
                    int begy = int (lp.yc - lp.lyT);

                    int zone = 0;
                    float localFactor = 1.f;
                    calcTransition (lox, loy, ach, lp, zone, localFactor);

                    if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        continue;
                    }

#ifdef __SSE2__
                    float rhue = atan2Buffer[x];
                    float rchro = sqrtBuffer[x];
#else
                    float rhue = xatan2f (original->b[y][x], original->a[y][x]);
                    float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif
                    float rL = original->L[y][x] / 327.68f;

                    float cli = 1.f;
                    float clc = 1.f;

                    //                       if (lp.curvact == true) {
                    cli = (buflight[loy - begy][lox - begx]);
                    clc = (bufchro[loy - begy][lox - begx]);

                    //                      } else {
                    //                          cli = lp.str;
                    //                           clc = params->locallab.chrrt;
                    //                       }

                    float aplus = (1.f - cli) / delhu;
                    float bplus = 1.f - aplus * hueplus;
                    float amoins = (cli - 1.f) / delhu;
                    float bmoins = 1.f - amoins * huemoins;

                    float aplusch = (1.f - clc) / delhu;
                    float bplusch = 1.f - aplusch * hueplus;
                    float amoinsch = (clc - 1.f) / delhu;
                    float bmoinsch = 1.f - amoinsch * huemoins;

                    float realstr = 1.f;
                    float realstrch = 1.f;
                    //prepare shape detection
                    float deltachro = fabs (rchro - chromaref);
                    float deltahue = fabs (rhue - hueref);

                    if (deltahue > rtengine::RT_PI) {
                        deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                    }

                    float deltaE = 20.f * deltahue + deltachro; //between 0 and 280
                    float deltaL = fabs (lumaref - rL); //between 0 and 100

                    float kch = 1.f;
                    float khu = 0.f;
                    float fach = 1.f;
                    float falu = 1.f;

                    if (deltachro < 160.f * SQR (lp.sensh / 100.f)) {
                        kch = 1.f;
                    } else {
                        float ck = 160.f * SQR (lp.sensh / 100.f);
                        float ak = 1.f / (ck - 160.f);
                        float bk = -160.f * ak;
                        kch = ak * deltachro + bk;
                    }

                    if (lp.sensh < 40.f ) {
                        kch = pow (kch, pa * lp.sensh + pb);   //increase under 40
                    }

                    bool kzon = false;

                    //transition = difficult to avoid artifact with scope on flat area (sky...)
                    //hue detection
                    if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;

                        } else if (rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = cli;
                            khu = 1.f;
                            realstrch = clc;

                        }

                        kzon = true;
                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = cli;
                            khu = 1.f;
                            realstrch = clc;

                        }

                        kzon = true;
                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = cli;
                            khu = 1.f;
                            realstrch = clc;

                        }

                        kzon = true;
                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = cli;
                            khu = 1.f;
                            realstrch = clc;

                        }

                        kzon = true;
                    }

                    //shape detection for hue chroma and luma
                    if (lp.sensh <= 20.f) { //to try...

                        if (deltaE <  2.8f * lp.sensh) {
                            fach = khu;
                        } else {
                            fach = khu * (ahu * deltaE + bhu);
                        }

                        float kcr = 10.f;

                        if (rchro < kcr) {
                            fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                        }

                        if (lp.qualmet >= 1) {
                        } else {
                            fach = 1.f;
                        }

                        if (deltaL <  lp.sensh) {
                            falu = 1.f;
                        } else {
                            falu = alum * deltaL + blum;
                        }

                    }


                    //    float kdiff = 0.f;
                    // I add these functions...perhaps not good
                    if (kzon) {
                        if (lp.sensh < 60.f) { //arbitrary value
                            if (hueref < -1.1f && hueref > -2.8f) { // detect blue sky
                                if (chromaref > 0.f && chromaref < 35.f * multchro) { // detect blue sky
                                    if ( (rhue > -2.79f && rhue < -1.11f) && (rchro < 35.f * multchro)) {
                                        realstr *= 0.9f;
                                    } else {
                                        realstr = 1.f;
                                    }
                                }
                            } else {
                                realstr = cli;
                            }

                            if (lp.sensh < 50.f) { //&& lp.chro > 0.f
                                if (hueref > -0.1f && hueref < 1.6f) { // detect skin
                                    if (chromaref > 0.f && chromaref < 55.f * multchroskin) { // detect skin
                                        if ( (rhue > -0.09f && rhue < 1.59f) && (rchro < 55.f * multchroskin)) {
                                            realstr *= 0.7f;
                                        } else {
                                            realstr = 1.f;
                                        }
                                    }
                                } else {
                                    realstr = cli;
                                }
                            }
                        }

                    }

                    float kcr = 100.f * lp.thr;
                    float falL = 1.f;

                    if (rchro < kcr && chromaref > kcr) { // reduce artifacts in grey tones near hue spot and improve algorithm
                        falL *= pow (rchro / kcr, lp.iterat / 10.f);
                    }

                    //                     int zone;
                    //                     float localFactor;
                    //                     calcTransition (lox, loy, ach, lp, zone, localFactor);

                    if (rL > 0.1f) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                        switch (zone) {
                            case 0: { // outside selection and outside transition zone => no effect, keep original values
                                if (chro == 0) {
                                    transformed->L[y][x] = original->L[y][x];
                                }

                                if (chro == 1) {
                                    transformed->a[y][x] = original->a[y][x];
                                    transformed->b[y][x] = original->b[y][x];
                                }

                                break;
                            }

                            case 1: { // inside transition zone
                                float factorx = localFactor;

                                if (chro == 0) {
                                    float difL;

                                    difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                                    difL *= factorx * (100.f + realstr * falL) / 100.f;
                                    difL *= kch * fach;

                                    transformed->L[y][x] = original->L[y][x] + difL;
                                }

                                if (chro == 1) {
                                    float difa, difb;

                                    difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                                    difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                                    difa *= factorx * (100.f + realstrch * falu * falL) / 100.f;
                                    difb *= factorx * (100.f + realstrch * falu * falL) / 100.f;
                                    transformed->a[y][x] = CLIPC (original->a[y][x] + difa);
                                    transformed->b[y][x] = CLIPC (original->b[y][x] + difb);

                                }

                                break;

                            }

                            case 2: { // inside selection => full effect, no transition
                                if (chro == 0) {
                                    float difL;

                                    difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                                    difL *= (100.f + realstr * falL) / 100.f;
                                    difL *= kch * fach;
                                    transformed->L[y][x] = original->L[y][x] + difL;

                                }

                                if (chro == 1) {
                                    float difa, difb;

                                    difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                                    difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                                    difa *= (100.f + realstrch * falu * falL) / 100.f;
                                    difb *= (100.f + realstrch * falu * falL) / 100.f;
                                    transformed->a[y][x] = CLIPC (original->a[y][x] + difa);
                                    transformed->b[y][x] = CLIPC (original->b[y][x] + difb);

                                }
                            }
                        }

                        //}
                    }
                }
            }
        }

    }
}

void ImProcFunctions::InverseBlurNoise_Local (const struct local_params & lp, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy)
{
    // BENCHFUN
//inverse local blur and noise
    float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);

            switch (zone) {
                case 0: { // outside selection and outside transition zone => full effect, no transition
                    transformed->L[y][x] = tmp1->L[y][x];

                    if (!lp.actsp) {
                        transformed->a[y][x] = tmp1->a[y][x];
                        transformed->b[y][x] = tmp1->b[y][x];
                    }

                    break;
                }

                case 1: { // inside transition zone
                    float difL = tmp1->L[y][x] - original->L[y][x];
                    float difa = tmp1->a[y][x] - original->a[y][x];
                    float difb = tmp1->b[y][x] - original->b[y][x];

                    float factorx = 1.f - localFactor;
                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;

                    transformed->L[y][x] = original->L[y][x] + difL;

                    if (!lp.actsp) {

                        transformed->a[y][x] = original->a[y][x] + difa;
                        transformed->b[y][x] = original->b[y][x] + difb;
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

struct local_contra {
    float alsup, blsup;
    float alsup2, blsup2;
    float alsup3, blsup3;
    float alinf;
    float aDY;
    float aa;
    float bb;
    float aaa, bbb;
    float ccc;
    float dx, dy;
    float ah, bh;
    float al, bl;
};

void ImProcFunctions::Contrast_Local (int call, float ave, LabImage * bufcontorig, float ** buflightc, float moy, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, float pm, struct local_contra & lco, float lumaref, float av, const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
    BENCHFUN
// contrast - perhaps for 4 areas   if need
// I tried shmap adaptaed to Lab, but no real gain and artifacts
    const float localtype = lumaref; // always spot area
    const float ach = (float)lp.trans / 100.f;
    float reducac;

    //constant and variable to prepare shape detection
    if (lp.sens < 30.f) {
        reducac = 0.2f * (lp.sens / 100.f);
    } else {
        float areduc = 0.6285714f; //0.44f/0.7f;
        float breduc = 0.5f - areduc;
        reducac = areduc * (lp.sens / 100.f) + breduc;
    }

    const float realcox = lco.dx, realcoy = lco.dy;

    lco.alsup = (-realcox) / (localtype / 2.f);
    lco.blsup = -lco.alsup * localtype;
    lco.alsup2 = (realcoy) / (50.f - localtype / 2.f);
    lco.blsup2 = -lco.alsup2 * localtype;
    lco.alsup3 = (realcoy) / (localtype / 2.f - 50.f);
    lco.blsup3 = -lco.alsup3 * 100.f;
    lco.aDY = realcoy;


    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.sens - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.sens;

    lco.alinf = realcox / (localtype / 2.f);
    const float vi = (localtype / 2.f) / 100.f;
    const float vinf = (50.f + localtype / 2.f) / 100.f;
    ImProcFunctions::secondeg_begin (reducac, vi, lco.aa, lco.bb);//parabolic
    ImProcFunctions::secondeg_end (reducac, vinf, lco.aaa, lco.bbb, lco.ccc);//parabolic

    if (call <= 3) {
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            //Todo optimization in this first part with something equivalent to bufcolorig and bufcoltra in colorlight_local
#ifdef __SSE2__
            float atan2Buffer[transformed->W] ALIGNED16;
            float sqrtBuffer[transformed->W] ALIGNED16;
            vfloat c327d68v = F2V (327.68f);
#endif

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int y = 0; y < transformed->H; y++)
            {
                const int loy = cy + y;
                const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

                if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                    continue;
                }

#ifdef __SSE2__
                int i = 0;

                for (; i < transformed->W - 3; i += 4) {
                    vfloat av = LVFU (original->a[y][i]);
                    vfloat bv = LVFU (original->b[y][i]);
                    STVF (atan2Buffer[i], xatan2f (bv, av));
                    STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                    sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
                }

#endif

                for (int x = 0; x < transformed->W; x++) {
                    const int lox = cx + x;

                    float rL;

                    if (lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && (rL = original->L[y][x]) > 3.2768f) {
                        // rL > 3.2768f to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                        int zone = 0;

                        float localFactor = 1.f;
                        calcTransition (lox, loy, ach, lp, zone, localFactor);

                        if (zone == 0) {
                            continue;
                        }

#ifdef __SSE2__
                        float rhue = atan2Buffer[x];
                        float rchro = sqrtBuffer[x];
#else
                        float rhue = xatan2f (original->b[y][x], original->a[y][x]);
                        float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif
                        //prepare shape detection
                        float khu = 0.f;
                        float kch = 1.f;
                        float fach = 1.f;

                        float cli = 1.f;

                        const int begx = lp.xc - lp.lxL;
                        const int begy = lp.yc - lp.lyT;

                        if (lp.curvact) {

                            cli = buflightc[loy - begy][lox - begx];

                            if (cli ==  0.0f) {
                                cli = 0.01f;
                            }
                        }

                        //parameters for linear interpolation in function of real hue
                        float apluscligh = (1.f - cli) / delhu;
                        float bpluscligh = 1.f - apluscligh * hueplus;
                        float amoinscligh = (cli - 1.f) / delhu;
                        float bmoinscligh = 1.f - amoinscligh * huemoins;
                        float realcligh = 1.f;
                        float deltachro = fabs (rchro - chromaref);
                        float deltahue = fabs (rhue - hueref);

                        if (deltahue > rtengine::RT_PI) {
                            deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                        }

                        float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                        //kch to modulate action with chroma
                        if (deltachro < 160.f * SQR (lp.sens / 100.f)) { // TODOPRECOMPUTE
                            kch = 1.f;
                        } else {
                            float ck = 160.f * SQR (lp.sens / 100.f);
                            float ak = 1.f / (ck - 160.f);
                            float bk = -160.f * ak;
                            kch = ak * deltachro + bk;

                            if (lp.sens < 40.f ) {
                                kch = pow_F (kch, pa * lp.sens + pb);   //increase under 40
                            }
                        }

                        // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                        if (lp.sens < 100.f && lp.qualmet >= 1) { //to try...
                            //hue detection
                            if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                                if (rhue >= hueplus - delhu )  {
                                    realcligh = apluscligh * rhue + bpluscligh;

                                    khu  = apl * rhue + bpl;
                                } else if (rhue < huemoins + delhu)  {
                                    realcligh = amoinscligh * rhue + bmoinscligh;

                                    khu = amo * rhue + bmo;
                                } else {
                                    realcligh = cli;

                                    khu = 1.f;
                                }


                            } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                                if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                    realcligh = apluscligh * rhue + bpluscligh;

                                    khu  = apl * rhue + bpl;
                                } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                    realcligh = amoinscligh * rhue + bmoinscligh;

                                    khu = amo * rhue + bmo;
                                } else {
                                    realcligh = cli;

                                    khu = 1.f;
                                }

                            }

                            if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins ) {
                                if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                    realcligh = apluscligh * rhue + bpluscligh;

                                    khu  = apl * rhue + bpl;
                                } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                    realcligh = amoinscligh * rhue + bmoinscligh;
                                    khu = amo * rhue + bmo;
                                } else {
                                    realcligh = cli;

                                    khu = 1.f;
                                }

                            } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                                if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                    realcligh = apluscligh * rhue + bpluscligh;

                                    khu  = apl * rhue + bpl;
                                } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                    realcligh = amoinscligh * rhue + bmoinscligh;

                                    khu = amo * rhue + bmo;
                                } else {
                                    realcligh = cli;

                                    khu = 1.f;
                                }

                            }

                            if (deltaE <  2.8f * lp.sens) {
                                fach = khu;
                            } else {
                                fach = khu * (ahu * deltaE + bhu);
                            }

                            constexpr float kcr = 10.f;

                            if (rchro < kcr) {
                                fach *= SQR (rchro) / SQR (kcr);
                            }
                        }

                        float kcr = 100.f * lp.thr;
                        float falL = 1.f;

                        if (rchro < kcr && chromaref > kcr) { // reduce artifacts in grey tones near hue spot and improve algorithm
                            falL *= pow_F (rchro / kcr, lp.iterat / 10.f);
                        }

                        float modu =  1.f ;//realclig / cli;
                        float localty = localtype;

                        switch (zone) {

                            case 1: { // inside transition zone
                                if (!lp.curvact) {
                                    modu = 1.f;
                                } else {
                                    modu = CLIP1 (realcligh / cli);

                                }


                                if (original->L[y][x] < 32768.f) {
                                    float factorx = localFactor;
                                    float prov100 = original->L[y][x] / 32768.f;
                                    float prov = prov100 * 100.f;

                                    if (prov > localty) {
                                        if (prov >= localty && prov < 50.f + localty / 2.f) {
                                            float core = (lco.alsup2 * prov + lco.blsup2) ;
                                            core *= factorx;

                                            transformed->L[y][x] = 327.68f * (prov + pm * (prov - localty) * (core) * kch * fach * falL * modu);
                                        } else {
                                            float core = lco.aDY * (lco.aaa * prov100 * prov100 + lco.bbb * prov100 + lco.ccc);

                                            core *= factorx;

                                            transformed->L[y][x] = 327.68f * (prov + pm * (prov - localty) * (core) * kch * fach * falL * modu);
                                        }
                                    } else  { //inferior
                                        if (2.f * prov > localty && prov < localty)  {
                                            float core = (lco.alsup * prov + lco.blsup) ;
                                            core *= factorx;

                                            transformed->L[y][x] = 327.68f * (prov - pm * (localty - prov) * core * kch * fach * falL * modu);
                                        } else if (2.f * prov <= localtype) {
                                            float core = prov * lco.alinf * (lco.aa * prov100 * prov100 + lco.bb * prov100);

                                            core *= factorx;

                                            transformed->L[y][x] = 327.68f * (prov - pm * (localty - prov) * core * kch * fach * falL * modu);
                                        }
                                    }
                                }

                                break;
                            }

                            case 2: { // inside selection => full effect, no transition
                                if (!lp.curvact) {
                                    modu = 1.f;
                                } else {
                                    //   modu = realcligh / (cli + 0.001f);
                                    //  printf("mo=%f", modu);
                                    modu = CLIP1 (realcligh / cli);

                                }

                                if (original->L[y][x] < 32768.f) {
                                    float prov100 = original->L[y][x] / 32768.f;
                                    float prov = prov100 * 100.f;

                                    if (prov > localty  ) {
                                        if (prov >= localty && prov < 50.f + localty / 2.f) {
                                            float core = (lco.alsup2 * prov + lco.blsup2) ;
                                            transformed->L[y][x] = 327.68f * (prov + pm * (prov - localty) * core * kch * fach * falL * modu);
                                        } else {
                                            float core = lco.aDY * (lco.aaa * prov100 * prov100 + lco.bbb * prov100 + lco.ccc);
                                            transformed->L[y][x] = 327.68f * (prov + pm * (prov - localty) * core * kch * fach * falL * modu);
                                        }
                                    } else  { //inferior
                                        if (2.f * prov > localty && prov < localty)  {
                                            float core = (lco.alsup * prov + lco.blsup) ;
                                            transformed->L[y][x] = 327.68f * (prov - pm * (localty - prov) * core * kch * fach * falL * modu);
                                        } else if (2.f * prov <= localtype) {
                                            float core = prov * lco.alinf * (lco.aa * prov100 * prov100 + lco.bb * prov100);
                                            transformed->L[y][x] = 327.68f * (prov - pm * (localty - prov) * core * kch * fach * falL * modu);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void ImProcFunctions::InverseContrast_Local (float ave, const local_contra & lco, const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
    //  BENCHFUN
    float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);

            switch (zone) {
                case 0: { // outside selection and outside transition zone => full effect, no transition
                    if (original->L[y][x] < 32768.f) {
                        float prov = original->L[y][x];

                        if (original->L[y][x] > ave) {
                            float kh = lco.ah * (original->L[y][x] / 327.68f) + lco.bh;
                            original->L[y][x] = ave + kh * (original->L[y][x] - ave);
                        } else  {
                            float kl = lco.al * (original->L[y][x] / 327.68f) + 1.f;
                            original->L[y][x] = ave - kl * (ave - original->L[y][x]);
                        }

                        float diflc = original->L[y][x] - prov;
                        transformed->L[y][x] =  prov + diflc;
                    } else {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    break;
                }

                case 1: { // inside transition zone
                    if (original->L[y][x] < 32768.f) {
                        float factorx = localFactor;
                        factorx = 1.f - factorx;
                        float prov = original->L[y][x];

                        if (original->L[y][x] > ave) {
                            float kh = lco.ah * (original->L[y][x] / 327.68f) + lco.bh;
                            original->L[y][x] = ave + kh * (original->L[y][x] - ave);
                        } else  {
                            float kl = lco.al * (original->L[y][x] / 327.68f) + 1.f;
                            original->L[y][x] = ave - kl * (ave - original->L[y][x]);
                        }

                        float diflc = original->L[y][x] - prov;
                        diflc *= factorx;
                        transformed->L[y][x] =  prov + diflc;

                    } else {
                        transformed->L[y][x] = original->L[y][x];
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

static void calclight (float lum, float  koef, float & lumnew, bool inv)
//replace L-curve that does not work in local or bad
{

    float blac = 0.3f;

    if (inv == false) {
        blac = 0.99f;
    } else {
        if (koef < -90.f) {
            blac = -0.069f * koef - 5.91f;
        }
    }

    if (koef >= 0.f) {
        lumnew = lum + 0.2f * (33000.f - lum) * koef / 100.f;
    }

    if (koef < 0.f) {
        lumnew = lum + blac * lum * koef / 100.f;//0.999 instead of 0.2

        if (lumnew < 0.f) {
            float kc = lum / (lum - lumnew);
            lumnew = lum + kc * 0.2f * lum * koef / 100.f;

        }

        //    if (inv == false && koef == -100.f) {
        if (koef == -100.f) {
            lumnew = 0.f;
        }

    }

    lumnew = CLIPLOC (lumnew);

}

void ImProcFunctions::InverseSharp_Local (int sp, float **loctemp, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
//local sharp
    //  BENCHFUN
    const float ach = (float)lp.trans / 100.f;
    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.senssha - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.senssha;

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V (327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < transformed->H; y++) {
#ifdef __SSE2__
            int i = 0;

            for (; i < transformed->W - 3; i += 4) {
                vfloat av = LVFU (original->a[y][i]);
                vfloat bv = LVFU (original->b[y][i]);
                STVF (atan2Buffer[i], xatan2f (bv, av));
                STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
            }

#endif

            int loy = cy + y;

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;
#ifdef __SSE2__
                float rhue = atan2Buffer[x];
                float rchro = sqrtBuffer[x];
#else
                float rhue = xatan2f (original->b[y][x], original->a[y][x]);
                float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif
                int zone;
                float localFactor = 1.f;
                calcTransition (lox, loy, ach, lp, zone, localFactor);
                //prepare shape detection
                float khu = 0.f;
                float kch = 1.f;
                float fach = 1.f;
                float deltachro = fabs (rchro - chromaref);
                float deltahue = fabs (rhue - hueref);

                if (deltahue > rtengine::RT_PI) {
                    deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                //kch to modulate action with chroma
                if (deltachro < 160.f * SQR (lp.senssha / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR (lp.senssha / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;
                }

                if (lp.senssha < 40.f ) {
                    kch = pow (kch, pa * lp.senssha + pb);   //increase under 40
                }


                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if (lp.senssha < 20.f) { //to try...
                    //hue detection
                    if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu )  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins ) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    }

                    if (deltaE <  2.8f * lp.senssha) {
                        fach = khu;
                    } else {
                        fach = khu * (ahu * deltaE + bhu);
                    }


                    float kcr = 10.f;

                    if (rchro < kcr) {
                        fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                    }

                    if (lp.qualmet >= 1) {
                    } else {
                        fach = 1.f;
                    }

                    //fach = khu ;

                } else {
                    /*
                        float kcr = 8.f;
                        if(lp.senssha > 30.f){
                        if (rchro < kcr) {
                            fach *= (1.f / (kcr)) * rchro;

                        }
                        }
                        */
                }



                switch (zone) {
                    case 0: { // outside selection and outside transition zone => full effect, no transition
                        float difL = loctemp[y][x] - original->L[y][x];
                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;

                        break;
                    }

                    case 1: { // inside transition zone
                        float difL = loctemp[y][x] - original->L[y][x];

                        float factorx = 1.f - localFactor;
                        difL *= factorx;

                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;
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


void ImProcFunctions::Sharp_Local (int call, int sp, float **loctemp, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
    BENCHFUN
    const float ach = (float)lp.trans / 100.f;
    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.senssha - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.senssha;

    const bool detectHue = lp.senssha < 20.f && lp.qualmet >= 1;
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V (327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < transformed->H; y++) {

            const int loy = cy + y;
            const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

            if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                for (int x = 0; x < transformed->W; x++) {
                    transformed->L[y][x] = original->L[y][x];
                }

                continue;
            }

#ifdef __SSE2__
            int i = 0;

            if (detectHue) {
                for (; i < transformed->W - 3; i += 4) {
                    vfloat av = LVFU (original->a[y][i]);
                    vfloat bv = LVFU (original->b[y][i]);
                    STVF (atan2Buffer[i], xatan2f (bv, av));
                    STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                    sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
                }
            } else {
                for (; i < transformed->W - 3; i += 4) {
                    vfloat av = LVFU (original->a[y][i]);
                    vfloat bv = LVFU (original->b[y][i]);
                    STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
                }
            }

#endif

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;
                int zone = 0;
                float localFactor = 1.f;
                calcTransition (lox, loy, ach, lp, zone, localFactor);

                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    continue;
                }

#ifdef __SSE2__
                float rchro = sqrtBuffer[x];
#else
                float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif
                //prepare shape detection
                float kch = 1.f;
                float fach = 1.f;
                float deltachro = fabs (rchro - chromaref);

                //kch to modulate action with chroma
                if (deltachro < 160.f * SQR (lp.senssha / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR (lp.senssha / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;

                    if (lp.senssha < 40.f ) {
                        kch = pow_F (kch, pa * lp.senssha + pb);   //increase under 40
                    }
                }

                if (lp.senssha >= 99.f) {
                    kch = 1.f;
                }

                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if (detectHue) { //to try...
#ifdef __SSE2__
                    float rhue = atan2Buffer[x];
#else
                    float rhue = xatan2f (original->b[y][x], original->a[y][x]);
#endif
                    float khu = 0.f;
                    float deltahue = fabs (rhue - hueref);

                    if (deltahue > rtengine::RT_PI) {
                        deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                    }

                    //hue detection
                    if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu )  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }


                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins ) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    }

                    float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                    if (deltaE <  2.8f * lp.senssha) {
                        fach = khu;
                    } else {
                        fach = khu * (ahu * deltaE + bhu);
                    }


                    float kcr = 10.f;

                    if (rchro < kcr) {
                        fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                    }

                }

                int begx = int (lp.xc - lp.lxL);
                int begy = int (lp.yc - lp.lyT);

                switch (zone) {

                    case 1: { // inside transition zone
                        float factorx = localFactor;
                        float difL;

                        if (call == 2) {
                            difL = loctemp[loy - begy][lox - begx] - original->L[y][x];
                        } else {
                            difL = loctemp[y][x] - original->L[y][x];

                        }

                        difL *= factorx;
                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;

                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        float difL;

                        if (call == 2) {
                            difL = loctemp[loy - begy][lox - begx] - original->L[y][x];
                        } else  {
                            difL = loctemp[y][x] - original->L[y][x];
                        }

                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;
                    }
                }
            }
        }
    }
}



void ImProcFunctions::ColorLight_Local (int call, LabImage * bufcolorig, float ** buflight, float ** bufchro, float ** buflightslid, int sp, float moy, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, bool locallutili, LUTf & lllocalcurve, const LocLHCurve & loclhCurve, LUTf & cclocalcurve, float chprov, float cligh, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
    BENCHFUN
// chroma and lightness
    const float ach = (float)lp.trans / 100.f;

    //chroma
    constexpr float amplchsens = 2.5f;
    constexpr float achsens = (amplchsens - 1.f) / (100.f - 20.f); //20. default locallab.sensi
    constexpr float bchsens = 1.f - 20.f * achsens;
    const float multchro = lp.sens * achsens + bchsens;

    //luma
    constexpr float ampllumsens = 2.f;
    constexpr float alumsens = (ampllumsens - 1.f) / (100.f - 20.f); //20. default locallab.sensi
    constexpr float blumsens = 1.f - 20.f * alumsens;
    const float multlum = lp.sens * alumsens + blumsens;

    //skin
    constexpr float amplchsensskin = 1.6f;
    constexpr float achsensskin = (amplchsensskin - 1.f) / (100.f - 20.f); //20. default locallab.sensi
    constexpr float bchsensskin = 1.f - 20.f * achsensskin;
    const float multchroskin = lp.sens * achsensskin + bchsensskin;

    //transition = difficult to avoid artifact with scope on flat area (sky...)
    constexpr float delhu = 0.1f; //between 0.05 and 0.2 ==> minima for scope
    //constexpr float delhu2 = 0.03f; //between 0.05 and 0.2

    const float aplus = (1.f - lp.chro) / delhu;
    const float bplus = 1.f - aplus * hueplus;
    const float amoins = (lp.chro - 1.f) / delhu;
    const float bmoins = 1.f - amoins * huemoins;

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.sens - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.sens;

    //luma
    constexpr float lumdelta = 11.f; //11
    float modlum = lumdelta * multlum;

    // constant and variables to prepare shape detection
    if (lumaref + modlum >= 100.f) {
        modlum = (100.f - lumaref) / 2.f;
    }

    if (lumaref - modlum <= 0.f) {
        modlum = (lumaref) / 2.f;
    }

    float aa, bb, aaa, bbb, ccc;
    float reducac = settings->reduchigh;//0.85f;
    float reducac2 = settings->reduclow;//0.2f;

    float vinf = (lumaref + modlum) / 100.f;
    float vi = (lumaref - modlum) / 100.f;
    ImProcFunctions::secondeg_begin (reducac, vi, aa, bb);//parabolic
    ImProcFunctions::secondeg_end (reducac, vinf, aaa, bbb, ccc);//parabolic
    float vinf2 = (lumaref + modlum) / 100.f;
    float vi2 = (lumaref - modlum) / 100.f;
    float aaaa, bbbb, cccc, aO, bO;
    ImProcFunctions::secondeg_end (reducac2, vinf2, aaaa, bbbb, cccc);//parabolic
    ImProcFunctions::secondeg_begin (reducac2, vi2, aO, bO);//parabolic

    if (call <= 3) {
        //Todo optimization in this first part with bufcolorig and bufcoltra

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
#ifdef __SSE2__
            float atan2Buffer[transformed->W] ALIGNED16;
            float sqrtBuffer[transformed->W] ALIGNED16;
            vfloat c327d68v = F2V (327.68f);
#endif

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int y = 0; y < transformed->H; y++)
            {
                const int loy = cy + y;
                const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

                if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                    continue;
                }

#ifdef __SSE2__
                int i = 0;

                for (; i < transformed->W - 3; i += 4) {
                    vfloat av = LVFU (original->a[y][i]);
                    vfloat bv = LVFU (original->b[y][i]);
                    STVF (atan2Buffer[i], xatan2f (bv, av));
                    STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                    sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
                }

#endif


                for (int x = 0; x < transformed->W; x++) {
                    const int lox = cx + x;
                    const int begx = int (lp.xc - lp.lxL);
                    const int begy = int (lp.yc - lp.lyT);

                    int zone = 0;

                    float localFactor = 1.f;
                    calcTransition (lox, loy, ach, lp, zone, localFactor);

                    if (zone == 0) {
                        continue;
                    }


#ifdef __SSE2__
                    float rhue = atan2Buffer[x];
                    float rchro = sqrtBuffer[x];
#else

                    float rhue = xatan2f (original->b[y][x], original->a[y][x]);

                    float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif

                    float rL = original->L[y][x] / 327.68f;
                    float rLL = original->L[y][x] / 327.68f;

                    if (fabs (original->b[y][x]) < 0.01f) {
                        original->b[y][x] = 0.01f;
                    }

                    float eps = 0.f;

                    if (fabs (original->b[y][x]) < 0.001f) {
                        eps = 0.01f;
                    }

                    //retriev data curve lightness
                    float cli = (buflight[loy - begy][lox - begx]);
                    //parameters for linear interpolation in function of real hue
                    float apluscligh = (1.f - cli) / delhu;
                    float bpluscligh = 1.f - apluscligh * hueplus;
                    float amoinscligh = (cli - 1.f) / delhu;
                    float bmoinscligh = 1.f - amoinscligh * huemoins;

                    float cchro = (bufchro[loy - begy][lox - begx]);
                    float apluscurv = (1.f - cchro) / delhu;
                    float bpluscurv = 1.f - apluscurv * hueplus;
                    float amoinscurv = (cchro - 1.f) / delhu;
                    float bmoinscurv = 1.f - amoinscurv * huemoins;

                    float clisl = (buflightslid[loy - begy][lox - begx]);
                    //parameters for linear interpolation in function of real hue
                    float aplusclighsl = (1.f - clisl) / delhu;
                    float bplusclighsl = 1.f - aplusclighsl * hueplus;
                    float amoinsclighsl = (clisl - 1.f) / delhu;
                    float bmoinsclighsl = 1.f - amoinsclighsl * huemoins;

                    float kab = (original->a[y][x] / (original->b[y][x] + eps));

                    //prepare shape detection
                    // real... = coefficient to apply at lightness, chroma,...
                    float realchro = 1.f;
                    float realcurv = 1.f;
                    float realcligh = 1.f;
                    float realclighsl = 1.f;

                    //evaluate delta Hue and delta Chro
                    float deltachro = fabs (rchro - chromaref);

                    float deltahue = fabs (rhue - hueref);

                    if (deltahue > rtengine::RT_PI) {
                        deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                    }

                    //pseudo deltaE
                    float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280
                    float deltaL = fabs (lumaref - rL); //between 0 and 100

                    float kch = 1.f;
                    float khu = 0.f;
                    float fach = 1.f;
                    float falu = 1.f;

                    //kch acts on luma
                    if (deltachro < 160.f * SQR (lp.sens / 100.f)) {
                        kch = 1.f;
                    } else {
                        float ck = 160.f * SQR (lp.sens / 100.f);
                        float ak = 1.f / (ck - 160.f);
                        float bk = -160.f * ak;
                        kch = ak * deltachro + bk;
                    }

                    if (lp.sens < 40.f ) {
                        kch = pow (kch, pa * lp.sens + pb);   //increase under 40
                    }

                    bool kzon = false;

                    //transition = difficult to avoid artifact with scope on flat area (sky...)
                    //hue detection
                    //for each quart calculate realchro, realcligh,... in function of Hue pixel
                    if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu)  {
                            realchro = aplus * rhue + bplus;
                            realcurv = apluscurv * rhue + bpluscurv;
                            realcligh = apluscligh * rhue + bpluscligh;
                            realclighsl = aplusclighsl * rhue + bplusclighsl;
                            khu  = apl * rhue + bpl;

                        } else if (rhue < huemoins + delhu)  {
                            realchro = amoins * rhue + bmoins;
                            realcurv = amoinscurv * rhue + bmoinscurv;
                            realcligh = amoinscligh * rhue + bmoinscligh;
                            realclighsl = amoinsclighsl * rhue + bmoinsclighsl;

                            khu = amo * rhue + bmo;

                        } else {
                            realchro = lp.chro;
                            realcurv = cchro;
                            realcligh = cli;
                            realclighsl = clisl;

                            khu = 1.f;

                        }

                        kzon = true;
                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realchro = aplus * rhue + bplus;
                            realcurv = apluscurv * rhue + bpluscurv;
                            realcligh = apluscligh * rhue + bpluscligh;
                            realclighsl = aplusclighsl * rhue + bplusclighsl;

                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realchro = amoins * rhue + bmoins;
                            realcurv = amoinscurv * rhue + bmoinscurv;
                            realcligh = amoinscligh * rhue + bmoinscligh;
                            realclighsl = amoinsclighsl * rhue + bmoinsclighsl;

                            khu = amo * rhue + bmo;

                        } else {
                            realchro = lp.chro;

                            realcurv = cchro;
                            realcligh = cli;
                            realclighsl = clisl;

                            khu = 1.f;

                        }

                        kzon = true;
                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realchro = aplus * rhue + bplus;
                            realcurv = apluscurv * rhue + bpluscurv;
                            realcligh = apluscligh * rhue + bpluscligh;
                            realclighsl = aplusclighsl * rhue + bplusclighsl;

                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realchro = amoins * rhue + bmoins;
                            realcurv = amoinscurv * rhue + bmoinscurv;
                            realcligh = amoinscligh * rhue + bmoinscligh;
                            realclighsl = amoinsclighsl * rhue + bmoinsclighsl;

                            khu = amo * rhue + bmo;

                        } else {
                            realchro = lp.chro;

                            realcurv = cchro;
                            realcligh = cli;
                            realclighsl = clisl;

                            khu = 1.f;

                        }

                        kzon = true;
                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realchro = aplus * rhue + bplus;
                            realcurv = apluscurv * rhue + bpluscurv;
                            realcligh = apluscligh * rhue + bpluscligh;
                            realclighsl = aplusclighsl * rhue + bplusclighsl;

                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realchro = amoins * rhue + bmoins;
                            realcurv = amoinscurv * rhue + bmoinscurv;
                            realcligh = amoinscligh * rhue + bmoinscligh;
                            realclighsl = amoinsclighsl * rhue + bmoinsclighsl;

                            khu = amo * rhue + bmo;

                        } else {
                            realchro = lp.chro;

                            realcurv = cchro;
                            realcligh = cli;
                            realclighsl = clisl;

                            khu = 1.f;

                        }

                        kzon = true;
                    }


                    //detection of deltaE and deltaL
                    if (lp.sens <= 20.f) { //to try...
                        //fach and kch acts on luma
                        if (deltaE <  2.8f * lp.sens) {
                            fach = khu;
                        } else {
                            fach = khu * (ahu * deltaE + bhu);
                        }

                        float kcr = 10.f;

                        if (rchro < kcr) {
                            fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                        }

                        //fach = 1.f;//to avoid artifacts in some cases
                        //can be probably improved
                        if (lp.qualmet >= 1) {
                        } else {
                            fach = 1.f;
                        }

                        //falu acts on chroma
                        if (deltaL <  lp.sens) {
                            falu = 1.f;
                        } else {
                            falu = 1.f;// alum * deltaL + blum;
                        }

                    }

                    if (kzon) {
                        if (lp.sens < 60.f) { //arbitrary value
                            if (hueref < -1.1f && hueref > -2.8f) { // detect blue sky
                                if (chromaref > 0.f && chromaref < 35.f * multchro) { // detect blue sky
                                    if ( (rhue > -2.79f && rhue < -1.11f) && (rchro < 35.f * multchro)) {
                                        realchro *= 0.9f;
                                        realcurv *= 0.9f;
                                    } else {
                                        realchro = 1.f;
                                        realcurv = 1.f;

                                    }
                                }
                            } else {
                                realchro = lp.chro;
                                realcurv = cchro;

                            }

                            if (lp.sens < 50.f && lp.chro > 0.f) {
                                if (hueref > -0.1f && hueref < 1.6f) { // detect skin
                                    if (chromaref > 0.f && chromaref < 55.f * multchroskin) { // detect skin
                                        if ( (rhue > -0.09f && rhue < 1.59f) && (rchro < 55.f * multchroskin)) {
                                            realchro *= 0.9f;
                                            realcurv *= 0.9f;

                                        } else {
                                            realchro = 1.f;
                                            realcurv = 1.f;

                                        }
                                    }
                                } else {
                                    realchro = lp.chro;
                                    realcurv = cchro;

                                }
                            }
                        }

                    }

                    float kLinf = rLL / (100.f);
                    float kLsup = kLinf;

                    float kdiff = 1.f;

                    if (kzon) { ///rhue < hueplus && rhue > huemoins

                        if ( (rLL > (lumaref - modlum) && rLL < (lumaref + modlum))) {
                            kdiff = 1.f;
                        } else if (rLL > 0.f && rLL <= (lumaref - modlum)) {
                            kdiff = (aa * kLinf * kLinf + bb * kLinf);   //parabolic

                            if (kdiff < 0.01f) {
                                kdiff = 0.01f;
                            }
                        } else if (rLL <= 100.f && rLL >= (lumaref + modlum)) {

                            kdiff = (aaa * kLsup * kLsup + bbb * kLsup + ccc);   //parabolic

                            if (kdiff < 0.01f) {
                                kdiff = 0.01f;
                            }

                        }

                        //end luma
                    } else {
                        float ktes = 1.f;

                        if ( (rLL > (lumaref - modlum) && rLL < (lumaref + modlum))) {
                            kdiff = ktes;
                        } else if (rLL > 0.f && rLL <= (lumaref - modlum)) {

                            kdiff = (ktes * (aO * kLinf * kLinf + bO * kLinf));    //parabolic

                            if (kdiff < 0.01f) {
                                kdiff = 0.01f;
                            }

                        } else if (rLL <= 100.f && rLL >= (lumaref + modlum)) {

                            kdiff = (ktes * (aaaa * kLsup * kLsup + bbbb * kLsup + cccc));    //parabolic

                            if (kdiff < 0.01f) {
                                kdiff = 0.01f;
                            }

                        }

                    }

                    float kcr = 100.f * lp.thr;
                    float falL = 1.f;

                    if (rchro < kcr && chromaref > kcr) { // reduce artifacts in grey tones near hue spot and improve algorithm
                        falL *= pow (rchro / kcr, lp.iterat / 10.f);
                    }


                    //     int zone;
                    //     float localFactor;
                    //     calcTransition (lox, loy, ach, lp, zone, localFactor);
                    float th_r = 0.01f;

                    if (rL > th_r) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                        switch (zone) {
                            case 0: { // outside selection and outside transition zone => no effect, keep original values
                                transformed->L[y][x] = original->L[y][x];
                                transformed->a[y][x] = original->a[y][x];
                                transformed->b[y][x] = original->b[y][x];
                                break;
                            }

                            case 1: { // inside transition zone
                                float lumnew = bufcolorig->L[loy - begy][lox - begx];

                                float lightcont;

                                if (lp.qualcurvemet == 1) {

                                    if (lllocalcurve) {
                                        float lumprov = lllocalcurve[lumnew * 1.9f];
                                        float lumred = 0.526316f * lumprov; //0.526316f
                                        lumnew = lumnew + (lumred - lumnew) / 4.f;//reduce sensibility

                                    }

                                    if (loclhCurve) {
                                        float l_r;//Luminance Lab in 0..1
                                        l_r = lumnew / 32768.f;
                                        {
                                            float khu = 1.9f; //in reserve in case of!

                                            float valparam = float ((loclhCurve[500.f * Color::huelab_to_huehsv2 (rhue)] - 0.5f)); //get l_r=f(H)
                                            float valparamneg;
                                            valparamneg = valparam;

                                            if (valparam > 0.f) {
                                                l_r = (1.f - valparam) * l_r + valparam * (1.f - SQR (((SQR (1.f - min (l_r, 1.0f))))));
                                            } else
                                                //for negative
                                            {
                                                l_r *= (1.f + khu * valparamneg);
                                            }
                                        }

                                        lumnew = l_r * 32768.f;
                                    }

                                }

                                if (lp.ligh != 0.f && lp.curvact == false) {
                                    calclight (lumnew, lp.ligh , lumnew, true);//replace L-curve
                                    lightcont = lumnew;

                                } else {
                                    lightcont = lumnew;
                                }

                                float factorx = localFactor;
                                float fli = 1.f;
                                float flisl = 1.f;

                                if (lp.curvact && lp.ligh != 0.f) {
                                    flisl = ((100.f + realclighsl * falL ) / 100.f);//luma transition
                                }

                                if (lp.qualcurvemet == 2) {
                                    fli = ((100.f + realcligh * falL ) / 100.f);//luma transition
                                }

                                float flicur = 1.f;

                                if (lp.qualcurvemet != 0) {
                                    flicur = ((100.f + realcurv * factorx * falu * falL) / 100.f);
                                }

                                float fac = flicur *  (100.f + factorx * realchro * falu * falL) / 100.f; //chroma factor transition
                                //if(fac < 0.2f) fac = 0.2f;
                                float diflc = lightcont * fli * flisl - original->L[y][x];
                                kdiff *= fach * kch;
                                diflc *= kdiff ;

                                diflc *= factorx; //transition lightness
                                transformed->L[y][x] = CLIPL (1.f * (original->L[y][x] + diflc));


                                if (fabs (kab) > 1.f) {
                                    transformed->a[y][x] = CLIPC (original->a[y][x] * fac) ;
                                    transformed->b[y][x] = CLIPC (original->a[y][x] * fac) / kab;
                                } else {
                                    transformed->b[y][x] = CLIPC (original->b[y][x] * fac);
                                    transformed->a[y][x] = CLIPC (original->b[y][x] * fac) * kab ;

                                }

                                break;
                            }

                            case 2: { // inside selection => full effect, no transition
                                float lumnew = bufcolorig->L[loy - begy][lox - begx];
                                float lightcont;

                                if (lp.qualcurvemet == 1) {

                                    if (lllocalcurve) {
                                        float lumprov = lllocalcurve[lumnew * 1.9f];
                                        float lumred = 0.526316 * lumprov; // 0.526316f
                                        lumnew = lumnew + (lumred - lumnew) / 4.f;//reduce sensibility
                                    }

                                    if (loclhCurve) {
                                        float l_r;//Luminance Lab in 0..1
                                        l_r = lumnew / 32768.f;
                                        {
                                            float khu = 1.9f;

                                            float valparam = float ((loclhCurve[500.f * Color::huelab_to_huehsv2 (rhue)] - 0.5f)); //get l_r=f(H)
                                            float valparamneg;
                                            valparamneg = valparam;

                                            if (valparam > 0.f) {
                                                l_r = (1.f - valparam) * l_r + valparam * (1.f - SQR (((SQR (1.f - min (l_r, 1.0f))))));
                                            } else
                                                //for negative
                                            {
                                                l_r *= (1.f + khu * valparamneg);
                                            }
                                        }

                                        lumnew = l_r * 32768.f;
                                    }

                                }


                                if (lp.ligh != 0.f && lp.curvact == false) {
                                    calclight (lumnew, lp.ligh , lumnew, true);//replace L-curve
                                    lightcont = lumnew;

                                } else {
                                    lightcont = lumnew;
                                }

                                float fli = 1.f;
                                float flisl = 1.f;

                                if (lp.curvact && lp.ligh != 0.f) {
                                    flisl = ((100.f + realclighsl * falL ) / 100.f);//luma transition
                                }

                                if (lp.qualcurvemet == 2) {
                                    fli = ((100.f + realcligh * falL) / 100.f);//luma transition
                                }

                                float flicur = 1.f;

                                if (lp.qualcurvemet != 0) {
                                    flicur = ((100.f + realcurv * falu * falL) / 100.f);
                                }

                                float fac = flicur * (100.f + realchro * falu * falL) / 100.f; //chroma factor transition7
                                //if(fac < 0.2f) fac = 0.2f;

                                float diflc = lightcont * fli * flisl - original->L[y][x];

                                kdiff *= fach * kch;
                                diflc *= kdiff ;
                                transformed->L[y][x] = CLIPL (1.f * (original->L[y][x] + diflc));

                                if (fabs (kab) > 1.f) {
                                    transformed->a[y][x] = CLIPC (original->a[y][x] * fac) ;
                                    transformed->b[y][x] = CLIPC (original->a[y][x] * fac) / kab;
                                } else {
                                    transformed->b[y][x] = CLIPC (original->b[y][x] * fac);
                                    transformed->a[y][x] = CLIPC (original->b[y][x] * fac) * kab;
                                }

                            }

                        }
                    }

                    //    }
                }
            }
        }


    }

}

void ImProcFunctions::InverseColorLight_Local (const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
    // BENCHFUN
    float ach = (float)lp.trans / 100.f;
    const float facc = (100.f + lp.chro) / 100.f; //chroma factor transition

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);

            switch (zone) {
                case 0: { // outside selection and outside transition zone => no effect, keep original values
                    float lumnew = original->L[y][x];

                    if (lp.ligh != 0.f) {
                        calclight (original->L[y][x], lp.ligh , lumnew, false);
                    }

                    float lightcont = lumnew ; //original->L[y][x] + (lp.ligh /100.f)*original->L[y][x] ; //apply lightness



                    transformed->L[y][x] = lightcont; //localcurve[original->L[y][x]];  //apply lightness
                    transformed->a[y][x] = original->a[y][x] * facc;
                    transformed->b[y][x] = original->b[y][x] * facc;
                    break;
                }

                case 1: { // inside transition zone
                    float factorx = 1.f - localFactor;
                    float fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition
                    float lumnew = original->L[y][x];

                    if (lp.ligh != 0.f) {
                        calclight (original->L[y][x], lp.ligh , lumnew, false);
                    }

                    float lightcont = lumnew ; //apply lightness

                    float diflc = lightcont - original->L[y][x];
                    diflc *= factorx;
                    transformed->L[y][x] = original->L[y][x] + diflc;
                    transformed->a[y][x] = original->a[y][x] * fac;
                    transformed->b[y][x] = original->b[y][x] * fac;
                    break;
                }

                case 2: { // inside selection => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                }
            }
        }
    }

}
void ImProcFunctions::calc_ref (int call, int sp, float** shbuffer, LabImage * original, LabImage * transformed, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh, bool locutili, int sk, const LocretigainCurve & locRETgainCcurve, bool locallutili, LUTf & lllocalcurve, const LocLHCurve & loclhCurve, LUTf & cclocalcurve, double & hueref, double & chromaref, double & lumaref)
{
    if (params->locallab.enabled) {
        //always calculate hueref, chromaref, lumaref  before others operations use in normal mode for all modules exceprt denoise
        struct local_params lp;
        calcLocalParams (oW, oH, params->locallab, lp);

// double precision for large summations
        double aveA = 0.;
        double aveB = 0.;
        double aveL = 0.;
        double aveChro = 0.;
// int precision for the counters
        int nab = 0;
// single precision for the result
        float avA, avB, avL;
        int spotSize = 0.88623f * max (1,  lp.cir / sk); //18
        //O.88623 = sqrt(PI / 4) ==> sqare equal to circle

        // very small region, don't use omp here
        for (int y = max (cy, (int) (lp.yc - spotSize)); y < min (transformed->H + cy, (int) (lp.yc + spotSize + 1)); y++) {
            for (int x = max (cx, (int) (lp.xc - spotSize)); x < min (transformed->W + cx, (int) (lp.xc + spotSize + 1)); x++) {
                aveL += original->L[y - cy][x - cx];
                aveA += original->a[y - cy][x - cx];
                aveB += original->b[y - cy][x - cx];
                aveChro += sqrtf (SQR (original->b[y - cy][x - cx]) + SQR (original->a[y - cy][x - cx]));

                nab++;
            }
        }

        aveL = aveL / nab;
        aveA = aveA / nab;
        aveB = aveB / nab;
        aveChro = aveChro / nab;
        aveChro /= 327.68f;
        avA = aveA / 327.68f;
        avB = aveB / 327.68f;
        avL = aveL / 327.68f;
        hueref = xatan2f (avB, avA);   //mean hue
        chromaref = aveChro;
        lumaref = avL;
    }
}

void ImProcFunctions::Lab_Local (int call, int sp, float** shbuffer, LabImage * original, LabImage * transformed, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh, bool locutili, int sk, const LocretigainCurve & locRETgainCcurve, bool locallutili, LUTf & lllocalcurve, const LocLHCurve & loclhCurve, LUTf & cclocalcurve, double & hueref, double & chromaref, double & lumaref)
{
    //general call of others functions : important return hueref, chromaref, lumaref
    if (params->locallab.enabled) {
        BENCHFUN
#ifdef _DEBUG
        MyTime t1e, t2e;
        t1e.set();
// init variables to display Munsell corrections
        MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif

        int del = 3; // to avoid crash with [loy - begy] and [lox - begx] and bfh bfw  // with gtk2 [loy - begy-1] [lox - begx -1 ] and del = 1
        float moy = 0.f;

        struct local_params lp;
        calcLocalParams (oW, oH, params->locallab, lp);

        const float radius = lp.rad / (sk * 1.4f); //0 to 70 ==> see skip

        double ave = 0.;
        int n = 0;
        float av = 0;
        int levred;
        bool noiscfactiv = false;

        if (lp.qualmet == 2) { //suppress artifacts with quality enhanced
            levred = 4;
            noiscfactiv = true;
        }    else {
            levred = 7;
            noiscfactiv = false;
        }

        if (lp.inv || lp.invret) { //exterior || lp.curvact
            ave = 0.f;
            n = 0;
            #pragma omp parallel for reduction(+:ave,n)

            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if (lox >= lp.xc && lox < lp.xc + lp.lx && loy >= lp.yc && loy < lp.yc + lp.ly) {
                    } else if (lox >= lp.xc && lox < lp.xc + lp.lx && loy < lp.yc && loy > lp.yc - lp.lyT) {
                    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy <= lp.yc && loy > lp.yc - lp.lyT) {
                    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy > lp.yc && loy < lp.yc + lp.ly) {
                    } else {
                        ave += original->L[y][x];
                        n++;
                    }
                }
            }

            if (n == 0) {
                ave = 15000.f;
                n = 1;
            }

            ave = ave / n;
            av = ave / 327.68f;
        }

        //printf ("call= %i sp=%i hueref=%f chromaref=%f lumaref=%f\n", call, sp, hueref, chromaref, lumaref);
        struct local_contra lco;

// we must here detect : general case, skin, sky,...foliages ???
// delta dhue, luminance and chroma
        constexpr float ared = (rtengine::RT_PI - 0.05f) / 100.f;

        constexpr float bred = 0.05f;

        float dhue = ared * lp.sens + bred; //delta hue lght chroma

        float dhueret = ared * lp.sensh + bred; //delta hue retinex

        constexpr float maxh = 3.5f; // 3.5 amplification contrast above mean

        constexpr float maxl = 2.5f; // 3 reductio contrast under mean

        float multh = (float) fabs (lp.cont) * (maxh - 1.f) / 100.f + 1.f;

        float mult = (float)fabs (lp.cont) * (maxl - 1.f) / 100.f + 1.f;

        lco.dx = 1.f - 1.f / mult;

        lco.dy = 1.f - 1.f / multh;


//Blur and noise

        if (((radius >= 1.5 * GAUSS_SKIP && lp.rad > 1.) || lp.stren > 0.1)  && lp.blurena) { // radius < GAUSS_SKIP means no gauss, just copy of original image
            LabImage *tmp1 = nullptr;
            LabImage *bufgb = nullptr;
            int GW = transformed->W;
            int GH = transformed->H;
            //  printf ("rad=%f gaus=%f call=%i skip=%i\n", radius, GAUSS_SKIP, call, sk);

            if (call == 2  && !lp.invrad) { //simpleprocess
                int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                int bfw = int (lp.lx + lp.lxL) + del;
                bufgb = new LabImage (bfw, bfh);

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufgb->L[ir][jr] = 0.f;
                        bufgb->a[ir][jr] = 0.f;
                        bufgb->b[ir][jr] = 0.f;
                    }

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufgb->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufgb->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufgb->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas
                        }
                    }

                /*
                There is a bug in this calculation ==> out of limits ==> crash
                                int yStart = lp.yc - lp.lyT - cy;
                                int yEnd = lp.yc + lp.ly - cy;
                                int xStart = lp.xc - lp.lxL - cx;
                                int xEnd = lp.xc + lp.lx - cx;
                                int begy = lp.yc - lp.lyT;
                                int begx = lp.xc - lp.lxL;
                #ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16)
                #endif

                                for (int y = yStart; y < yEnd ; y++) {
                                    int loy = cy + y;

                                    for (int x = xStart, lox = cx + x; x < xEnd; x++, lox++) {
                                        bufgb->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                                        bufgb->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                                        bufgb->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas
                                    }
                                }
                */
                tmp1 = new LabImage (bfw, bfh);
#ifdef _OPENMP
                #pragma omp parallel
#endif
                {
                    gaussianBlur (bufgb->L, tmp1->L, bfw, bfh, radius);
                    gaussianBlur (bufgb->a, tmp1->a, bfw, bfh, radius);
                    gaussianBlur (bufgb->b, tmp1->b, bfw, bfh, radius);

                }

            } else {
                tmp1 = new LabImage (transformed->W, transformed->H);;

#ifdef _OPENMP
                #pragma omp parallel
#endif
                {
                    gaussianBlur (original->L, tmp1->L, GW, GH, radius);
                    gaussianBlur (original->a, tmp1->a, GW, GH, radius);
                    gaussianBlur (original->b, tmp1->b, GW, GH, radius);

                }
            }

            if (lp.stren > 0.1f) {
                float mean = 0.f;//0 best result
                float variance = lp.stren ; //(double) SQR(lp.stren)/sk;
                addGaNoise (tmp1, tmp1, mean, variance, sk) ;
            }

            if (!lp.invrad) { //blur and noise (center)
                //         BlurNoise_Local(call, lp, original, transformed, tmp1, cx, cy);
                float hueplus = hueref + dhue;
                float huemoins = hueref - dhue;

                if (hueplus > rtengine::RT_PI) {
                    hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
                }

                if (huemoins < -rtengine::RT_PI) {
                    huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
                }

                BlurNoise_Local (call, sp, tmp1, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, original, transformed, cx, cy);

            } else {

                InverseBlurNoise_Local (lp, original, transformed, tmp1, cx, cy);

            }

            if (call == 2  && !lp.invrad) {
                delete bufgb;
            }

            delete tmp1;
        }

        //      }

//local denoise
        //all these variables are to prevent use of denoise when non necessary
        // but with qualmet = 2 (default for best quality) we must denoise chroma with little values to prevent artifacts due to variations of Hue
        // but if user select volontary denoise, it is that choice the good (prioritary)
        bool execdenoi = false ;
        bool execcolor = (lp.chro != 0.f || lp.ligh != 0.f || lp.cont != 0.f); // only if one slider ore more is engaged
        bool execbdl = (lp.mulloc[0] != 1.f || lp.mulloc[1] != 1.f || lp.mulloc[2] != 1.f || lp.mulloc[3] != 1.f || lp.mulloc[4] != 1.f) ;//only if user want cbdl
        execdenoi = noiscfactiv && ((lp.colorena && execcolor) || (lp.tonemapena && lp.strengt != 0.f) || (lp.cbdlena && execbdl) || (lp.sharpena && lp.shrad > 0.42) || (lp.retiena  && lp.str > 0.f));

        if (((lp.noiself > 0.f || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f) && lp.denoiena) || execdenoi) {
            StopWatch Stop1 ("locallab Denoise called");

            if (lp.noisecf > 0.1f || lp.noisecc > 0.1f) {
                noiscfactiv = false;
                levred = 7;
            }

#ifdef _OPENMP
            const int numThreads = omp_get_max_threads();
#else
            const int numThreads = 1;

#endif

            if (call == 1) {
                LabImage tmp1 (transformed->W, transformed->H);
                int GW = transformed->W;
                int GH = transformed->H;

                for (int ir = 0; ir < GH; ir++)
                    for (int jr = 0; jr < GW; jr++) {
                        tmp1.L[ir][jr] = original->L[ir][jr];
                        tmp1.a[ir][jr] = original->a[ir][jr];
                        tmp1.b[ir][jr] = original->b[ir][jr];
                    }

                int DaubLen = 6;

                int levwavL = levred;
                int skip = 1;

                wavelet_decomposition Ldecomp (tmp1.L[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, DaubLen);
                wavelet_decomposition adecomp (tmp1.a[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, DaubLen);
                wavelet_decomposition bdecomp (tmp1.b[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, DaubLen);

                float madL[8][3];
                int edge = 2;

                if (!Ldecomp.memoryAllocationFailed) {
                    #pragma omp parallel for collapse(2) schedule(dynamic,1)

                    for (int lvl = 0; lvl < levred; lvl++) {
                        for (int dir = 1; dir < 4; dir++) {
                            int Wlvl_L = Ldecomp.level_W (lvl);
                            int Hlvl_L = Ldecomp.level_H (lvl);

                            float ** WavCoeffs_L = Ldecomp.level_coeffs (lvl);

                            madL[lvl][dir - 1] = SQR (Mad (WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                        }
                    }

                    float vari[levred];

                    if (levred == 7) {
                        edge = 2;
                        vari[0] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[1] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));

                        vari[3] = 8.f * SQR ((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[4] = 8.f * SQR ((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[5] = 8.f * SQR ((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[6] = 8.f * SQR ((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    } else if (levred == 4) {
                        edge = 3;
                        vari[0] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[1] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[3] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiselc / 25.0));

                    }

                    if (( lp.noiself > 0.1f ||  lp.noiselc > 0.1f)) {
                        vari[0] = max (0.0001f, vari[0]);
                        vari[1] = max (0.0001f, vari[1]);
                        vari[2] = max (0.0001f, vari[2]);
                        vari[3] = max (0.0001f, vari[3]);

                        if (levred == 7) {
                            vari[4] = max (0.0001f, vari[4]);
                            vari[5] = max (0.0001f, vari[5]);
                            vari[6] = max (0.0001f, vari[6]);
                        }

                        float* noisevarlum = nullptr;  // we need a dummy to pass it to WaveletDenoiseAllL

                        WaveletDenoiseAllL (Ldecomp, noisevarlum, madL, vari, edge, numThreads);
                    }
                }

                float variC[levred];

                if (!adecomp.memoryAllocationFailed && !bdecomp.memoryAllocationFailed) {
                    if (levred == 7) {
                        edge = 2;
                        variC[0] = SQR (lp.noisecf / 10.0);
                        variC[1] = SQR (lp.noisecf / 10.0);
                        variC[2] = SQR (lp.noisecf / 10.0);

                        variC[3] = SQR (lp.noisecf / 10.0);
                        variC[4] = SQR (lp.noisecf / 10.0);
                        variC[5] = SQR (lp.noisecc / 10.0);
                        variC[6] = SQR (lp.noisecc / 10.0);
                    } else if (levred == 4) {
                        edge = 3;
                        variC[0] = SQR (lp.noisecf / 10.0);
                        variC[1] = SQR (lp.noisecf / 10.0);
                        variC[2] = SQR (lp.noisecf / 10.0);
                        variC[3] = SQR (lp.noisecf / 10.0);
                    }


                    if (( lp.noisecf > 0.1f ||  lp.noisecc > 0.1f  || noiscfactiv)) {
                        float minic = 0.0001f;

                        if (noiscfactiv) {
                            minic = 0.01f;//only for artifact shape detection
                        }

                        variC[0] = max (minic, variC[0]);
                        variC[1] = max (minic, variC[1]);
                        variC[2] = max (minic, variC[2]);
                        variC[3] = max (minic, variC[3]);

                        if (levred == 7) {

                            variC[4] = max (0.0001f, variC[4]);
                            variC[5] = max (0.0001f, variC[5]);
                            variC[6] = max (0.0001f, variC[6]);
                        }

                        float* noisevarchrom = new float[GH * GW];

                        for (int q = 0; q < GH * GW; q++) {
                            noisevarchrom[q] = 1.f;
                        }

                        float noisevarab_r = 100.f; //SQR(lp.noisecc / 10.0);
                        WaveletDenoiseAllAB (Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, false, false, false, numThreads);
                        WaveletDenoiseAllAB (Ldecomp, bdecomp, noisevarchrom, madL, variC, edge, noisevarab_r, false, false, false, numThreads);
                        delete[] noisevarchrom;

                    }
                }

                if (!Ldecomp.memoryAllocationFailed) {

                    Ldecomp.reconstruct (tmp1.L[0]);
                }

                if (!adecomp.memoryAllocationFailed) {

                    adecomp.reconstruct (tmp1.a[0]);
                }

                if (!bdecomp.memoryAllocationFailed) {

                    bdecomp.reconstruct (tmp1.b[0]);
                }

                DeNoise_Local (call, lp, original, transformed, tmp1, cx, cy);

            } else if (call == 2) { //simpleprocess

                int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                int bfw = int (lp.lx + lp.lxL) + del;
                LabImage bufwv (bfw, bfh);
                bufwv.clear (true);

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufwv.L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufwv.a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufwv.b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas
                        }

                    }

                int DaubLen = 6;

                int levwavL = levred;
                int skip = 1;
                wavelet_decomposition Ldecomp (bufwv.L[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, DaubLen);
                wavelet_decomposition adecomp (bufwv.a[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, DaubLen);
                wavelet_decomposition bdecomp (bufwv.b[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, DaubLen);

                float madL[8][3];
                int edge = 2;

                if (!Ldecomp.memoryAllocationFailed) {
                    #pragma omp parallel for collapse(2) schedule(dynamic,1)

                    for (int lvl = 0; lvl < levred; lvl++) {
                        for (int dir = 1; dir < 4; dir++) {
                            int Wlvl_L = Ldecomp.level_W (lvl);
                            int Hlvl_L = Ldecomp.level_H (lvl);

                            float ** WavCoeffs_L = Ldecomp.level_coeffs (lvl);

                            madL[lvl][dir - 1] = SQR (Mad (WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                        }
                    }

                    float vari[levred];

                    if (levred == 7) {
                        edge = 2;
                        vari[0] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[1] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));

                        vari[3] = 8.f * SQR ((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[4] = 8.f * SQR ((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[5] = 8.f * SQR ((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[6] = 8.f * SQR ((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    } else if (levred == 4) {
                        edge = 3;
                        vari[0] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[1] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[3] = 8.f * SQR ((lp.noiself / 125.0) * (1.0 + lp.noiselc / 25.0));

                    }


                    if (( lp.noiself > 0.1f ||  lp.noiselc > 0.1f)) {
                        vari[0] = max (0.0001f, vari[0]);
                        vari[1] = max (0.0001f, vari[1]);
                        vari[2] = max (0.0001f, vari[2]);
                        vari[3] = max (0.0001f, vari[3]);

                        if (levred == 7) {

                            vari[4] = max (0.0001f, vari[4]);
                            vari[5] = max (0.0001f, vari[5]);
                            vari[6] = max (0.0001f, vari[6]);
                        }

                        float* noisevarlum = nullptr;  // we need a dummy to pass it to WaveletDenoiseAllL

                        WaveletDenoiseAllL (Ldecomp, noisevarlum, madL, vari, edge, numThreads);
                    }
                }


                float variC[levred];

                if (!adecomp.memoryAllocationFailed && !bdecomp.memoryAllocationFailed) {

                    if (levred == 7) {
                        edge = 2;
                        variC[0] = SQR (lp.noisecf / 10.0);
                        variC[1] = SQR (lp.noisecf / 10.0);
                        variC[2] = SQR (lp.noisecf / 10.0);

                        variC[3] = SQR (lp.noisecf / 10.0);
                        variC[4] = SQR (lp.noisecf / 10.0);
                        variC[5] = SQR (lp.noisecc / 10.0);
                        variC[6] = SQR (lp.noisecc / 10.0);
                    } else if (levred == 4) {
                        edge = 3;
                        variC[0] = SQR (lp.noisecf / 10.0);
                        variC[1] = SQR (lp.noisecf / 10.0);
                        variC[2] = SQR (lp.noisecf / 10.0);
                        variC[3] = SQR (lp.noisecf / 10.0);
                    }

                    if (( lp.noisecf > 0.1f ||  lp.noisecc > 0.1f  || noiscfactiv)) {
                        float minic = 0.0001f;

                        if (noiscfactiv) {
                            minic = 0.01f;//only for artifact shape detection
                        }

                        variC[0] = max (minic, variC[0]);
                        variC[1] = max (minic, variC[1]);
                        variC[2] = max (minic, variC[2]);
                        variC[3] = max (minic, variC[3]);

                        if (levred == 7) {

                            variC[4] = max (0.0001f, variC[4]);
                            variC[5] = max (0.0001f, variC[5]);
                            variC[6] = max (0.0001f, variC[6]);
                        }

                        float* noisevarchrom = new float[bfh * bfw];

                        for (int q = 0; q < bfh * bfw; q++) {
                            noisevarchrom[q] = 1.f;
                        }


                        float noisevarab_r = 100.f; //SQR(lp.noisecc / 10.0);
                        WaveletDenoiseAllAB (Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, false, false, false, numThreads);
                        WaveletDenoiseAllAB (Ldecomp, bdecomp, noisevarchrom, madL, variC, edge, noisevarab_r, false, false, false, numThreads);
                        delete[] noisevarchrom;

                    }
                }

                if (!Ldecomp.memoryAllocationFailed) {

                    Ldecomp.reconstruct (bufwv.L[0]);
                }

                if (!adecomp.memoryAllocationFailed) {

                    adecomp.reconstruct (bufwv.a[0]);
                }

                if (!bdecomp.memoryAllocationFailed) {

                    bdecomp.reconstruct (bufwv.b[0]);
                }

                DeNoise_Local (call, lp, original, transformed, bufwv, cx, cy);
            }

        }

//local color and light
        if (!lp.inv  && (lp.chro != 0 || lp.ligh != 0.f || lp.qualcurvemet != 0) && lp.colorena) { // || lllocalcurve)) { //interior ellipse renforced lightness and chroma  //locallutili
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            //printf("hueplus=%f huemoins=%f dhu=%f\n", hueplus, huemoins, dhue);
            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            LabImage *bufcolorig = nullptr;
            float chprov = 1.f;
            float chpro = 1.f;
            float cligh = 1.f;
            float clighL = 1.f;
            float clighmax ;
            float **buflight = nullptr;
            float **bufchro = nullptr;
            float **buflightslid = nullptr;

            int bfh = 0.f, bfw = 0.f;


            float adjustr = 1.0f;

//adapt chroma to working profile
            if      (params->icm.working == "ProPhoto")   {
                adjustr = 1.2f;   // 1.2 instead 1.0 because it's very rare to have C>170..
            } else if (params->icm.working == "Adobe RGB")  {
                adjustr = 1.8f;
            } else if (params->icm.working == "sRGB")       {
                adjustr = 2.0f;
            } else if (params->icm.working == "WideGamut")  {
                adjustr = 1.2f;
            } else if (params->icm.working == "Beta RGB")   {
                adjustr = 1.4f;
            } else if (params->icm.working == "BestRGB")    {
                adjustr = 1.4f;
            } else if (params->icm.working == "BruceRGB")   {
                adjustr = 1.8f;
            }




            if (call <= 3) { //simpleprocess, dcrop, improccoordinator
                bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                bfw = int (lp.lx + lp.lxL) + del;
                bufcolorig = new LabImage (bfw, bfh);//buffer for data in zone limit

                buflight   = new float*[bfh];//for lightness curve

                for (int i = 0; i < bfh; i++) {
                    buflight[i] = new float[bfw];
                }

                bufchro   = new float*[bfh];//for chroma curve

                for (int i = 0; i < bfh; i++) {
                    bufchro[i] = new float[bfw];
                }

                buflightslid   = new float*[bfh];//for chroma curve

                for (int i = 0; i < bfh; i++) {
                    buflightslid[i] = new float[bfw];
                }


#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufcolorig->L[ir][jr] = 0.f;
                        bufcolorig->a[ir][jr] = 0.f;
                        bufcolorig->b[ir][jr] = 0.f;
                        bufchro[ir][jr] = 0.f;
                        buflightslid[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                    }

                clighmax = 0.f;
                /*
                                int yStart = lp.yc - lp.lyT - cy;
                                int yEnd = lp.yc + lp.ly - cy;
                                int xStart = lp.xc - lp.lxL - cx;
                                int xEnd = lp.xc + lp.lx - cx;

                     //   int begx = int (lp.xc - lp.lxL);
                     //   int begy = int (lp.yc - lp.lyT);

                */

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;
                /*
                                            int yStart = lp.yc - lp.lyT - cy;
                                            int yEnd = lp.yc + lp.ly - cy;
                                            int xStart = lp.xc - lp.lxL - cx;
                                            int xEnd = lp.xc + lp.lx - cx;
                                            //there is a bug in calculation==> outof limits ==> crash
                                            printf("cy=%i cx=%i begy=%i begx=%i yS=%i yE=%i xS=%i xE=%i tH=%i tW=%i\n", cy, cx, begy, begx, yStart, yEnd, xStart, xEnd, transformed->H, transformed->W );
                                            int ymax = min(transformed->H, yEnd);
                                            int xmax = min(transformed->W, xEnd);

                            #ifdef _OPENMP
                                            #pragma omp parallel for schedule(dynamic,16)
                            #endif

                                            for (int y = yStart; y < ymax ; y++) {
                                                int loy = cy + y;

                                                for (int x = xStart, lox = cx + x; x < xmax; x++, lox++) {
                                                    bufcolorig->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                                                    bufcolorig->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                                                    bufcolorig->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas
                                                }
                                            }
                */


#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufcolorig->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufcolorig->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufcolorig->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas

                            chprov = 0.f;
                            chpro = 0.f;

                            //Chroma curve
                            if (cclocalcurve  && lp.qualcurvemet != 0) { // C=f(C) curve
                                float chromat = sqrt (SQR (bufcolorig->a[loy - begy][lox - begx]) +  SQR (bufcolorig->b[loy - begy][lox - begx]));
                                float ch;
                                float ampli = 25.f;
                                ch = (cclocalcurve[chromat * adjustr ])  / ((chromat + 0.00001f) * adjustr); //ch between 0 and 0 50 or more

                                if (ch <= 1.f) {//convert data curve near values of slider -100 + 100, to be used after to detection shape
                                    chpro = 99.f * ch - 99.f;
                                } else {
                                    chpro = CLIPCHRO (ampli * ch - ampli); //ampli = 25.f arbitrary empirical coefficient between 5 and 50
                                }

                                bufchro[loy - begy][lox - begx] = chpro;

                            }


                            //slider lightness
                            clighL = 0.f;

                            if (lp.ligh != 0.f  && lp.curvact) {
                                float lL;
                                float lighLnew;
                                float amplil = 140.f;
                                float lighL = bufcolorig->L[loy - begy][lox - begx];
                                calclight (lighL, lp.ligh , lighLnew, true);//replace L-curve
                                lL = lighLnew / lighL;

                                if (lL <= 1.f) {//convert data curve near values of slider -100 + 100, to be used after to detection shape
                                    clighL = 99.f * lL - 99.f;
                                } else {
                                    clighL = CLIPLIG (amplil * lL - amplil); //ampli = 25.f arbitrary empirical coefficient between 5 and 150
                                }

                                buflightslid[loy - begy][lox - begx] = clighL;

                            }

                            cligh = 0.f;

                            //luma curve
                            if (lllocalcurve  && lp.qualcurvemet == 2) {// L=f(L) curve enhanced
                                float lh;
                                float amplil = 25.f;
                                float lighn = bufcolorig->L[loy - begy][lox - begx];
                                lh = (lllocalcurve[lighn * 1.9f]) / ((lighn + 0.00001f) * 1.9f) ; // / ((lighn) / 1.9f) / 3.61f; //lh between 0 and 0 50 or more

                                if (lh <= 1.f) {//convert data curve near values of slider -100 + 100, to be used after to detection shape
                                    cligh = 0.3f * (100.f * lh - 100.f);//0.3 reduce sensibility
                                } else {
                                    cligh = CLIPLIG (amplil * lh - amplil);
                                }

                                buflight[loy - begy][lox - begx] = cligh;

                            }

                        }
                    }

            }


            ColorLight_Local (call, bufcolorig, buflight, bufchro, buflightslid, sp, moy, hueplus, huemoins, hueref, dhue, chromaref, lumaref, locallutili, lllocalcurve, loclhCurve, cclocalcurve, chprov, clighmax, lp, original, transformed, cx, cy);

            if (call <= 3) {

                delete bufcolorig;

                // delete bufcoltra;

                for (int i = 0; i < bfh; i++) {
                    delete [] buflight[i];
                }

                delete [] buflight;

                for (int i = 0; i < bfh; i++) {
                    delete [] bufchro[i];
                }

                delete [] bufchro;

                for (int i = 0; i < bfh; i++) {
                    delete [] buflightslid[i];
                }


                delete [] buflightslid;
            }
        }
//inverse
        else if (lp.inv  && (lp.chro != 0 || lp.ligh != 0.f) && lp.colorena) {

            InverseColorLight_Local (lp, original, transformed, cx, cy);
        }


        if (!lp.inv  && lp.cont != 0 && lp.colorena) {  //contrast interior ellipse
            const float pm = lp.cont < 0.f ? -1.f : 1.f;
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            LabImage *bufcontorig = nullptr;
            float **buflightc = nullptr;
            int bfh = 0, bfw = 0;
            float clighc = 0.f;
            const float localtype = lumaref;
            //         const float localtype = ave;
            float reducac;
            float corered;

            if (lp.sens < 30.f) {
                reducac = 0.2f * (lp.sens / 100.f);
            } else {
                float areduc = 0.6285714f; //0.44f/0.7f;
                float breduc = 0.5f - areduc;
                reducac = areduc * (lp.sens / 100.f) + breduc;
            }

            const float realcox = lco.dx, realcoy = lco.dy;

            lco.alsup = (-realcox) / (localtype / 2.f);
            lco.blsup = -lco.alsup * localtype;
            lco.alsup2 = (realcoy) / (50.f - localtype / 2.f);
            lco.blsup2 = -lco.alsup2 * localtype;
            lco.alsup3 = (realcoy) / (localtype / 2.f - 50.f);
            lco.blsup3 = -lco.alsup3 * 100.f;
            lco.aDY = realcoy;

            lco.alinf = realcox / (localtype / 2.f);
            const float vi = (localtype / 2.f) / 100.f;
            const float vinf = (50.f + localtype / 2.f) / 100.f;
            ImProcFunctions::secondeg_begin (reducac, vi, lco.aa, lco.bb);//parabolic
            ImProcFunctions::secondeg_end (reducac, vinf, lco.aaa, lco.bbb, lco.ccc);//parabolic

            if (call <= 3) { //simpleprocess, dcrop, improccoordinator
                bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                bfw = int (lp.lx + lp.lxL) + del;
                bufcontorig = new LabImage (bfw, bfh);//buffer for data in zone limit

                buflightc   = new float*[bfh];//for lightness curve

                for (int i = 0; i < bfh; i++) {
                    buflightc[i] = new float[bfw];
                }



#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufcontorig->L[ir][jr] = 0.f;
                        //   bufcontorig->a[ir][jr] = 0.f;
                        //   bufcontorig->b[ir][jr] = 0.f;
                        buflightc[ir][jr] = 0.f;


                    }

                float localty;
                localty = localtype;

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;
                //  float maxc = -10000.f;
                //  float minc = 100000.f;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufcontorig->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas

                            //slider contrast
                            clighc = 1.f;
                            corered = 0.f;


                            if (lp.cont != 0.f  && lp.curvact) {

                                float cL;
                                float amplil = 150.f;
                                float prov100 = bufcontorig->L[loy - begy][lox - begx] / 32768.f;
                                float prov = prov100 * 100.f;
                                cL = 1.f;

                                if (prov > localty) {
                                    if (prov >= localty && prov < 50.f + localty / 2.f) {
                                        float core = (lco.alsup2 * prov + lco.blsup2) ;
                                        corered = prov + pm * (prov - localty) * (core);
                                    } else {
                                        float core = lco.aDY * (lco.aaa * prov100 * prov100 + lco.bbb * prov100 + lco.ccc);
                                        corered = prov + pm * (prov - localty) * (core);

                                    }
                                }  else {
                                    if (2.f * prov > localty && prov < localty)  {
                                        float core = (lco.alsup * prov + lco.blsup) ;
                                        corered = prov - pm * (localty - prov) * core;
                                    } else if (2.f * prov <= localty) {
                                        float core = prov * lco.alinf * (lco.aa * prov100 * prov100 + lco.bb * prov100);
                                        corered = prov - pm * (localty - prov) * core;

                                    }
                                }

                                cL = corered / prov;

                                if (cL <= 1.f) {//convert data curve near values of slider -100 + 100, to be used after to detection shape
                                    clighc = 99.f * cL - 99.f;
                                } else {
                                    clighc = CLIPLIG (amplil * cL - amplil); //arbitrary empirical coefficient between 5 and 150
                                }

                                /*

                                                                if (clighc > maxc) {
                                                                    maxc = clighc;
                                                                }

                                                                if (clighc < minc) {
                                                                    minc = clighc;
                                                                }
                                */
                                buflightc[loy - begy][lox - begx] = clighc;

                            }
                        }
                    }

                //       printf ("min=%2.2f max=%2.2f", minc, maxc);


            }

            Contrast_Local (call, ave, bufcontorig, buflightc, moy, hueplus, huemoins, hueref, dhue, chromaref, pm, lco, lumaref, av, lp, original, transformed, cx, cy);

            if (call <= 3) {

                delete bufcontorig;


                for (int i = 0; i < bfh; i++) {
                    delete [] buflightc[i];
                }

                delete [] buflightc;
            }


        } else if (lp.inv && lp.cont != 0 && lp.colorena) {

            float multL = (float)lp.cont * (maxl - 1.f) / 100.f + 1.f;
            float multH = (float) lp.cont * (maxh - 1.f) / 100.f + 1.f;

            lco.ah = (multH - 1.f) / (av - 100.f); //av ==> lumaref
            lco.bh = 1.f - 100.f * lco.ah;
            lco.al = (multL - 1.f) / av;
            lco.bl = 1.f;

            InverseContrast_Local (ave, lco, lp, original, transformed, cx, cy);
        }


// end contrast interior and exterior

//Tone mapping

//&& lp.tonemapena
        if (lp.strengt != 0.f  && lp.tonemapena) {
            LabImage *tmp1 = nullptr;
            float **buflight = nullptr;

            LabImage *bufgb = nullptr;
            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;

            if (call <= 3) { //simpleprocess dcrop improcc

                bufgb = new LabImage (bfw, bfh);
                buflight   = new float*[bfh];//for lightness reti

                for (int i = 0; i < bfh; i++) {
                    buflight[i] = new float[bfw];
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufgb->L[ir][jr] = 0.f;
                        bufgb->a[ir][jr] = 0.f;
                        bufgb->b[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                    }

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufgb->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufgb->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufgb->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas
                        }
                    }

                tmp1 = new LabImage (bfw, bfh);
                ImProcFunctions::EPDToneMaplocal (bufgb, tmp1, 5 , 1);
            } /*else { //stay here in case of

                tmp = new LabImage (transformed->W, transformed->H);
                tmp->CopyFrom (original);
                tmp1 = new LabImage (transformed->W, transformed->H);
                ImProcFunctions::EPDToneMaplocal (tmp, tmp1, 5 , sk);
                delete tmp;
            }
*/
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            int begy = lp.yc - lp.lyT;
            int begx = lp.xc - lp.lxL;
            int yEn = lp.yc + lp.ly;
            int xEn = lp.xc + lp.lx;

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int y = 0; y < transformed->H ; y++) //{
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {

                        float rL = CLIPRET ((tmp1->L[loy - begy][lox - begx] - original->L[y][x]) / 400.f);
                        /*
                                                if (rL > maxc) {
                                                    maxc = rL;
                                                }

                                                if (rL < minc) {
                                                    minc = rL;
                                                }
                        */

                        buflight[loy - begy][lox - begx]  = rL;

                    }
                }

//            printf ("min=%2.2f max=%2.2f", minc, maxc);

            TM_Local (call, sp, tmp1, buflight, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, original, transformed, cx, cy);

            if (call <= 3) {
                delete bufgb;

                for (int i = 0; i < bfh; i++) {
                    delete [] buflight[i];
                }

                delete [] buflight;

            }

            delete tmp1;


        }

//begin cbdl
        if ((lp.mulloc[0] != 1.f || lp.mulloc[1] != 1.f || lp.mulloc[2] != 1.f || lp.mulloc[3] != 1.f || lp.mulloc[4] != 1.f) && lp.cbdlena) {
            float **bufsh = nullptr;//buffer por square zone
            float **loctemp = nullptr;
            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;
            float b_l = -5.f;
            float t_l = 25.f;
            float t_r = 120.f;
            float b_r = 170.f;
            double skinprot = 0.;
            int choice = 0;
            float **buflight = nullptr;


            if (call <= 3) { //call from simpleprocess dcrop improcc
                bufsh   = new float*[bfh];

                for (int i = 0; i < bfh; i++) {
                    bufsh[i] = new float[bfw];
                }

                buflight   = new float*[bfh];//for lightness reti

                for (int i = 0; i < bfh; i++) {
                    buflight[i] = new float[bfw];
                }


#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufsh[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;

                    }


                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufsh[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                        }
                    }

                loctemp = new float*[bfh];//allocate temp

                for (int i = 0; i < bfh; i++) {
                    loctemp[i] = new float[bfw];
                }

                ImProcFunctions::cbdl_local_temp (bufsh, bufsh, loctemp, bfw, bfh, lp.mulloc, lp.threshol, skinprot, false,  b_l, t_l, t_r, b_r, choice, sk);


#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            float rL;
                            rL = CLIPRET ((loctemp[loy - begy][lox - begx] - original->L[y][x]) / 330.f);
                            /*
                                                        if (rL > maxc) {
                                                            maxc = rL;
                                                        }

                                                        if (rL < minc) {
                                                            minc = rL;
                                                        }
                            */

                            buflight[loy - begy][lox - begx]  = rL;

                        }
                    }

//                printf ("min=%2.2f max=%2.2f", minc, maxc);


            }


            // I initialize these variable in case of !
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            cbdl_Local (call, sp, buflight, loctemp, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, original, transformed, cx, cy);

            if (call <=  3) {
                for (int i = 0; i < bfh; i++) {
                    delete [] loctemp[i];
                }

                delete [] loctemp;

                for (int i = 0; i < bfh; i++) {
                    delete [] bufsh[i];
                }

                delete [] bufsh;

                for (int i = 0; i < bfh; i++) {
                    delete [] buflight[i];
                }

                delete [] buflight;

            } /* else {

                for (int i = 0; i < GH; i++) {
                    delete [] loctemp[i];
                }

                delete [] loctemp;

            }
*/


        }

//       }

//end cbdl
        if (!lp.invshar && lp.shrad > 0.42 && call < 3  && lp.sharpena) { //interior ellipse for sharpening, call = 1 and 2 only with Dcrop and simpleprocess
            int bfh = call == 2 ? int (lp.ly + lp.lyT) + del : original->H; //bfw bfh real size of square zone
            int bfw = call == 2 ? int (lp.lx + lp.lxL) + del : original->W;
            const JaggedArray<float> loctemp (bfw, bfh);

            if (call == 2) { //call from simpleprocess
                const JaggedArray<float> bufsh (bfw, bfh, true);
                const JaggedArray<float> hbuffer (bfw, bfh);

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufsh[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                        }
                    }

                //   }

                //sharpen only square area instaed of all image
                ImProcFunctions::deconvsharpeningloc (bufsh, hbuffer, bfw, bfh, loctemp, params->locallab.shardamping, (double)params->locallab.sharradius / 100., params->locallab.shariter, params->locallab.sharamount);
            } else { //call from dcrop.cc

                ImProcFunctions::deconvsharpeningloc (original->L, shbuffer, bfw, bfh, loctemp, params->locallab.shardamping, (double)params->locallab.sharradius / 100., params->locallab.shariter, params->locallab.sharamount);

            }

            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            //sharpen ellipse and transition
            Sharp_Local (call, sp, loctemp, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, original, transformed, cx, cy);

        } else if (lp.invshar && lp.shrad > 0.42 && call < 3 && lp.sharpena) {
            int GW = original->W;
            int GH = original->H;
            const JaggedArray<float> loctemp (GW, GH);

            ImProcFunctions::deconvsharpeningloc (original->L, shbuffer, GW, GH, loctemp, params->locallab.shardamping, (double)params->locallab.sharradius / 100., params->locallab.shariter, params->locallab.sharamount);

            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            InverseSharp_Local (sp, loctemp, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, original, transformed, cx, cy);
        }

        //      }
//&& lp.retiena
        if (lp.str > 0.f  && lp.retiena) {
            int GW = transformed->W;
            int GH = transformed->H;

            LabImage *bufreti = nullptr;
            float **buflight = nullptr;
            float **bufchro = nullptr;
            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;

            float hueplus = hueref + dhueret;
            float huemoins = hueref - dhueret;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhueret - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhueret + 2.f * rtengine::RT_PI;
            }

            int Hd, Wd;
            Hd = GH;
            Wd = GW;

            if (!lp.invret && call <= 3) {

                Hd = bfh;
                Wd = bfw;
                bufreti = new LabImage (bfw, bfh);
                buflight   = new float*[bfh];//for lightness reti

                for (int i = 0; i < bfh; i++) {
                    buflight[i] = new float[bfw];
                }

                bufchro   = new float*[bfh];//for chroma reti

                for (int i = 0; i < bfh; i++) {
                    bufchro[i] = new float[bfw];
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufreti->L[ir][jr] = 0.f;
                        bufreti->a[ir][jr] = 0.f;
                        bufreti->b[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufchro[ir][jr] = 0.f;
                    }

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufreti->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufreti->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufreti->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas
                        }
                    }



            }

            float *orig[Hd] ALIGNED16;
            float *origBuffer = new float[Hd * Wd];

            for (int i = 0; i < Hd; i++) {
                orig[i] = &origBuffer[i * Wd];
            }

            float *orig1[Hd] ALIGNED16;
            float *origBuffer1 = new float[Hd * Wd];

            for (int i = 0; i < Hd; i++) {
                orig1[i] = &origBuffer1[i * Wd];
            }


            LabImage *tmpl = nullptr;

            if (!lp.invret && call <= 3) {


#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int ir = 0; ir < Hd; ir += 1)
                    for (int jr = 0; jr < Wd; jr += 1) {
                        orig[ir][jr] = bufreti->L[ir][jr];
                        orig1[ir][jr] = bufreti->L[ir][jr];
                    }

                tmpl = new LabImage (Wd, Hd);

            } /* else {

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int ir = 0; ir < Hd; ir += 1)
                    for (int jr = 0; jr < Wd; jr += 1) {
                        orig[ir][jr] = original->L[ir][jr];
                        orig1[ir][jr] = transformed->L[ir][jr];
                    }

                tmpl = new LabImage (transformed->W, transformed->H);


            }
*/
            float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            ImProcFunctions::MSRLocal (orig, tmpl->L, orig1, Wd, Hd, params->locallab, sk, locRETgainCcurve, 0, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int ir = 0; ir < Hd; ir += 1)
                for (int jr = 0; jr < Wd; jr += 1) {
                    tmpl->L[ir][jr] = orig[ir][jr];

                    if (!lp.invret) {
                        float rL;
                        rL = CLIPRET ((tmpl->L[ir][jr] - bufreti->L[ir][jr]) / 328.f);
                        /*
                                                if (rL > maxc) {
                                                    maxc = rL;
                                                }

                                                if (rL < minc) {
                                                    minc = rL;
                                                }
                        */
                        buflight[ir][jr] = rL;
                    }
                }

//            printf ("min=%2.2f max=%2.2f", minc, maxc);

//new shape detection


            if (!lp.invret) {

                Reti_Local (call, buflight, bufchro, hueplus, huemoins, hueref, dhueret, chromaref, lumaref, lp, original, transformed, tmpl, cx, cy, 0);
            } else {
                InverseReti_Local (lp, original, transformed, tmpl, cx, cy, 0);
            }

            if (params->locallab.chrrt > 0) {

                if (!lp.invret && call <= 3) {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int ir = 0; ir < Hd; ir += 1)
                        for (int jr = 0; jr < Wd; jr += 1) {

                            orig[ir][jr] = sqrt (SQR (bufreti->a[ir][jr]) + SQR (bufreti->b[ir][jr]));
                            orig1[ir][jr] = sqrt (SQR (bufreti->a[ir][jr]) + SQR (bufreti->b[ir][jr]));
                        }

                } /* else {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int ir = 0; ir < GH; ir += 1)
                        for (int jr = 0; jr < GW; jr += 1) {
                            orig[ir][jr] = sqrt (SQR (original->a[ir][jr]) + SQR (original->b[ir][jr]));
                            orig1[ir][jr] = sqrt (SQR (transformed->a[ir][jr]) + SQR (transformed->b[ir][jr]));
                        }
                }
*/
                ImProcFunctions::MSRLocal (orig, tmpl->L, orig1, Wd, Hd, params->locallab, sk, locRETgainCcurve, 1, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);

                if (!lp.invret && call <= 3) {


#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int ir = 0; ir < Hd; ir += 1)
                        for (int jr = 0; jr < Wd; jr += 1) {
                            float Chprov = orig1[ir][jr];
                            float2 sincosval;
                            sincosval.y = Chprov == 0.0f ? 1.f : bufreti->a[ir][jr] / Chprov;
                            sincosval.x = Chprov == 0.0f ? 0.f : bufreti->b[ir][jr] / Chprov;
                            tmpl->a[ir][jr] = orig[ir][jr] * sincosval.y;
                            tmpl->b[ir][jr] = orig[ir][jr] * sincosval.x;


                            if (!lp.invret) {

                                float ra;
                                ra = CLIPRET ((sqrt (SQR (tmpl->a[ir][jr]) + SQR (tmpl->b[ir][jr])) - Chprov) / 300.f);
                                /*
                                                                if (ra > maxch) {
                                                                    maxch = ra;
                                                                }

                                                                if (ra < minch) {
                                                                    minch = ra;
                                                                }
                                */
                                bufchro[ir][jr] = ra;
                            }

                        }

//                    printf ("minch=%2.2f maxch=%2.2f", minch, maxch);


                } /* else {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int ir = 0; ir < Hd; ir += 1)
                        for (int jr = 0; jr < Wd; jr += 1) {
                            float Chprov = orig1[ir][jr];
                            float2 sincosval;
                            sincosval.y = Chprov == 0.0f ? 1.f : transformed->a[ir][jr] / Chprov;
                            sincosval.x = Chprov == 0.0f ? 0.f : transformed->b[ir][jr] / Chprov;
                            tmpl->a[ir][jr] = orig[ir][jr] * sincosval.y;
                            tmpl->b[ir][jr] = orig[ir][jr] * sincosval.x;

                        }
                }
*/

                if (!lp.invret) {

                    Reti_Local (call, buflight, bufchro, hueplus, huemoins, hueref, dhueret, chromaref, lumaref, lp, original, transformed, tmpl, cx, cy, 1);
                } else {
                    InverseReti_Local (lp, original, transformed, tmpl, cx, cy, 1);
                }

            }

            delete tmpl;
            delete [] origBuffer;
            delete [] origBuffer1;

            if (!lp.invret && call <= 3) {

                delete  bufreti;

                for (int i = 0; i < bfh; i++) {
                    delete [] buflight[i];
                }

                delete [] buflight;

                for (int i = 0; i < bfh; i++) {
                    delete [] bufchro[i];
                }

                delete [] bufchro;

            }
        }


// Gamut and Munsell control - very important do not desactivated to avoid crash
        if (params->locallab.avoid) {
            TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.working);
            float wip[3][3] = {
                {static_cast<float> (wiprof[0][0]), static_cast<float> (wiprof[0][1]), static_cast<float> (wiprof[0][2])},
                {static_cast<float> (wiprof[1][0]), static_cast<float> (wiprof[1][1]), static_cast<float> (wiprof[1][2])},
                {static_cast<float> (wiprof[2][0]), static_cast<float> (wiprof[2][1]), static_cast<float> (wiprof[2][2])}
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
                vfloat c327d68v = F2V (327.68f);
                vfloat onev = F2V (1.f);
#endif

#ifdef _OPENMP
#ifdef _DEBUG
                #pragma omp for schedule(dynamic,16) firstprivate(MunsDebugInfo)
#else
                #pragma omp for schedule(dynamic,16)
#endif
#endif

                for (int y = 0; y < transformed->H; y++) {
#ifdef __SSE2__
                    int i = 0;

                    for (; i < transformed->W - 3; i += 4) {
                        vfloat av = LVFU (transformed->a[y][i]);
                        vfloat bv = LVFU (transformed->b[y][i]);

                        if (needHH) { // only do expensive atan2 calculation if needed
                            STVF (atan2Buffer[i], xatan2f (bv, av));
                        }

                        vfloat Chprov1v = vsqrtf (SQRV (bv) + SQRV (av));
                        STVF (sqrtBuffer[i], Chprov1v / c327d68v);
                        vfloat sincosyv = av / Chprov1v;
                        vfloat sincosxv = bv / Chprov1v;
                        vmask selmask = vmaskf_eq (Chprov1v, ZEROV);
                        sincosyv = vself (selmask, onev, sincosyv);
                        sincosxv = vselfnotzero (selmask, sincosxv);
                        STVF (sincosyBuffer[i], sincosyv);
                        STVF (sincosxBuffer[i], sincosxv);
                    }

                    for (; i < transformed->W; i++) {
                        float aa = transformed->a[y][i];
                        float bb = transformed->b[y][i];

                        if (needHH) { // only do expensive atan2 calculation if needed
                            atan2Buffer[i] = xatan2f (bb, aa);
                        }

                        float Chprov1 = sqrtf (SQR (bb) + SQR (aa));
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
                        float Lprov1 = transformed->L[y][x] / 327.68f;
                        float2 sincosval;
#ifdef __SSE2__
                        float HH = atan2Buffer[x]; // reading HH from line buffer even if line buffer is not filled is faster than branching
                        float Chprov1 = sqrtBuffer[x];
                        sincosval.y = sincosyBuffer[x];
                        sincosval.x = sincosxBuffer[x];
                        float chr = 0.f;

#else
                        float aa = transformed->a[y][x];
                        float bb = transformed->b[y][x];
                        float HH = 0.f, chr = 0.f;

                        if (needHH) { // only do expensive atan2 calculation if needed
                            HH = xatan2f (bb, aa);
                        }

                        float Chprov1 = sqrtf (SQR (aa) + SQR (bb)) / 327.68f;

                        if (Chprov1 == 0.0f) {
                            sincosval.y = 1.f;
                            sincosval.x = 0.0f;
                        } else {
                            sincosval.y = aa / (Chprov1 * 327.68f);
                            sincosval.x = bb / (Chprov1 * 327.68f);
                        }

#endif

#ifdef _DEBUG
                        bool neg = false;
                        bool more_rgb = false;
// Color::pregamutlab (Lprov1, HH, chr);
                        Chprov1 = min (Chprov1, chr);

                        Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.92f, neg, more_rgb);
#else
                        Color::pregamutlab (Lprov1, HH, chr);
                        Chprov1 = min (Chprov1, chr);
                        Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.92f);
#endif

                        transformed->L[y][x] = Lprov1 * 327.68f;
                        transformed->a[y][x] = 327.68f * Chprov1 * sincosval.y;
                        transformed->b[y][x] = 327.68f * Chprov1 * sincosval.x;

                        if (needHH) {
                            float Lprov2 = original->L[y][x] / 327.68f;
                            float correctionHue = 0.f; // Munsell's correction
                            float correctlum = 0.f;
                            float memChprov = sqrtf (SQR (original->a[y][x]) + SQR (original->b[y][x])) / 327.68f;
                            float Chprov = sqrtf (SQR (transformed->a[y][x]) + SQR (transformed->b[y][x])) / 327.68f;
#ifdef _DEBUG
                            Color::AllMunsellLch (true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum, MunsDebugInfo);
#else
                            Color::AllMunsellLch (true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum);
#endif

                            if (fabs (correctionHue) < 0.015f) {
                                HH += correctlum;    // correct only if correct Munsell chroma very little.
                            }

                            float2 sincosval = xsincosf (HH + correctionHue);

                            transformed->a[y][x] = 327.68f * Chprov * sincosval.y; // apply Munsell
                            transformed->b[y][x] = 327.68f * Chprov * sincosval.x;
                        }
                    }
                }
            }
        }


#ifdef _DEBUG

        if (settings->verbose) {
            t2e.set();
            printf ("Color::AllMunsellLch (correction performed in %d usec):\n", t2e.etime (t1e));
            //  printf("   Munsell chrominance: MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhue[0],    MunsDebugInfo->maxdhue[1],    MunsDebugInfo->maxdhue[2],    MunsDebugInfo->maxdhue[3],    MunsDebugInfo->depass);
            //  printf("   Munsell luminance  : MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhuelum[0], MunsDebugInfo->maxdhuelum[1], MunsDebugInfo->maxdhuelum[2], MunsDebugInfo->maxdhuelum[3], MunsDebugInfo->depassLum);
        }

        delete MunsDebugInfo;
#endif

    }

}

}
