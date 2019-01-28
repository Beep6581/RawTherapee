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
#include <fftw3.h>

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
#include "../rtgui/thresholdselector.h"
#include "median.h"

#include "cplx_wavelet_dec.h"
#include "ciecam02.h"

#define BENCHMARK
#include "StopWatch.h"
#include "rt_algo.h"
#include "guidedfilter.h"

#define cliploc( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )

#define TS 64       // Tile size
#define offset 25   // shift between tiles
#define fTS ((TS/2+1))  // second dimension of Fourier tiles
#define blkrad 1    // radius of block averaging

#define epsilon 0.001f/(TS*TS) //tolerance


#define CLIPC(a) ((a)>-42000?((a)<42000?(a):42000):-42000)  // limit a and b  to 130 probably enough ?
#define CLIPL(x) LIM(x,0.f,40000.f) // limit L to about L=120 probably enough ?
#define CLIPLOC(x) LIM(x,0.f,32767.f)
#define CLIPLIG(x) LIM(x,-99.5f, 99.5f)
#define CLIPCHRO(x) LIM(x,0.f, 140.f)
#define CLIPRET(x) LIM(x,-99.5f, 99.5f)
#define CLIP1(x) LIM(x, 0.f, 1.f)
#pragma GCC diagnostic warning "-Wextra"


namespace
{

float calcLocalFactor(const float lox, const float loy, const float lcx, const float dx, const float lcy, const float dy, const float ach)
{
//elipse x2/a2 + y2/b2=1
//transition elipsoidal
//x==>lox y==>loy
// a==> dx  b==>dy

    float kelip = dx / dy;
    float belip = sqrt((rtengine::SQR((lox - lcx) / kelip) + rtengine::SQR(loy - lcy)));    //determine position ellipse ==> a and b
    float aelip = belip * kelip;
    float degrad = aelip / dx;
    float ap = rtengine::RT_PI / (1.f - ach);
    float bp = rtengine::RT_PI - ap;
    return 0.5f * (1.f + xcosf(degrad * ap + bp));  //trigo cos transition

}
float calcLocalFactorrect(const float lox, const float loy, const float lcx, const float dx, const float lcy, const float dy, const float ach)
{
    float eps = 0.0001f;
    float krap = fabs(dx / dy);
    float kx = (lox - lcx);
    float ky = (loy - lcy);
    float ref = 0.f;

    if (fabs(kx / (ky + eps)) < krap) {
        ref = sqrt(rtengine::SQR(dy) * (1.f + rtengine::SQR(kx / (ky + eps))));
    } else {
        ref = sqrt(rtengine::SQR(dx) * (1.f + rtengine::SQR(ky / (kx + eps))));
    }

    float rad = sqrt(rtengine::SQR(kx) + rtengine::SQR(ky));
    float coef = rad / ref;
    float ac = 1.f / (ach - 1.f);
    float fact = ac * (coef - 1.f);
    return fact;

}

/*
float calcLocalFactorinv (const float lox, const float loy, const float lcx, const float dx, const float lcy, const float dy, const float ach)
{
//elipse x2/a2 + y2/b2=1
//transition elipsoidal
//x==>lox y==>loy
// a==> dx  b==>dy

    float kelip = dx / dy;
    float belip = sqrt ((rtengine::SQR ((lox - lcx) / kelip) + rtengine::SQR (loy - lcy))); //determine position ellipse ==> a and b
    float aelip = belip * kelip;
    float degrad = aelip / dx;
    float ap = rtengine::RT_PI / (ach);
//    float bp = rtengine::RT_PI - ap;
    return 0.5f * (1.f + xcosf (degrad * ap)); //trigo cos transition

}

*/
}

namespace rtengine


{
extern MyMutex *fftwMutex;


using namespace procparams;



extern const Settings* settings;


struct local_params {
    float yc, xc;
    float ycbuf, xcbuf;
    float lx, ly;
    float lxL, lyT;
    float dxx, dyy;
    float iterat;
    int cir;
    float thr;
    int prox;
    int chro, cont, sens, sensh, senscb, sensbn, senstm, sensex, sensexclu, sensden, senslc, senssf;
    float struco;
    float blendmacol;
    float blendmaexp;
    float struexp;
    float blurexp;
    float blurcol;
    float ligh;
    int shamo, shdamp, shiter, senssha, sensv;
    float neig;
    float strng;
    float lcamount;
    double shrad;
    double shblurr;
    double rad;
    double stren;
    int trans;
    int dehaze;
    bool inv;
    bool curvact;
    bool invrad;
    bool invret;
    bool invshar;
    bool showmaskexpo;
    bool actsp;
    float str;
    int qualmet;
    int qualcurvemet;
    int showmaskcolmet;
    int showmaskexpmet;
    int blurmet;
    float noiself;
    float noiseldetail;
    int noiselequal;
    float noisechrodetail;
    float bilat;
    float noiselc;
    float noisecf;
    float noisecc;
    float mulloc[5];
    float threshol;
    float chromacb;
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
    bool lcena;
    bool sfena;
    bool cbdlena;
    bool denoiena;
    bool expvib;
    bool exposena;
    bool cut_past;
    float past;
    float satur;
    int blac;
    int shcomp;
    int hlcomp;
    int hlcompthr;
    double expcomp;
    int excmet;
    int strucc;
    int war;
    float adjch;
    int shapmet;
    bool enaColorMask;
    bool enaExpMask;
    bool enaColMask;
};

static void SobelCannyLuma(float **sobelL, float **luma, int bfw, int bfh, float radius)
{
    //base of the process to detect shape in complement of deltaE
    //use for calcualte Spot reference
    // and for structure of the shape
    // actually , as thr program don't use these function, I just create a simple "Canny" near of Sobel. This can be completed after with teta, etc.
    float *tmLBuffer = new float[bfh * bfw];
    float *tmL[bfh];

    for (int i = 0; i < bfh; i++) {
        tmL[i] = &tmLBuffer[i * bfw];
    }


    int GX[3][3];
    int GY[3][3];
    float SUML;

    float sumXL, sumYL;

    //Sobel Horizontal
    GX[0][0] = 1;
    GX[0][1] = 0;
    GX[0][2] = -1;
    GX[1][0] = 2;
    GX[1][1] = 0;
    GX[1][2] = -2;
    GX[2][0] = 1;
    GX[2][1] = 0;
    GX[2][2] = -1;

    //Sobel Vertical
    GY[0][0] =  1;
    GY[0][1] =  2;
    GY[0][2] =  1;
    GY[1][0] =  0;
    GY[1][1] =  0;
    GY[1][2] =  0;
    GY[2][0] = -1;
    GY[2][1] = -2;
    GY[2][2] = -1;
    //inspired from Chen Guanghua Zhang Xiaolong
//  gaussianBlur (luma, tmL, bfw, bfh, radius);


    {
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw ; x++) {
                sobelL[y][x] = 0.f;
                tmL[y][x] = luma[y][x];
            }
        }

        if (radius > 0.f) {
            radius /= 2.f;

            if (radius < 0.5f) {
                radius = 0.5f;
            }

            gaussianBlur(luma, tmL, bfw, bfh, radius);
        }

//}
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw ; x++) {
                sumXL    = 0.f;
                sumYL    = 0.f;

                if (y == 0 || y == bfh - 1) {
                    SUML = 0.f;
                } else if (x == 0 || x == bfw - 1) {
                    SUML = 0.f;
                } else {
                    for (int i = -1; i < 2; i++) {
                        for (int j = -1; j < 2; j++) {
                            sumXL += GX[j + 1][i + 1] * tmL[y + i][x + j];
                        }
                    }

                    for (int i = -1; i < 2; i++) {
                        for (int j = -1; j < 2; j++) {
                            sumYL += GY[j + 1][i + 1] * tmL[y + i][x + j];
                        }
                    }

                    //Edge strength
                    SUML = sqrt(SQR(sumXL) + SQR(sumYL));
                    //we can add if need teta = atan2 (sumYr, sumXr)
                }

                SUML = CLIPLOC(SUML);

                sobelL[y][x] = SUML;
            }
        }

    }

    delete [] tmLBuffer;
    /*
        //mean to exclude litlle values
        for (int y = 1; y < bfh - 1 ; y++) {
            for (int x = 1; x < bfw - 1 ; x++) {
                sobelL[y][x] = (sobelL[y - 1][x - 1] + sobelL[y - 1][x] + sobelL[y - 1][x + 1] + sobelL[y][x - 1] + sobelL[y][x] + sobelL[y][x + 1] + sobelL[y + 1][x - 1] + sobelL[y + 1][x] + sobelL[y + 1][x + 1]) / 9;
            }
        }
    */

}



static void calcLocalParams(int sp, int oW, int oH, const LocallabParams& locallab, struct local_params& lp, int llColorMask, int llExpMask)
{
    int w = oW;
    int h = oH;
    int circr = locallab.spots.at(sp).circrad;
    float streng = ((float)locallab.spots.at(sp).stren) / 100.f;
    float gam = ((float)locallab.spots.at(sp).gamma) / 100.f;
    float est = ((float)locallab.spots.at(sp).estop) / 100.f;
    float scal_tm = ((float)locallab.spots.at(sp).scaltm) / 10.f;
    float rewe = ((float)locallab.spots.at(sp).rewei);
    float strlight = ((float)locallab.spots.at(sp).streng) / 100.f;
    float thre = locallab.spots.at(sp).thresh / 100.f;
    double local_x = locallab.spots.at(sp).locX / 2000.0;
    double local_y = locallab.spots.at(sp).locY / 2000.0;
    double local_xL = locallab.spots.at(sp).locXL / 2000.0;
    double local_yT = locallab.spots.at(sp).locYT / 2000.0;
    double local_center_x = locallab.spots.at(sp).centerX / 2000.0 + 0.5;
    double local_center_y = locallab.spots.at(sp).centerY / 2000.0 + 0.5;
    double local_center_xbuf = 0.0; // Provision
    double local_center_ybuf = 0.0; // Provision
    double local_dxx = locallab.spots.at(sp).iter / 8000.0; //for proxi = 2==> # 1 pixel
    double local_dyy = locallab.spots.at(sp).iter / 8000.0;
    float iterati = (float) locallab.spots.at(sp).iter;
//    double local_dyy = locallab.proxi;
    float neigh = float (locallab.spots.at(sp).neigh);
    float chromaPastel = float (locallab.spots.at(sp).pastels)   / 100.0f;
    float chromaSatur  = float (locallab.spots.at(sp).saturated) / 100.0f;
    int local_sensiv = locallab.spots.at(sp).sensiv;
    int local_sensiex = locallab.spots.at(sp).sensiex;

    if (locallab.spots.at(sp).qualityMethod == "std") {
        lp.qualmet = 0;
    } else if (locallab.spots.at(sp).qualityMethod == "enh") {
        lp.qualmet = 1;
    } else if (locallab.spots.at(sp).qualityMethod == "enhden") {
        lp.qualmet = 2;
    }

    if (locallab.spots.at(sp).qualitycurveMethod == "none") {
        lp.qualcurvemet = 0;
    } else if (locallab.spots.at(sp).qualitycurveMethod == "std") {
        lp.qualcurvemet = 1;
    } /* else if (locallab.spots.at(sp).qualitycurveMethod == "enh") {

        lp.qualcurvemet = 2;
    } else if (locallab.spots.at(sp).qualitycurveMethod == "enhsup") {
        lp.qualcurvemet = 3;
    } else if (locallab.spots.at(sp).qualitycurveMethod == "contr") {
        lp.qualcurvemet = 4;
    } else if (locallab.spots.at(sp).qualitycurveMethod == "sob2") {
        lp.qualcurvemet = 5;
    }
*/
    lp.showmaskcolmet = llColorMask;
    lp.showmaskexpmet = llExpMask;
    lp.enaColorMask = locallab.spots.at(sp).enaColorMask && llColorMask == 0 && llExpMask == 0; // Color & Light mask is deactivated if Exposure mask is visible
    lp.enaExpMask = locallab.spots.at(sp).enaExpMask && llExpMask == 0 && llColorMask == 0; // Exposure mask is deactivated if Color & Light mask is visible

    if (locallab.spots.at(sp).blurMethod == "norm") {
        lp.blurmet = 0;
    } else if (locallab.spots.at(sp).blurMethod == "inv") {
        lp.blurmet = 1;
    } else if (locallab.spots.at(sp).blurMethod == "sym") {
        lp.blurmet = 2;
    }

    if (locallab.spots.at(sp).spotMethod == "norm") {
        lp.excmet = 0;
    } else if (locallab.spots.at(sp).spotMethod == "exc") {
        lp.excmet = 1;
    }

    if (locallab.spots.at(sp).shape == "ELI") {
        lp.shapmet = 0;
    } else if (locallab.spots.at(sp).shape == "RECT") {
        lp.shapmet = 1;
    }

    float local_noiself = (float)locallab.spots.at(sp).noiselumf;
    float local_noiselc = (float)locallab.spots.at(sp).noiselumc;
    float local_noiseldetail = (float)locallab.spots.at(sp).noiselumdetail;
    int local_noiselequal = locallab.spots.at(sp).noiselequal;
    float local_noisechrodetail = (float)locallab.spots.at(sp).noisechrodetail;
    int local_sensiden = locallab.spots.at(sp).sensiden;

    float local_noisecf = ((float)locallab.spots.at(sp).noisechrof) / 10.f;
    float local_noisecc = ((float)locallab.spots.at(sp).noisechroc) / 10.f;
    float multi[5];

    for (int y = 0; y < 5; y++) {
        multi[y] = ((float) locallab.spots.at(sp).mult[y]) / 100.f;
    }

    float thresho = ((float)locallab.spots.at(sp).threshold) / 100.f;
    float chromcbdl = (float)locallab.spots.at(sp).chromacbdl;

    int local_chroma = locallab.spots.at(sp).chroma;
    int local_sensi = locallab.spots.at(sp).sensi;
    int local_sensibn = locallab.spots.at(sp).sensibn;
    int local_sensitm = locallab.spots.at(sp).sensitm;
    int local_sensiexclu = locallab.spots.at(sp).sensiexclu;
    int local_sensilc = locallab.spots.at(sp).sensilc;
    int local_struc = locallab.spots.at(sp).struc;
    int local_warm = locallab.spots.at(sp).warm;
    int local_sensih = locallab.spots.at(sp).sensih;
    int local_dehaze = locallab.spots.at(sp).dehaz;
    int local_sensicb = locallab.spots.at(sp).sensicb;
    int local_contrast = locallab.spots.at(sp).contrast;
    float local_lightness = (float) locallab.spots.at(sp).lightness;
    float structcolor = (float) locallab.spots.at(sp).structcol;
    float blendmaskcolor = ((float) locallab.spots.at(sp).blendmaskcol) / 100.f ;
    float blendmaskexpo = ((float) locallab.spots.at(sp).blendmaskexp) / 100.f ;
    float structexpo = (float) locallab.spots.at(sp).structexp;
    float blurexpo = (float) locallab.spots.at(sp).blurexpde;
    float blurcolor = (float) locallab.spots.at(sp).blurcolde;
    int local_transit = locallab.spots.at(sp).transit;
    double radius = (double) locallab.spots.at(sp).radius;
    double sharradius = ((double) locallab.spots.at(sp).sharradius) / 100. ;
    double lcamount = ((double) locallab.spots.at(sp).lcamount) / 100. ;
    double sharblurr = ((double) locallab.spots.at(sp).sharblur) / 100. ;
    int local_sensisha = locallab.spots.at(sp).sensisha;
    int local_sharamount = locallab.spots.at(sp).sharamount;
    int local_shardamping = locallab.spots.at(sp).shardamping;
    int local_shariter = locallab.spots.at(sp).shariter;
    bool inverse = locallab.spots.at(sp).invers;
    bool curvacti = locallab.spots.at(sp).curvactiv;
    bool acti = locallab.spots.at(sp).activlum;
    bool cupas = false; // Provision
    int local_sensisf = locallab.spots.at(sp).sensisf;

    bool inverserad = false; // Provision
    bool inverseret = locallab.spots.at(sp).inversret;
    bool inversesha = locallab.spots.at(sp).inverssha;
    double strength = (double) locallab.spots.at(sp).strength;
    float str = (float)locallab.spots.at(sp).str;
//   lp.showmaskexpo = showmexp;

    lp.cir = circr;
    lp.actsp = acti;
    lp.xc = w * local_center_x;
    lp.yc = h * local_center_y;
    lp.xcbuf = w * local_center_xbuf;
    lp.ycbuf = h * local_center_ybuf;
    lp.yc = h * local_center_y;
    lp.lx = w * local_x;
    lp.ly = h * local_y;
    lp.lxL = w * local_xL;
    lp.lyT = h * local_yT;
    lp.chro = local_chroma;
    lp.struco = structcolor;
    lp.blendmacol = blendmaskcolor;
    lp.blendmaexp = blendmaskexpo;
    lp.struexp = structexpo;
    lp.blurexp = blurexpo;
    lp.blurcol = blurcolor;
    lp.sens = local_sensi;
    lp.sensh = local_sensih;
    lp.dehaze = local_dehaze;
    lp.senscb = local_sensicb;
    lp.cont = local_contrast;
    lp.ligh = local_lightness;
    lp.senssf = local_sensisf;
    lp.strng = strlight;
    lp.neig = neigh;

    if (lp.ligh >= -2.f && lp.ligh <= 2.f) {
        lp.ligh /= 5.f;
    }

    lp.trans = local_transit;
    lp.rad = radius;
    lp.stren = strength;
    lp.sensbn = local_sensibn;
    lp.sensexclu = local_sensiexclu;
    lp.senslc = local_sensilc;
    lp.lcamount = lcamount;
    lp.inv = inverse;
    lp.curvact = curvacti;
    lp.invrad = inverserad;
    lp.invret = inverseret;
    lp.invshar = inversesha;
    lp.str = str;
    lp.shrad = sharradius;
    lp.shblurr = sharblurr;
    lp.senssha = local_sensisha;
    lp.shamo = local_sharamount;
    lp.shdamp = local_shardamping;
    lp.shiter = local_shariter;
    lp.iterat = iterati;
    lp.dxx = w * local_dxx;
    lp.dyy = h * local_dyy;
    lp.thr = thre;
    lp.noiself = local_noiself;
    lp.noiseldetail = local_noiseldetail;
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

    for (int y = 0; y < 5; y++) {
        lp.mulloc[y] = multi[y];
    }

    lp.threshol = thresho;
    lp.chromacb = chromcbdl;
    lp.colorena = locallab.spots.at(sp).expcolor && llExpMask == 0; // Color & Light tool is deactivated if Exposure mask is visible
    lp.blurena = locallab.spots.at(sp).expblur;
    lp.tonemapena = locallab.spots.at(sp).exptonemap;
    lp.retiena = locallab.spots.at(sp).expreti;
    lp.sharpena = locallab.spots.at(sp).expsharp;
    lp.lcena = locallab.spots.at(sp).expcontrast;
    lp.sfena = locallab.spots.at(sp).expsoft;
    lp.cbdlena = locallab.spots.at(sp).expcbdl;
    lp.denoiena = locallab.spots.at(sp).expdenoi;
    lp.expvib = locallab.spots.at(sp).expvibrance;
    lp.sensv = local_sensiv;
    lp.past =  chromaPastel;
    lp.satur = chromaSatur;

    lp.exposena = locallab.spots.at(sp).expexpose && llColorMask == 0; // Exposure tool is deactivated if Color & Light mask is visible
    lp.cut_past = cupas;
    lp.blac = locallab.spots.at(sp).black;
    lp.shcomp = locallab.spots.at(sp).shcompr;
    lp.hlcomp = locallab.spots.at(sp).hlcompr;
    lp.hlcompthr = locallab.spots.at(sp).hlcomprthresh;
    lp.expcomp = locallab.spots.at(sp).expcomp / 100.;
    lp.sensex = local_sensiex;
    lp.strucc = local_struc;
    lp.war = local_warm;
}



static void calcTransitionrect(const float lox, const float loy, const float ach, const local_params& lp, int &zone, float &localFactor)
{
    zone = 0;

    if (lox >= lp.xc && lox < (lp.xc + lp.lx) && loy >= lp.yc && loy < lp.yc + lp.ly) {
        if (lox < (lp.xc + lp.lx * ach)  && loy < (lp.yc + lp.ly * ach)) {
            zone = 2;
        } else {
            zone = 1;
            localFactor = calcLocalFactorrect(lox, loy, lp.xc, lp.lx, lp.yc, lp.ly, ach);
        }

    } else if (lox >= lp.xc && lox < lp.xc + lp.lx && loy < lp.yc && loy > lp.yc - lp.lyT) {
        if (lox < (lp.xc + lp.lx * ach) && loy > (lp.yc - lp.lyT * ach)) {
            zone = 2;
        } else {
            zone = 1;
            localFactor = calcLocalFactorrect(lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
        }


    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy <= lp.yc && loy > lp.yc - lp.lyT) {
        if (lox > (lp.xc - lp.lxL * ach) && loy > (lp.yc - lp.lyT * ach)) {
            zone = 2;
        } else {
            zone = 1;
            localFactor = calcLocalFactorrect(lox, loy, lp.xc, lp.lxL, lp.yc, lp.lyT, ach);
        }

    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy > lp.yc && loy < lp.yc + lp.ly) {
        if (lox > (lp.xc - lp.lxL * ach) && loy < (lp.yc + lp.ly * ach)) {
            zone = 2;
        } else {
            zone = 1;
            localFactor = calcLocalFactorrect(lox, loy, lp.xc, lp.lxL, lp.yc, lp.ly, ach);
        }

    }

}



static void calcTransition(const float lox, const float loy, const float ach, const local_params& lp, int &zone, float &localFactor)
{
    // returns the zone (0 = outside selection, 1 = transition zone between outside and inside selection, 2 = inside selection)
    // and a factor to calculate the transition in case zone == 1

    zone = 0;

    if (lox >= lp.xc && lox < (lp.xc + lp.lx) && loy >= lp.yc && loy < lp.yc + lp.ly) {
        float zoneVal = SQR((lox - lp.xc) / (ach * lp.lx)) + SQR((loy - lp.yc) / (ach * lp.ly));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR((lox - lp.xc) / (lp.lx)) + SQR((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor(lox, loy, lp.xc, lp.lx, lp.yc, lp.ly, ach);
            }
        }
    } else if (lox >= lp.xc && lox < lp.xc + lp.lx && loy < lp.yc && loy > lp.yc - lp.lyT) {
        float zoneVal = SQR((lox - lp.xc) / (ach * lp.lx)) + SQR((loy - lp.yc) / (ach * lp.lyT));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR((lox - lp.xc) / (lp.lx)) + SQR((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor(lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
            }
        }
    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy <= lp.yc && loy > lp.yc - lp.lyT) {
        float zoneVal = SQR((lox - lp.xc) / (ach * lp.lxL)) + SQR((loy - lp.yc) / (ach * lp.lyT));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR((lox - lp.xc) / (lp.lxL)) + SQR((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor(lox, loy, lp.xc, lp.lxL, lp.yc, lp.lyT, ach);
            }
        }
    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy > lp.yc && loy < lp.yc + lp.ly) {
        float zoneVal = SQR((lox - lp.xc) / (ach * lp.lxL)) + SQR((loy - lp.yc) / (ach * lp.ly));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR((lox - lp.xc) / (lp.lxL)) + SQR((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor(lox, loy, lp.xc, lp.lxL, lp.yc, lp.ly, ach);
            }
        }
    }
}

/*
static void calcTransitioninv (const float lox, const float loy, const float ach, const local_params& lp, int &zone, float &localFactor)
{
    // returns the zone (0 = outside selection, 1 = zone between outside and inside selection, 2 = inside selection with transition)
    // and a factor to calculate the transition in case zone == 2

    zone = 0;

    if (lox >= lp.xc && lox < (lp.xc + lp.lx) && loy >= lp.yc && loy < lp.yc + lp.ly) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lx)) + SQR ((loy - lp.yc) / (ach * lp.ly));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lx)) + SQR ((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;
        }

        if (zone == 2) {
            localFactor = calcLocalFactorinv (lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
        }

    } else if (lox >= lp.xc && lox < lp.xc + lp.lx && loy < lp.yc && loy > lp.yc - lp.lyT) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lx)) + SQR ((loy - lp.yc) / (ach * lp.lyT));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lx)) + SQR ((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;
        }

        if (zone == 2) {
            localFactor = calcLocalFactorinv (lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
        }

    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy <= lp.yc && loy > lp.yc - lp.lyT) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lxL)) + SQR ((loy - lp.yc) / (ach * lp.lyT));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lxL)) + SQR ((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;
        }

        if (zone == 2) {
            localFactor = calcLocalFactorinv (lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
        }

    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy > lp.yc && loy < lp.yc + lp.ly) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lxL)) + SQR ((loy - lp.yc) / (ach * lp.ly));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lxL)) + SQR ((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;
        }

        if (zone == 2) {
            localFactor = calcLocalFactorinv (lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
        }

    }
}
*/

void ImProcFunctions::ciecamloc_02float(int sp, LabImage* lab, LabImage* dest)
{
    //be carefull quasi duplicate with branch cat02wb
    BENCHFUN
#ifdef _DEBUG
    MyTime t1e, t2e;
    t1e.set();
#endif

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


    xwd = 100.f * Xwout;
    zwd = 100.f * Zwout;
    ywd = 100.f;

    xws = 100.f * Xwsc;
    zws = 100.f * Zwsc;
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
    xw = 100.0f * Xw;
    yw = 100.0f * Yw;
    zw = 100.0f * Zw;
    float xw1 = 0.f, yw1 = 0.f, zw1 = 0.f, xw2 = 0.f, yw2 = 0.f, zw2 = 0.f;
// free temp and green
    xw1 = xws;
    yw1 = yws;
    zw1 = zws;
    xw2 = xwd;
    yw2 = ywd;
    zw2 = zwd;

    float cz, wh, pfl;
    Ciecam02::initcam1float(yb, pilot, f, la, xw, yw, zw, n, d, nbb, ncb, cz, aw, wh, pfl, fl, c);
//   const float chr = 0.f;
    const float pow1 = pow_F(1.64f - pow_F(0.29f, n), 0.73f);
    float nj, nbbj, ncbj, czj, awj, flj;
    Ciecam02::initcam2float(yb2, pilotout, f2,  la2,  xw2,  yw2,  zw2, nj, dj, nbbj, ncbj, czj, awj, flj);
    const float reccmcz = 1.f / (c2 * czj);
    const float pow1n = pow_F(1.64f - pow_F(0.29f, nj), 0.73f);
//    const float QproFactor = (0.4f / c) * (aw + 4.0f) ;
    const bool LabPassOne = true;

#ifdef __SSE2__
    int bufferLength = ((width + 3) / 4) * 4; // bufferLength has to be a multiple of 4
#endif
#ifndef _DEBUG
    #pragma omp parallel
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
#ifndef _DEBUG
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
                x = (float)x1 / 655.35f;
                y = (float)y1 / 655.35f;
                z = (float)z1 / 655.35f;
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

                if (LabPassOne) {
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
                    float x, y, z;
                    x = xx * 655.35f;
                    y = yy * 655.35f;
                    z = zz * 655.35f;
                    float Ll, aa, bb;
                    //convert xyz=>lab
                    Color::XYZ2Lab(x,  y,  z, Ll, aa, bb);
                    dest->L[i][j] = Ll;
                    dest->a[i][j] = aa;
                    dest->b[i][j] = bb;

#endif
                }

                //    }
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

                dest->L[i][j] = Ll;
                dest->a[i][j] = aa;
                dest->b[i][j] = bb;
            }

#endif
        }

    }
#ifdef _DEBUG

    if (settings->verbose) {
        t2e.set();
        printf("CAT02 local 02 performed in %d usec:\n", t2e.etime(t1e));
    }

#endif
}




void ImProcFunctions::vibrancelocal(int sp, int bfw, int bfh, LabImage* lab,  LabImage* dest, bool & localskutili, LUTf & sklocalcurve)
{
    if (!((bool)params->locallab.spots.at(sp).expvibrance)) {
        return;
    }

    const int width = bfw;
    const int height = bfh;

#ifdef _DEBUG
    MyTime t1e, t2e;
    t1e.set();
    int negat = 0, moreRGB = 0, negsat = 0, moresat = 0;
#endif

    const float chromaPastel = float (params->locallab.spots.at(sp).pastels)   / 100.0f;
    const float chromaSatur  = float (params->locallab.spots.at(sp).saturated) / 100.0f;
    const float p00 = 0.07f;
    const float limitpastelsatur = (static_cast<float>(params->locallab.spots.at(sp).psthreshold.getTopLeft())    / 100.0f) * (1.0f - p00) + p00;
    const float maxdp = (limitpastelsatur - p00) / 4.0f;
    const float maxds = (1.0 - limitpastelsatur) / 4.0f;
    const float p0 = p00 + maxdp;
    const float p1 = p00 + 2.0f * maxdp;
    const float p2 = p00 + 3.0f * maxdp;
    const float s0 = limitpastelsatur + maxds;
    const float s1 = limitpastelsatur + 2.0f * maxds;
    const float s2 = limitpastelsatur + 3.0f * maxds;
    const float transitionweighting = static_cast<float>(params->locallab.spots.at(sp).psthreshold.getBottomLeft()) / 100.0f;
    float chromamean = 0.0f;

    if (chromaPastel != chromaSatur) {
        //if sliders pastels and saturated are different: transition with a double linear interpolation: between p2 and limitpastelsatur, and between limitpastelsatur and s0
        //modify the "mean" point in function of double threshold  => differential transition
        chromamean = maxdp * (chromaSatur - chromaPastel) / (s0 - p2) + chromaPastel;

        // move chromaMean up or down depending on transitionCtrl
        if (transitionweighting > 0.0f) {
            chromamean = (chromaSatur - chromamean) * transitionweighting + chromamean;
        } else if (transitionweighting < 0.0f) {
            chromamean = (chromamean - chromaPastel)  * transitionweighting + chromamean;
        }
    }

    const float chromaPastel_a = (chromaPastel - chromamean) / (p2 - limitpastelsatur);
    const float chromaPastel_b = chromaPastel - chromaPastel_a * p2;

    const float chromaSatur_a = (chromaSatur - chromamean) / (s0 - limitpastelsatur);
    const float chromaSatur_b = chromaSatur - chromaSatur_a * s0;

    const float dhue = 0.15f; //hue transition
    const float dchr = 20.0f; //chroma transition
    const float skbeg = -0.05f; //begin hue skin
    const float skend = 1.60f; //end hue skin
    const float xx = 0.5f; //soft : between 0.3 and 1.0
    const float ask = 65535.0f / (skend - skbeg);
    const float bsk = -skbeg * ask;


    const bool highlight = params->toneCurve.hrenabled;//Get the value if "highlight reconstruction" is activated
    const bool protectskins = params->locallab.spots.at(sp).protectskins;
    const bool avoidcolorshift = params->locallab.spots.at(sp).avoidcolorshift;

    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);
    //inverse matrix user select
    const double wip[3][3] = {
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };


#ifdef _DEBUG
    MunsellDebugInfo* MunsDebugInfo = nullptr;

    if (avoidcolorshift) {
        MunsDebugInfo = new MunsellDebugInfo();
    }

    #pragma omp parallel default(shared) firstprivate(lab, dest, MunsDebugInfo) reduction(+: negat, moreRGB, negsat, moresat) if (multiThread)
#else
    #pragma omp parallel default(shared) if (multiThread)
#endif
    {

        float sathue[5], sathue2[4]; // adjust sat in function of hue


#ifdef _OPENMP

        if (settings->verbose && omp_get_thread_num() == 0) {
#else

        if (settings->verbose) {
#endif
            printf("vibrance:  p0=%1.2f  p1=%1.2f  p2=%1.2f  s0=%1.2f s1=%1.2f s2=%1.2f\n", p0, p1, p2, s0, s1, s2);
            printf("           pastel=%f   satur=%f   limit= %1.2f   chromamean=%0.5f\n", 1.0f + chromaPastel, 1.0f + chromaSatur, limitpastelsatur, chromamean);
        }

        #pragma omp for schedule(dynamic, 16)

        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++) {
                float LL = lab->L[i][j] / 327.68f;
                float CC = sqrt(SQR(lab->a[i][j]) + SQR(lab->b[i][j])) / 327.68f;
                float HH = xatan2f(lab->b[i][j], lab->a[i][j]);

                float satredu = 1.0f; //reduct sat in function of skin

                if (protectskins) {
                    Color::SkinSat(LL, HH, CC, satredu); // for skin colors
                }

                // here we work on Chromaticity and Hue
                // variation of Chromaticity  ==> saturation via RGB
                // Munsell correction, then conversion to Lab
                float Lprov = LL;
                float Chprov = CC;
                float R, G, B;
                float2 sincosval;

                if (CC == 0.0f) {
                    sincosval.y = 1.f;
                    sincosval.x = 0.0f;
                } else {
                    sincosval.y = lab->a[i][j] / (CC * 327.68f);
                    sincosval.x = lab->b[i][j] / (CC * 327.68f);
                }

#ifdef _DEBUG
                bool neg = false;
                bool more_rgb = false;
                //gamut control : Lab values are in gamut
                Color::gamutLchonly(HH, sincosval, Lprov, Chprov, R, G, B, wip, highlight, 0.15f, 0.98f, neg, more_rgb);

                if (neg) {
                    negat++;
                }

                if (more_rgb) {
                    moreRGB++;
                }

#else
                //gamut control : Lab values are in gamut
                Color::gamutLchonly(HH, sincosval, Lprov, Chprov, R, G, B, wip, highlight, 0.15f, 0.98f);
#endif

                if (Chprov > 6.0f) {
                    const float saturation = SAT(R, G, B);

                    if (saturation > 0.0f) {
                        if (satredu != 1.0f) {
                            // for skin, no differentiation
                            sathue [0] = sathue [1] = sathue [2] = sathue [3] = sathue[4] = 1.0f;
                            sathue2[0] = sathue2[1] = sathue2[2] = sathue2[3]          = 1.0f;
                        } else {
                            //double pyramid: LL and HH
                            //I try to take into account: Munsell response (human vision) and Gamut..(less response for red): preferably using Prophoto or WideGamut
                            //blue: -1.80 -3.14  green = 2.1 3.14   green-yellow=1.4 2.1  red:0 1.4  blue-purple:-0.7  -1.4   purple: 0 -0.7
                            //these values allow a better and differential response
                            if (LL < 20.0f) { //more for blue-purple, blue and red modulate
                                if (/*HH> -3.1415f &&*/ HH < -1.5f) {
                                    sathue[0] = 1.3f;    //blue
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.05f;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f) {
                                    sathue[0] = 1.6f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.15f;
                                    sathue2[2] = 1.1f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=1.1f;sathue[1]=1.1f;sathue[2]=1.1f;sathue[3]=1.0f ;sathue[4]=0.4f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f) {
                                    sathue[0] = 1.3f;    //red   0.8 0.7
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f) {
                                    sathue[0] = 1.0f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.4f;    //green
                                    sathue[1] = 1.3f;
                                    sathue[2] = 1.2f;
                                    sathue[3] = 1.15f;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.15f;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                }
                            } else if (LL < 50.0f) { //more for blue and green, less for red and green-yellow
                                if (/*HH> -3.1415f &&*/ HH < -1.5f) {
                                    sathue[0] = 1.5f;    //blue
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f) {
                                    sathue[0] = 1.3f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.05f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f;sathue[3]=0.8f ;sathue[4]=0.4f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 0.8f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f) {
                                    sathue[0] = 1.1f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.1f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 0.9f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.7f ;
                                    sathue2[3] = 0.6f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.5f;    //green
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                }

                            } else if (LL < 80.0f) { //more for green, less for red and green-yellow
                                if (/*HH> -3.1415f &&*/ HH < -1.5f) {
                                    sathue[0] = 1.3f;    //blue
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.15f;
                                    sathue[3] = 1.1f ;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.1f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f) {
                                    sathue[0] = 1.3f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.15f;
                                    sathue[3] = 1.1f ;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.1f ;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f ;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f ;sathue[3]=0.8f ;sathue[4]=0.3f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f ;
                                    sathue[3] = 0.8f ;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 0.8f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f) {
                                    sathue[0] = 1.3f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f ;
                                    sathue[3] = 1.05f;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 0.9f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.7f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.6f;    //green - even with Prophoto green are too "little"  1.5 1.3
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f ;
                                    sathue[3] = 1.25f;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.25f;
                                    sathue2[1] = 1.2f ;
                                    sathue2[2] = 1.15f;
                                    sathue2[3] = 1.05f;
                                }
                            } else { /*if (LL>=80.0f)*/ //more for green-yellow, less for red and purple
                                if (/*HH> -3.1415f &&*/ HH < -1.5f) {
                                    sathue[0] = 1.0f;    //blue
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f) {
                                    sathue[0] = 1.0f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 0.9f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 0.9f;
                                    sathue2[1] = 0.9f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f;sathue[3]=0.8f;sathue[4]=0.2f;sathue2[0]=0.8f;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f) {
                                    sathue[0] = 1.6f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.5f;
                                    sathue[2] = 1.4f;
                                    sathue[3] = 1.2f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 1.1f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.4f;    //green
                                    sathue[1] = 1.3f;
                                    sathue[2] = 1.2f;
                                    sathue[3] = 1.1f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 1.1f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                }
                            }
                        }

                        float chmodpastel = 0.f, chmodsat = 0.f;
                        // variables to improve transitions
                        float pa, pb;// transition = pa*saturation + pb
                        float chl00 = chromaPastel * satredu * sathue[4];
                        float chl0  = chromaPastel * satredu * sathue[0];
                        float chl1  = chromaPastel * satredu * sathue[1];
                        float chl2  = chromaPastel * satredu * sathue[2];
                        float chl3  = chromaPastel * satredu * sathue[3];
                        float chs0  = chromaSatur * satredu * sathue2[0];
                        float chs1  = chromaSatur * satredu * sathue2[1];
                        float chs2  = chromaSatur * satredu * sathue2[2];
                        float chs3  = chromaSatur * satredu * sathue2[3];
                        float s3    = 1.0f;

                        // We handle only positive values here ;  improve transitions
                        if (saturation < p00) {
                            chmodpastel = chl00 ;    //neutral tones
                        } else if (saturation < p0)               {
                            pa = (chl00 - chl0) / (p00 - p0);
                            pb = chl00 - pa * p00;
                            chmodpastel = pa * saturation + pb;
                        } else if (saturation < p1)                {
                            pa = (chl0 - chl1) / (p0 - p1);
                            pb = chl0 - pa * p0;
                            chmodpastel = pa * saturation + pb;
                        } else if (saturation < p2)                {
                            pa = (chl1 - chl2) / (p1 - p2);
                            pb = chl1 - pa * p1;
                            chmodpastel = pa * saturation + pb;
                        } else if (saturation < limitpastelsatur)  {
                            pa = (chl2 - chl3) / (p2 - limitpastelsatur);
                            pb = chl2 - pa * p2;
                            chmodpastel = pa * saturation + pb;
                        } else if (saturation < s0)                {
                            pa = (chl3 - chs0) / (limitpastelsatur - s0) ;
                            pb = chl3 - pa * limitpastelsatur;
                            chmodsat    = pa * saturation + pb;
                        } else if (saturation < s1)                {
                            pa = (chs0 - chs1) / (s0 - s1);
                            pb = chs0 - pa * s0;
                            chmodsat    = pa * saturation + pb;
                        } else if (saturation < s2)                {
                            pa = (chs1 - chs2) / (s1 - s2);
                            pb = chs1 - pa * s1;
                            chmodsat    = pa * saturation + pb;
                        } else                                     {
                            pa = (chs2 - chs3) / (s2 - s3);
                            pb = chs2 - pa * s2;
                            chmodsat    = pa * saturation + pb;
                        }

                        if (chromaPastel != chromaSatur) {

                            // Pastels
                            if (saturation > p2 && saturation < limitpastelsatur) {
                                float newchromaPastel = chromaPastel_a * saturation + chromaPastel_b;
                                chmodpastel = newchromaPastel * satredu * sathue[3];
                            }

                            // Saturated
                            if (saturation < s0 && saturation >= limitpastelsatur) {
                                float newchromaSatur = chromaSatur_a * saturation + chromaSatur_b;
                                chmodsat = newchromaSatur * satredu * sathue2[0];
                            }
                        }// end transition

                        if (saturation <= limitpastelsatur) {
                            if (chmodpastel >  2.0f) {
                                chmodpastel = 2.0f;    //avoid too big values
                            } else if (chmodpastel < -0.93f) {
                                chmodpastel = -0.93f;    //avoid negative values
                            }

                            Chprov *= (1.0f + chmodpastel);

                            if (Chprov < 6.0f) {
                                Chprov = 6.0f;
                            }
                        } else { //if (saturation > limitpastelsatur)
                            if (chmodsat >  1.8f) {
                                chmodsat = 1.8f;    //saturated
                            } else if (chmodsat < -0.93f) {
                                chmodsat = -0.93f;
                            }

                            Chprov *= 1.0f + chmodsat;

                            if (Chprov < 6.0f) {
                                Chprov = 6.0f;
                            }
                        }
                    }
                }

                bool hhModified = false;

                // Vibrance's Skin curve
                if (sklocalcurve  && localskutili) {
                    if (HH > skbeg && HH < skend) {
                        if (Chprov < 60.0f) { //skin hue  : todo ==> transition
                            float HHsk = ask * HH + bsk;
                            float Hn = (sklocalcurve[HHsk] - bsk) / ask;
                            float Hc = (Hn * xx + HH * (1.0f - xx));
                            HH = Hc;
                            hhModified = true;
                        } else if (Chprov < (60.0f + dchr)) { //transition chroma
                            float HHsk = ask * HH + bsk;
                            float Hn = (sklocalcurve[HHsk] - bsk) / ask;
                            float Hc = (Hn * xx + HH * (1.0f - xx));
                            float aa = (HH - Hc) / dchr ;
                            float bb = HH - (60.0f + dchr) * aa;
                            HH = aa * Chprov + bb;
                            hhModified = true;
                        }
                    }
                    //transition hue
                    else if (HH > (skbeg - dhue) && HH <= skbeg && Chprov < (60.0f + dchr * 0.5f)) {
                        float HHsk = ask * skbeg + bsk;
                        float Hn = (sklocalcurve[HHsk] - bsk) / ask;
                        float Hcc = (Hn * xx + skbeg * (1.0f - xx));
                        float adh = (Hcc - (skbeg - dhue)) / (dhue);
                        float bdh = Hcc - adh * skbeg;
                        HH = adh * HH + bdh;
                        hhModified = true;
                    } else if (HH >= skend && HH < (skend + dhue) && Chprov < (60.0f + dchr * 0.5f)) {
                        float HHsk = ask * skend + bsk;
                        float Hn = (sklocalcurve[HHsk] - bsk) / ask;
                        float Hcc = (Hn * xx + skend * (1.0f - xx));
                        float adh = (skend + dhue - Hcc) / (dhue);
                        float bdh = Hcc - adh * skend;
                        HH = adh * HH + bdh;
                        hhModified = true;
                    }
                } // end skin hue

                //Munsell correction
                if (!avoidcolorshift && hhModified) {
                    sincosval = xsincosf(HH);
                }

                float aprovn, bprovn;
                bool inGamut;

                do {
                    inGamut = true;

                    if (avoidcolorshift) {
                        float correctionHue = 0.0f;
                        float correctlum = 0.0f;

#ifdef _DEBUG
                        Color::AllMunsellLch(false, Lprov, Lprov, HH, Chprov, CC, correctionHue, correctlum, MunsDebugInfo);
#else
                        Color::AllMunsellLch(false, Lprov, Lprov, HH, Chprov, CC, correctionHue, correctlum);
#endif

                        if (correctionHue != 0.f || hhModified) {
                            sincosval = xsincosf(HH + correctionHue);
                            hhModified = false;
                        }
                    }

                    aprovn = Chprov * sincosval.y;
                    bprovn = Chprov * sincosval.x;

                    float fyy = (Color::c1By116 * Lprov) + Color::c16By116;
                    float fxx = (0.002f * aprovn) + fyy;
                    float fzz = fyy - (0.005f * bprovn);
                    float xx_ = 65535.f * Color::f2xyz(fxx) * Color::D50x;
                    //  float yy_ = 65535.0f * Color::f2xyz(fyy);
                    float zz_ = 65535.f * Color::f2xyz(fzz) * Color::D50z;
                    float yy_ = 65535.f * ((Lprov > Color::epskap) ? fyy * fyy*fyy : Lprov / Color::kappa);

                    Color::xyz2rgb(xx_, yy_, zz_, R, G, B, wip);

                    if (R < 0.0f || G < 0.0f || B < 0.0f) {
#ifdef _DEBUG
                        negsat++;
#endif
                        Chprov *= 0.98f;
                        inGamut = false;
                    }

                    // if "highlight reconstruction" enabled don't control Gamut for highlights
                    if ((!highlight) && (R > 65535.0f || G > 65535.0f || B > 65535.0f)) {
#ifdef _DEBUG
                        moresat++;
#endif
                        Chprov *= 0.98f;
                        inGamut = false;
                    }
                } while (!inGamut);

                //put new values in Lab
                dest->L[i][j] = Lprov * 327.68f;
                dest->a[i][j] = aprovn * 327.68f;
                dest->b[i][j] = bprovn * 327.68f;
            }
    } // end of parallelization

#ifdef _DEBUG
    t2e.set();

    if (settings->verbose) {
        printf("Vibrance local (performed in %d usec):\n", t2e.etime(t1e));
        printf("   Gamut: G1negat=%iiter G165535=%iiter G2negsat=%iiter G265535=%iiter\n", negat, moreRGB, negsat, moresat);

        if (MunsDebugInfo) {
            printf("   Munsell chrominance: MaxBP=%1.2frad  MaxRY=%1.2frad  MaxGY=%1.2frad  MaxRP=%1.2frad  depass=%u\n", MunsDebugInfo->maxdhue[0], MunsDebugInfo->maxdhue[1], MunsDebugInfo->maxdhue[2], MunsDebugInfo->maxdhue[3], MunsDebugInfo->depass);
        }
    }

    if (MunsDebugInfo) {
        delete MunsDebugInfo;
    }

#endif

}

void ImProcFunctions::exlabLocal(const local_params& lp, int bfh, int bfw, LabImage* bufexporig, LabImage* lab,  LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve)
{
    //exposure local

    float maxran = 65536.f; //65536
    const float exp_scale = pow(2.0, lp.expcomp); //lp.expcomp
    const float comp = (max(0.0, lp.expcomp) + 1.0) * lp.hlcomp / 100.0;
    const float shoulder = ((maxran / max(1.0f, exp_scale)) * (lp.hlcompthr / 200.0)) + 0.1;
    const float hlrange = maxran - shoulder;


#define TSE 112

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        char *buffer;

        buffer = (char *) malloc(3 * sizeof(float) * TSE * TSE + 20 * 64 + 63);
        char *data;
        data = (char*)((uintptr_t (buffer) + uintptr_t (63)) / 64 * 64);

        float *Ltemp = (float (*))data;
        float *atemp = (float (*))((char*)Ltemp + sizeof(float) * TSE * TSE + 4 * 64);
        float *btemp = (float (*))((char*)atemp + sizeof(float) * TSE * TSE + 8 * 64);
        int istart;
        int jstart;
        int tW;
        int tH;

#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif

        for (int ii = 0; ii < bfh; ii += TSE)
            for (int jj = 0; jj < bfw; jj += TSE) {

                istart = ii;
                jstart = jj;
                tH = min(ii + TSE, bfh);
                tW = min(jj + TSE, bfw);


                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                        Ltemp[ti * TSE + tj] = bufexporig->L[i][j];
                        atemp[ti * TSE + tj] = bufexporig->a[i][j];
                        btemp[ti * TSE + tj] = bufexporig->b[i][j];;
                    }
                }


                //    float niv = maxran;

                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        float L = Ltemp[ti * TSE + tj];
                        float tonefactor = (2 * L < MAXVALF ? hltonecurve[2 * L] : CurveFactory::hlcurve(exp_scale, comp, hlrange, 2 * L)); // niv));
                        Ltemp[ti * TSE + tj] = L * tonefactor;
                    }
                }

                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        float L = Ltemp[ti * TSE + tj];
                        //shadow tone curve
                        float Y = L;
                        float tonefactor = shtonecurve[2 * Y];
                        Ltemp[ti * TSE + tj] = Ltemp[ti * TSE + tj] * tonefactor;
                    }
                }

                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        Ltemp[ti * TSE + tj] = tonecurve[Ltemp[ti * TSE + tj] ];

                    }
                }

                /*
                                if (lp.chro != 0) {
                                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                                            float satby100 = lp.chro / 100.f;
                                            float a = atemp[ti * TSE + tj];
                                            float b = btemp[ti * TSE + tj];

                                            atemp[ti * TSE + tj] = a * (1.f + satby100);
                                            btemp[ti * TSE + tj] = b * (1.f + satby100);
                                        }
                                    }
                                }
                */

                bool vasy = true;

                if (vasy) {
                    // ready, fill lab
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                            lab->L[i][j] = Ltemp[ti * TSE + tj];
                            lab->a[i][j] = atemp[ti * TSE + tj];
                            lab->b[i][j] = btemp[ti * TSE + tj];
                        }
                    }
                }
            }

        free(buffer);


    }


}


void ImProcFunctions::addGaNoise(LabImage *lab, LabImage *dst, const float mean, const float variance, const int sk)
{
//   BENCHFUN
//Box-Muller method.
// add luma noise to image

    srand(1);

    const float variaFactor = SQR(variance) / sk;
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

                float u1f = u1 * randFactor;
                float u2f = u2 * randFactor;

                float2 sincosval = xsincosf(2.f * rtengine::RT_PI * u2f);
                float factor = sqrtf(-2.f * xlogf(u1f));
                z0 = factor * sincosval.y;
                z1 = factor * sincosval.x;

                dst->L[y][x] = LIM(lab->L[y][x] + mean + varia * z0, 0.f, 32768.f);

            }
        }
    }
}


void ImProcFunctions::DeNoise_Local(int call, const struct local_params& lp,  int levred, float hueplus, float huemoins, float hueref, float dhueden, LabImage* original, LabImage* transformed, LabImage &tmp1, int cx, int cy, int sk)
{
    // local denoise and impulse
    //simple algo , perhaps we can improve as the others, but noise is here and not good for hue detection
    // BENCHFUN
    const float ach = (float)lp.trans / 100.f;
    constexpr float delhu = 0.1f; //between 0.05 and 0.2


    float factnoise1 = 1.f + (lp.noisecf) / 500.f;
    float factnoise2 = 1.f + (lp.noisecc) / 500.f;

    constexpr float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    constexpr float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;
    /*
        constexpr float pb = 4.f;
        constexpr float pa = (1.f - pb) / 40.f;
    */
    int GW = transformed->W;
    int GH = transformed->H;

    LabImage *origblur = nullptr;

    origblur = new LabImage(GW, GH);

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
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
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
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                }

                continue;
            }

#ifdef __SSE2__
            int i = 0;

            for (; i < transformed->W - 3; i += 4) {
                vfloat av = LVFU(origblur->a[y][i]);
                vfloat bv = LVFU(origblur->b[y][i]);
                STVF(atan2Buffer[i], xatan2f(bv, av));
                STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f(origblur->b[y][i], origblur->a[y][i]);
                sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
            }

#endif

            for (int x = 0, lox = cx + x; x < transformed->W; x++, lox++) {
                int zone = 0;
                // int lox = cx + x;
                int begx = int (lp.xc - lp.lxL);
                int begy = int (lp.yc - lp.lyT);

                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else if (lp.shapmet == 1) {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }


                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                    continue;
                }

#ifdef __SSE2__
                const float rhue = atan2Buffer[x];
                //    const float rchro = sqrtBuffer[x];
#else
                const float rhue = xatan2f(origblur->b[y][x], origblur->a[y][x]);
                //    const float rchro = sqrt(SQR(origblur->b[y][x]) + SQR(origblur->a[y][x])) / 327.68f;
#endif

                //    float rL = original->L[y][x] / 327.68f;

                float khu = 0.f;
                //   bool kzon = false;

                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if (levred == 7 && lp.sensden < 90.f) { // after 90 plein effect
                    //hue detection
                    if ((hueref + dhueden) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }


                    } else if ((hueref + dhueden) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }


                    }

                    if ((hueref - dhueden) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }


                    } else if ((hueref - dhueden) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }


                    }
                } else {
                    khu = 1.f;
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

                        if (call == 2  /*|| call == 1  || call == 3 */) { //simpleprocess
                            difL = tmp1.L[loy - begy][lox - begx] - original->L[y][x];
                            difa = tmp1.a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1.b[loy - begy][lox - begx] - original->b[y][x];
                        } else  { //dcrop
                            difL = tmp1.L[y][x] - original->L[y][x];
                            difa = tmp1.a[y][x] - original->a[y][x];
                            difb = tmp1.b[y][x] - original->b[y][x];

                        }

                        difL *= factorx * khu;
                        difa *= factorx * khu;
                        difb *= factorx * khu;
                        transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                        transformed->a[y][x] = CLIPC((original->a[y][x] + difa) * factnoise1 * factnoise2);
                        transformed->b[y][x] = CLIPC((original->b[y][x] + difb) * factnoise1 * factnoise2) ;
                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        float difL, difa, difb;

                        if (call == 2 /*|| call == 1 || call == 3 */) { //simpleprocess
                            difL = tmp1.L[loy - begy][lox - begx] - original->L[y][x];
                            difa = tmp1.a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1.b[loy - begy][lox - begx] - original->b[y][x];
                        } else  { //dcrop
                            difL = tmp1.L[y][x] - original->L[y][x];
                            difa = tmp1.a[y][x] - original->a[y][x];
                            difb = tmp1.b[y][x] - original->b[y][x];

                        }

                        difL *= khu;
                        difa *= khu;
                        difb *= khu;

                        transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                        transformed->a[y][x] = CLIPC((original->a[y][x] + difa) * factnoise1 * factnoise2);
                        transformed->b[y][x] = CLIPC((original->b[y][x] + difb) * factnoise1 * factnoise2);
                    }
                }

            }
        }
    }
    delete origblur;

}

void ImProcFunctions::TM_Local(float moddE, float powdE, LabImage * tmp1, float **buflight, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
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

    float refa = chromaref * cos(hueref);
    float refb = chromaref * sin(hueref);

    //  const float moddE = 2.f;

    const float ahu = 1.f / (2.8f * lp.senstm - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.senstm;

    const float alum = 1.f / (lp.senstm - 100.f);
    const float blum = 1.f - alum * lp.senstm;

    int GW = transformed->W;
    int GH = transformed->H;

    LabImage *origblur = nullptr;

    origblur = new LabImage(GW, GH);

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
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
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
                vfloat av = LVFU(origblur->a[y][i]);
                vfloat bv = LVFU(origblur->b[y][i]);
                STVF(atan2Buffer[i], xatan2f(bv, av));
                STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f(origblur->b[y][i], origblur->a[y][i]);
                sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
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

                    if (lp.shapmet == 0) {
                        calcTransition(lox, loy, ach, lp, zone, localFactor);
                    } else if (lp.shapmet == 1) {
                        calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                    }


                    if (zone == 0) {
                        continue;
                    }


#ifdef __SSE2__
                    float rhue = atan2Buffer[x];
                    float rchro = sqrtBuffer[x];
#else
                    float rhue = xatan2f(origblur->b[y][x], origblur->a[y][x]);
                    float rchro = sqrt(SQR(origblur->b[y][x]) + SQR(origblur->a[y][x])) / 327.68f;
#endif
                    // int zone;
                    float rL = origblur->L[y][x] / 327.68f;
                    float dE = sqrt(SQR(refa - origblur->a[y][x] / 327.68f) + SQR(refb - origblur->b[y][x] / 327.68f) + SQR(lumaref - rL));

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


                    //prepare shape detection
                    float khu = 0.f;
                    float kch = 1.f;
                    float kchchro = 1.f;
                    float fach = 1.f;
                    float falu = 1.f;
                    float deltachro = fabs(rchro - chromaref);
                    float deltahue = fabs(rhue - hueref);

                    if (deltahue > rtengine::RT_PI) {
                        deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                    }

                    float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280
                    float deltaL = fabs(lumaref - rL);  //between 0 and 100

                    /*
                                        //kch to modulate action with chroma
                                        if (deltachro < 160.f * SQR(lp.senstm / 100.f)) {
                                            kch = 1.f;
                                        } else {
                                            float ck = 160.f * SQR(lp.senstm / 100.f);
                                            float ak = 1.f / (ck - 160.f);
                                            float bk = -160.f * ak;
                                            kch = ak * deltachro + bk;
                                        }

                                        if (lp.senstm < 40.f) {
                                            kch = pow(kch, pa * lp.senstm + pb);    //increase under 40
                                        }
                    */
                    if (deltachro < 160.f * SQR(lp.senstm / 100.f)) {
                        kch = kchchro = 1.f;
                    } else {
                        float ck = 160.f * SQR(lp.senstm / 100.f);
                        float ak = 1.f / (ck - 160.f);
                        float bk = -160.f * ak;
                        // kch = ak * deltachro + bk;
                        kch  = kchchro = ak * deltachro + bk;

                    }

                    float kkch = 1.f;
                    kkch = pa * lp.senstm + pb;
                    float kch0 = kch;

                    if (lp.senstm < 40.f) {
                        kch = kchchro = pow(kch0, kkch);    //increase under 40
                        float dEsens = moddE * lp.senstm;
                        float kdE = 1.f;

                        if (dE > dEsens) {
                            kdE = 1.f;
                        } else {
                            kdE = SQR(SQR(dE / dEsens));
                        }

                        if (deltahue < 0.3f && settings->detectshape == true) {
                            kchchro = kch = pow(kch0, (1.f + kdE * moddE * (3.f - 10.f * deltahue)) * kkch);
                        }
                    }

                    float dEsensall = lp.senstm;
                    float kD = 1.f;

                    if (settings->detectshape == true) {
                        if (dE < dEsensall) {
                            kD = 1.f;
                        } else {
                            kD = pow(dEsensall / dE, powdE);
                        }
                    }

                    // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                    if (lp.senstm < 100.f) { //to try...
                        //hue detection
                        if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                            if (rhue >= hueplus - delhu)  {
                                realcligh = apluscligh * rhue + bpluscligh;

                                khu  = apl * rhue + bpl;
                            } else if (rhue < huemoins + delhu)  {
                                realcligh = amoinscligh * rhue + bmoinscligh;

                                khu = amo * rhue + bmo;
                            } else {
                                khu = 1.f;
                                realcligh = cli;

                            }


                        } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
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

                        }

                        if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
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

                        } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
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

                        }

                        realcligh *= kD;

                        if (settings->detectshape == false) {

                            if (lp.senstm <= 20.f) {

                                if (deltaE <  2.8f * lp.senstm) {
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

                                if (deltaL <  lp.senstm) {
                                    falu = 1.f;
                                } else {
                                    falu = alum * deltaL + blum;
                                }
                            }
                        }

                    } else {
                    }

                    float fli = ((100.f + realcligh) / 100.f);//luma transition
                    float kcr = 100.f * lp.thr;
                    float falL = 1.f;

                    if (rchro < kcr && chromaref > kcr) { // reduce artifacts in grey tones near hue spot and improve algorithm
                        falL *= pow(rchro / kcr, lp.iterat / 10.f);
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

                            transformed->L[y][x] = CLIP(original->L[y][x] + difL * kch * fach);


                            transformed->a[y][x] = CLIPC(original->a[y][x] + difa * kch * fach * falu);//same as Luma
                            transformed->b[y][x] = CLIPC(original->b[y][x] + difb * kch * fach * falu);//same as Luma

                            break;
                        }

                        case 2: { // inside selection => full effect, no transition
                            float difL, difa, difb;

                            difL = tmp1->L[loy - begy][lox - begx] * fli * falL - original->L[y][x];
                            difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];

                            transformed->L[y][x] = CLIP(original->L[y][x] + difL * kch * fach);


                            transformed->a[y][x] = CLIPC(original->a[y][x] + difa * kch * fach * falu);//same as Luma
                            transformed->b[y][x] = CLIPC(original->b[y][x] + difb * kch * fach * falu);//same as Luma
                        }
                    }

                }
            }
        }
    }
    delete origblur;
}


void ImProcFunctions::BlurNoise_Local(float moddE, float powdE, int call, LabImage * tmp1, LabImage * tmp2, float ** buflight, float ** bufchro, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
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
//   const float moddE = 2.f;

    const float alum = 1.f / (lp.sensbn - 100.f);
    const float blum = 1.f - alum * lp.sensbn;

    int GW = transformed->W;
    int GH = transformed->H;
    float refa = chromaref * cos(hueref);
    float refb = chromaref * sin(hueref);

    LabImage *origblur = nullptr;

    origblur = new LabImage(GW, GH);

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
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < transformed->H; y++) {
            const int loy = cy + y;

            const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

            if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                for (int x = 0; x < transformed->W; x++) {
                    if (lp.blurmet == 0) {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    if (lp.blurmet == 2) {
                        transformed->L[y][x] = tmp2->L[y][x];
                    }
                }

                if (!lp.actsp) {
                    for (int x = 0; x < transformed->W; x++) {
                        if (lp.blurmet == 0) {
                            transformed->a[y][x] = original->a[y][x];
                            transformed->b[y][x] = original->b[y][x];
                        }

                        if (lp.blurmet == 2) {
                            transformed->a[y][x] = tmp2->a[y][x];
                            transformed->b[y][x] = tmp2->b[y][x];
                        }
                    }
                }

                continue;
            }

#ifdef __SSE2__
            int i = 0;

            for (; i < transformed->W - 3; i += 4) {
                vfloat av = LVFU(origblur->a[y][i]);
                vfloat bv = LVFU(origblur->b[y][i]);
                STVF(atan2Buffer[i], xatan2f(bv, av));
                STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f(origblur->b[y][i], origblur->a[y][i]);
                sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
            }

#endif



            for (int x = 0, lox = cx + x; x < transformed->W; x++, lox++) {
                int zone = 0;
                // int lox = cx + x;
                int begx = int (lp.xc - lp.lxL);
                int begy = int (lp.yc - lp.lyT);

                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else if (lp.shapmet == 1) {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }


                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    if (lp.blurmet == 0) {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    if (lp.blurmet == 2) {
                        transformed->L[y][x] = tmp2->L[y][x];
                    }

                    if (!lp.actsp) {
                        if (lp.blurmet == 0) {
                            transformed->a[y][x] = original->a[y][x];
                            transformed->b[y][x] = original->b[y][x];
                        }

                        if (lp.blurmet == 2) {
                            transformed->a[y][x] = tmp2->a[y][x];
                            transformed->b[y][x] = tmp2->b[y][x];
                        }
                    }

                    continue;
                }

#ifdef __SSE2__
                const float rhue = atan2Buffer[x];
                const float rchro = sqrtBuffer[x];
#else
                const float rhue = xatan2f(origblur->b[y][x], origblur->a[y][x]);
                const float rchro = sqrt(SQR(origblur->b[y][x]) + SQR(origblur->a[y][x])) / 327.68f;
#endif

                float rL = origblur->L[y][x] / 327.68f;
                float dE = sqrt(SQR(refa - origblur->a[y][x] / 327.68f) + SQR(refb - origblur->b[y][x] / 327.68f) + SQR(lumaref - rL));

                float cli = 1.f;
                float clc = 1.f;

                cli = (buflight[loy - begy][lox - begx]);
                clc = (bufchro[loy - begy][lox - begx]);

                float aplus = (1.f - cli) / delhu;
                float bplus = 1.f - aplus * hueplus;
                float amoins = (cli - 1.f) / delhu;
                float bmoins = 1.f - amoins * huemoins;

                float aplusch = (1.f - clc) / delhu;
                float bplusch = 1.f - aplusch * hueplus;
                float amoinsch = (clc - 1.f) / delhu;
                float bmoinsch = 1.f - amoinsch * huemoins;

                //prepare shape detection
                float kchchro = 1.f;
                float kch = 1.f;
                float fach = 1.f;
                float falu = 1.f;
                float deltachro = fabs(rchro - chromaref);
                float deltahue = fabs(rhue - hueref);
                float realstr = 1.f;
                float realstrch = 1.f;
                float khu = 0.f;

                if (deltahue > rtengine::RT_PI) {
                    deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280
                float deltaL = fabs(lumaref - rL);  //between 0 and 100

                /*
                                //kch to modulate action with chroma
                                if (deltachro < 160.f * SQR(lp.sensbn / 100.f)) {
                                    kch = 1.f;
                                } else {
                                    float ck = 160.f * SQR(lp.sensbn / 100.f);
                                    float ak = 1.f / (ck - 160.f);
                                    float bk = -160.f * ak;
                                    kch = ak * deltachro + bk;
                                }

                                if (lp.sensbn < 40.f) {
                                    kch = pow(kch, pa * lp.sensbn + pb);    //increase under 90
                                }

                                */
                //kch and kgcchro acts on luma and chroma
                if (deltachro < 160.f * SQR(lp.sensbn / 100.f)) {
                    kch = kchchro = 1.f;
                } else {
                    float ck = 160.f * SQR(lp.sensbn / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch  = kchchro = ak * deltachro + bk;
                }

                float kkch = 1.f;
                kkch = pa * lp.sensbn + pb;
                float kch0 = kch;

                if (lp.sensbn < 40.f) {
                    kch = kchchro = pow(kch0, kkch);    //increase under 40
                    float dEsens = moddE * lp.sensbn;
                    float kdE = 1.f;

                    if (dE > dEsens) {
                        kdE = 1.f;
                    } else {
                        kdE = SQR(SQR(dE / dEsens));
                    }

                    if (deltahue < 0.3f && settings->detectshape == true) {
                        kchchro = kch = pow(kch0, (1.f + kdE * moddE * (3.f - 10.f * deltahue)) * kkch);
                    }

                    /*
                    if (deltahue < 0.3f && dE > (2.5f * lp.sensbn)  && settings->detectshape == true) {
                        float epsi = 0.000001f;
                        kchchro = kch = pow(kch0, (1.f + 5.f * (3.f - 10.f * deltahue)) * kkch);//increase more and more if deltahue very small and chrom diffrent
                    }
                    */
                }

                float dEsensall = lp.sensbn;
                float kD = 1.f;

                if (settings->detectshape == true) {
                    if (dE < dEsensall) {
                        kD = 1.f;
                    } else {
                        kD = pow(dEsensall / dE, powdE);
                    }
                }

                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if (lp.qualmet >= 1 && lp.sensbn < 50.f) { //to try...
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
                            realstrch = clc;
                            khu = 1.f;
                        }
                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
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
                            realstrch = clc;
                            khu = 1.f;
                        }
                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstrch = amoinsch * rhue + bmoinsch;
                            realstr = amoins * rhue + bmoins;
                            khu = amo * rhue + bmo;
                        } else {
                            realstr = cli;
                            realstrch = clc;
                            khu = 1.f;
                        }
                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;
                        } else {
                            realstrch = clc;
                            realstr = cli;
                            khu = 1.f;
                        }
                    }

                    realstr *= kD;
                    realstrch *= kD;

                    if (settings->detectshape == false) {

                        if (lp.sensbn <= 20.f) { //to try...

                            if (deltaE <  2.8f * lp.sensbn) {
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

                            if (deltaL <  lp.sensbn) {
                                falu = 1.f;
                            } else {
                                falu = alum * deltaL + blum;
                            }

                        }
                    }
                }


                switch (zone) {

                    case 1: { // inside transition zone
                        float difL, difa, difb;
                        float factorx = localFactor;

                        if (call <= 3) {
                            difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                            difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                        } else {
                            difL = tmp1->L[y][x] - original->L[y][x];
                            difa = tmp1->a[y][x] - original->a[y][x];
                            difb = tmp1->b[y][x] - original->b[y][x];


                        }

                        difL *= factorx * (100.f + realstr * kchchro) / 100.f;
                        difL *= kch * fach;

                        if (lp.blurmet == 0) {
                            transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                        }

                        if (lp.blurmet == 2) {
                            transformed->L[y][x] = CLIP(tmp2->L[y][x] - difL);
                        }

                        if (!lp.actsp) {
                            difa *= factorx * (100.f + realstrch * falu * kchchro) / 100.f;
                            difb *= factorx * (100.f + realstrch * falu * kchchro) / 100.f;
                            difa *= kch * fach;
                            difb *= kch * fach;

                            if (lp.blurmet == 0) {
                                transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                            }

                            if (lp.blurmet == 2) {
                                transformed->a[y][x] = CLIPC(tmp2->a[y][x] - difa);
                                transformed->b[y][x] = CLIPC(tmp2->b[y][x] - difb);
                            }
                        }

                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        float difL, difa, difb;

                        if (call <= 3) {
                            difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                            difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                            difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                        } else {
                            difL = tmp1->L[y][x] - original->L[y][x];
                            difa = tmp1->a[y][x] - original->a[y][x];
                            difb = tmp1->b[y][x] - original->b[y][x];

                        }

                        difL *= (100.f + realstr * kchchro) / 100.f;
                        difL *= kch * fach;

                        if (lp.blurmet == 0) {
                            transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                        }

                        if (lp.blurmet == 2) {
                            transformed->L[y][x] = CLIP(tmp2->L[y][x] - difL);
                        }

                        if (!lp.actsp) {
                            difa *= (100.f + realstrch * falu * kchchro) / 100.f;
                            difb *= (100.f + realstrch * falu * kchchro) / 100.f;
                            difa *= kch * fach;
                            difb *= kch * fach;

                            if (lp.blurmet == 0) {
                                transformed->a[y][x] = CLIPC(original->a[y][x] + difa);                                ;
                                transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                            }

                            if (lp.blurmet == 2) {
                                transformed->a[y][x] = CLIPC(tmp2->a[y][x] - difa);
                                transformed->b[y][x] = CLIPC(tmp2->b[y][x] - difb);
                            }

                        }
                    }
                }
            }
        }
    }
    delete origblur;
}

void ImProcFunctions::InverseReti_Local(const struct local_params & lp, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy, int chro)
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

            if (lp.shapmet == 0) {
                calcTransition(lox, loy, ach, lp, zone, localFactor);
            } else if (lp.shapmet == 1) {
                calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
            }


            switch (zone) {
                case 0: { // outside selection and outside transition zone => full effect, no transition
                    if (chro == 0) {
                        transformed->L[y][x] = CLIP(tmp1->L[y][x]);
                    }

                    if (chro == 1) {

                        transformed->a[y][x] = CLIPC(tmp1->a[y][x]);
                        transformed->b[y][x] = CLIPC(tmp1->b[y][x]);
                    }

                    break;
                }

                case 1: { // inside transition zone
                    float factorx = 1.f - localFactor;

                    if (chro == 0) {
                        float difL = tmp1->L[y][x] - original->L[y][x];
                        difL *= factorx;
                        transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                    }

                    if (chro == 1) {
                        float difa = tmp1->a[y][x] - original->a[y][x];
                        float difb = tmp1->b[y][x] - original->b[y][x];

                        difa *= factorx;
                        difb *= factorx;

                        transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                        transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
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




void ImProcFunctions::InverseBlurNoise_Local(const struct local_params & lp, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy)
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

            if (lp.shapmet == 0) {
                calcTransition(lox, loy, ach, lp, zone, localFactor);
            } else if (lp.shapmet == 1) {
                calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
            }


            switch (zone) {
                case 0: { // outside selection and outside transition zone => full effect, no transition
                    transformed->L[y][x] = CLIP(tmp1->L[y][x]);

                    if (!lp.actsp) {
                        transformed->a[y][x] = CLIPC(tmp1->a[y][x]);
                        transformed->b[y][x] = CLIPC(tmp1->b[y][x]);
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

                    transformed->L[y][x] = CLIP(original->L[y][x] + difL);

                    if (!lp.actsp) {

                        transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                        transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
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

static void calclight(float lum, float  koef, float & lumnew, LUTf & lightCurveloc)
//replace L-curve that does not work in local or bad
{
    if (koef >= 0.f) {
        lumnew = lightCurveloc[lum];
    }

    if (koef < 0.f) {
        lumnew = lightCurveloc[lum];

        if (koef == -100.f) {
            lumnew = 0.f;
        }

    }

    lumnew = CLIPLOC(lumnew);

}


void ImProcFunctions::InverseSharp_Local(float **loctemp, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
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
    int GW = transformed->W;
    int GH = transformed->H;

    LabImage *origblur = nullptr;

    origblur = new LabImage(GW, GH);

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
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < transformed->H; y++) {
#ifdef __SSE2__
            int i = 0;

            for (; i < transformed->W - 3; i += 4) {
                vfloat av = LVFU(origblur->a[y][i]);
                vfloat bv = LVFU(origblur->b[y][i]);
                STVF(atan2Buffer[i], xatan2f(bv, av));
                STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f(origblur->b[y][i], origblur->a[y][i]);
                sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
            }

#endif

            int loy = cy + y;

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;
#ifdef __SSE2__
                float rhue = atan2Buffer[x];
                float rchro = sqrtBuffer[x];
#else
                float rhue = xatan2f(origblur->b[y][x], origblur->a[y][x]);
                float rchro = sqrt(SQR(origblur->b[y][x]) + SQR(origblur->a[y][x])) / 327.68f;
#endif
                int zone;
                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else if (lp.shapmet == 1) {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                //prepare shape detection
                float khu = 0.f;
                float kch = 1.f;
                float fach = 1.f;
                float deltachro = fabs(rchro - chromaref);
                float deltahue = fabs(rhue - hueref);

                if (deltahue > rtengine::RT_PI) {
                    deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                //kch to modulate action with chroma
                if (deltachro < 160.f * SQR(lp.senssha / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR(lp.senssha / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;
                }

                if (lp.senssha < 40.f) {
                    kch = pow(kch, pa * lp.senssha + pb);    //increase under 40
                }


                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if (lp.senssha < 20.f) { //to try...
                    //hue detection
                    if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
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
                        transformed->L[y][x] = CLIP(original->L[y][x] + difL * kch * fach);

                        break;
                    }

                    case 1: { // inside transition zone
                        float difL = loctemp[y][x] - original->L[y][x];

                        float factorx = 1.f - localFactor;
                        difL *= factorx;

                        transformed->L[y][x] = CLIP(original->L[y][x] + difL * kch * fach);
                        break;
                    }

                    case 2: { // inside selection => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                    }
                }
            }
        }
    }
    delete origblur;
}


void ImProcFunctions::Sharp_Local(int call, float **loctemp,  int senstype, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{
    BENCHFUN
    const float ach = (float)lp.trans / 100.f;
    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;
    float varsens = 0.f;

    if (senstype == 0) {
        varsens = lp.senssha;

    } else if (senstype == 0) {
        varsens = lp.senslc;
    }

    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * varsens - 280.f);
    const float bhu = 1.f - ahu * 2.8f * varsens;

    const bool detectHue = varsens < 20.f && lp.qualmet >= 1;

    int GW = transformed->W;
    int GH = transformed->H;

    LabImage *origblur = nullptr;

    origblur = new LabImage(GW, GH);

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
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
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
                    vfloat av = LVFU(origblur->a[y][i]);
                    vfloat bv = LVFU(origblur->b[y][i]);
                    STVF(atan2Buffer[i], xatan2f(bv, av));
                    STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    atan2Buffer[i] = xatan2f(origblur->b[y][i], origblur->a[y][i]);
                    sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
                }
            } else {
                for (; i < transformed->W - 3; i += 4) {
                    vfloat av = LVFU(origblur->a[y][i]);
                    vfloat bv = LVFU(origblur->b[y][i]);
                    STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
                }
            }

#endif

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;
                int zone = 0;
                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else if (lp.shapmet == 1) {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    continue;
                }

#ifdef __SSE2__
                float rchro = sqrtBuffer[x];
#else
                float rchro = sqrt(SQR(origblur->b[y][x]) + SQR(origblur->a[y][x])) / 327.68f;
#endif
                //prepare shape detection
                float kch = 1.f;
                float fach = 1.f;
                float deltachro = fabs(rchro - chromaref);

                //kch to modulate action with chroma
                if (deltachro < 160.f * SQR(varsens / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR(varsens / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;

                    if (varsens < 40.f) {
                        kch = pow_F(kch, pa * varsens + pb);    //increase under 40
                    }
                }

                if (varsens >= 99.f) {
                    kch = 1.f;
                }

                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if (detectHue) { //to try...
#ifdef __SSE2__
                    float rhue = atan2Buffer[x];
#else
                    float rhue = xatan2f(origblur->b[y][x], origblur->a[y][x]);
#endif
                    float khu = 0.f;
                    float deltahue = fabs(rhue - hueref);

                    if (deltahue > rtengine::RT_PI) {
                        deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                    }

                    //hue detection
                    if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }


                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            khu  = apl * rhue + bpl;
                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            khu = amo * rhue + bmo;
                        } else {
                            khu = 1.f;
                        }

                    }

                    float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                    if (deltaE <  2.8f * varsens) {
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
                        transformed->L[y][x] = CLIP(original->L[y][x] + difL * kch * fach);

                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        float difL;

                        if (call == 2) {
                            difL = loctemp[loy - begy][lox - begx] - original->L[y][x];
                        } else  {
                            difL = loctemp[y][x] - original->L[y][x];
                        }

                        transformed->L[y][x] = CLIP(original->L[y][x] + difL * kch * fach);
                    }
                }
            }
        }
    }
    delete origblur;
}



void ImProcFunctions::Exclude_Local(int sen, float **deltaso, float **buflight, float **bufchro, const float hueref, const float chromaref, const float lumaref, const struct local_params & lp, LabImage * original, LabImage * transformed, LabImage * rsv, LabImage * reserv, int cx, int cy, int sk)
{

    BENCHFUN {
        const float ach = (float)lp.trans / 100.f;
        float varsens =  lp.sensexclu;

        if (sen == 1)
        {
            varsens =  lp.sensexclu;
        }

        int GW = transformed->W;
        int GH = transformed->H;

        float refa = chromaref * cos(hueref);
        float refb = chromaref * sin(hueref);


        LabImage *origblur = nullptr;

        origblur = new LabImage(GW, GH);

        float radius = 3.f / sk;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(reserv->L, origblur->L, GW, GH, radius);
            gaussianBlur(reserv->a, origblur->a, GW, GH, radius);
            gaussianBlur(reserv->b, origblur->b, GW, GH, radius);

        }

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
#ifdef __SSE2__
            float atan2Buffer[transformed->W] ALIGNED16;
            float sqrtBuffer[transformed->W] ALIGNED16;
            vfloat c327d68v = F2V(327.68f);
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
                    vfloat av = LVFU(origblur->a[y][i]);
                    vfloat bv = LVFU(origblur->b[y][i]);
                    STVF(atan2Buffer[i], xatan2f(bv, av));
                    STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    atan2Buffer[i] = xatan2f(origblur->b[y][i], origblur->a[y][i]);
                    sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
                }

#endif

                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int begx = int (lp.xc - lp.lxL);
                    int begy = int (lp.yc - lp.lyT);

                    int zone = 0;
                    float localFactor = 1.f;

                    if (lp.shapmet == 0) {
                        calcTransition(lox, loy, ach, lp, zone, localFactor);
                    } else if (lp.shapmet == 1) {
                        calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                    }


                    if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        continue;
                    }

#ifdef __SSE2__
//                   float rhue = atan2Buffer[x];
#else
//                    float rhue = xatan2f(origblur->b[y][x], origblur->a[y][x]);
#endif
                    float rL = origblur->L[y][x] / 327.68f;
                    float dE = sqrt(SQR(refa - origblur->a[y][x] / 327.68f) + SQR(refb - origblur->b[y][x] / 327.68f) + SQR(lumaref - rL));

                    float cli = 1.f;
                    float clc = 1.f;

                    cli = (buflight[loy - begy][lox - begx]);
                    clc = (bufchro[loy - begy][lox - begx]);

                    float reducdE = 0.f;
                    float mindE = 2.f + 0.05f * varsens;//between 2 and 7
                    float maxdE = 5.f + 1.5f * varsens; // between 5 and 150, we can chnage this values

                    float ar = 1.f / (mindE - maxdE);

                    float br = - ar * maxdE;

                    if (dE > maxdE) {
                        reducdE = 0.f;
                    }

                    if (dE > mindE && dE <= maxdE) {
                        reducdE = ar * dE + br;
                    }

                    if (dE <= mindE) {
                        reducdE = 1.f;
                    }


                    float realstrdE = reducdE * cli;
                    float realstrchdE = reducdE * clc;


                    if (rL > 0.1f) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                        switch (zone) {
                            case 0: { // outside selection and outside transition zone => no effect, keep original values
                                transformed->L[y][x] = original->L[y][x];
                                transformed->a[y][x] = original->a[y][x];
                                transformed->b[y][x] = original->b[y][x];

                                break;
                            }

                            case 1: { // inside transition zone
                                float factorx = localFactor;

                                float difL;
                                difL = rsv->L[loy - begy][lox - begx] - original->L[y][x];
                                //  difL *= factorx * (100.f + realstr * falL) / 100.f;
                                // difL *= kch * fach;
                                difL *= factorx * (100.f + realstrdE) / 100.f;

                                //   difL *= kch * fach;

                                if (deltaso[loy - begy][lox - begx] == 0.f) {
                                    transformed->L[y][x]  = original->L[y][x];
                                } else {
                                    transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                                }

                                float difa, difb;

                                difa = rsv->a[loy - begy][lox - begx] - original->a[y][x];
                                difb = rsv->b[loy - begy][lox - begx] - original->b[y][x];
                                //   difa *= factorx * (100.f + realstrch * falu * falL) / 100.f;
                                //  difb *= factorx * (100.f + realstrch * falu * falL) / 100.f;
                                difa *= factorx * (100.f +  realstrchdE) / 100.f;
                                difb *= factorx * (100.f +  realstrchdE) / 100.f;

                                //    difa *= kch * fach;
                                //    difb *= kch * fach;

                                if (deltaso[loy - begy][lox - begx] == 0.f) {
                                    transformed->a[y][x] = original->a[y][x]; //rsv->a[loy - begy][lox - begx];
                                    transformed->b[y][x] = original->b[y][x]; //rsv->b[loy - begy][lox - begx];
                                } else {
                                    transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                    transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                                }

                                break;

                            }

                            case 2: { // inside selection => full effect, no transition
                                float difL;

                                difL = rsv->L[loy - begy][lox - begx] - original->L[y][x];
                                difL *= (100.f + realstrdE) / 100.f;

                                //    difL *= kch * fach;

                                if (deltaso[loy - begy][lox - begx] == 0.f) {
                                    //  printf ("0");
                                    transformed->L[y][x]  = original->L[y][x]; //rsv->L[loy - begy][lox - begx];
                                } else {
                                    transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                                }

                                //  transformed->L[y][x] = original->L[y][x] + difL;
                                float difa, difb;

                                difa = rsv->a[loy - begy][lox - begx] - original->a[y][x];
                                difb = rsv->b[loy - begy][lox - begx] - original->b[y][x];
                                difa *= (100.f + realstrchdE) / 100.f;
                                difb *= (100.f + realstrchdE) / 100.f;

                                //   difa *= kch * fach;
                                //   difb *= kch * fach;

                                if (deltaso[loy - begy][lox - begx] == 0.f) {
                                    // printf ("0");
                                    transformed->a[y][x] = original->a[y][x]; //rsv->a[loy - begy][lox - begx];
                                    transformed->b[y][x] = original->b[y][x]; //rsv->b[loy - begy][lox - begx];
                                } else {
                                    // printf ("1");
                                    transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                    transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                                }

                            }
                        }

                    }

                }
            }

        }
        delete origblur;
    }
}


void ImProcFunctions::transit_shapedetect(int senstype, LabImage * bufexporig, LabImage * originalmask, float **buflight, float **bufchro, float **buf_a_cat, float ** buf_b_cat, float ** bufhh, bool HHutili, const float hueref, const float chromaref,  const float lumaref, float sobelref, float ** blend2, const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{

    BENCHFUN {
        const float ach = (float)lp.trans / 100.f;
        float varsens =  lp.sensex;

        if (senstype == 0) //Color and Light
        {
            varsens =  lp.sens;
        }

        if (senstype == 1) //exposure
        {
            varsens =  lp.sensex;
        }

        if (senstype == 2) //vibrance
        {
            varsens =  lp.sensv;
        }

        if (senstype == 3) //soft light
        {
            varsens =  lp.senssf;
        }

        if (senstype == 4 || senstype == 5) //retinex
        {
            varsens =  lp.sensh;
        }

        if (senstype == 6 || senstype == 7) //cbdl
        {
            varsens =  lp.senscb;
        }

        //  printf("varsen=%f \n", varsens);
        //sobel
        sobelref /= 100.;

        if (sobelref > 60.)
        {
            sobelref = 60.;
        }

        sobelref = log(1.f + sobelref);

        int GW = transformed->W;
        int GH = transformed->H;

        float refa = chromaref * cos(hueref);
        float refb = chromaref * sin(hueref);

        bool expshow = ((lp.showmaskexpmet == 1 || lp.showmaskexpmet == 2)  &&  senstype == 1);
        bool colshow = ((lp.showmaskcolmet == 1 || lp.showmaskcolmet == 2)  &&  senstype == 0);


        LabImage *origblur = nullptr;

        origblur = new LabImage(GW, GH);
        LabImage *origblurmask = nullptr;

        float radius = 3.f;

        if (senstype == 1)
        {
            radius = (2.f + 0.2f * lp.blurexp) / sk;
        }

        if (senstype == 0)
        {
            radius = (2.f + 0.2f * lp.blurcol) / sk;
        }

        bool usemask = (lp.showmaskexpmet == 2 || lp.enaExpMask) && senstype == 1;
        bool usemaskcol = (lp.showmaskcolmet == 2 || lp.enaColorMask) && senstype == 0;

        if (usemask  || usemaskcol)
        {
            origblurmask = new LabImage(GW, GH);

#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                gaussianBlur(originalmask->L, origblurmask->L, GW, GH, radius);
                gaussianBlur(originalmask->a, origblurmask->a, GW, GH, radius);
                gaussianBlur(originalmask->b, origblurmask->b, GW, GH, radius);
            }
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
#ifdef __SSE2__
            float atan2Buffer[transformed->W] ALIGNED16;
            float sqrtBuffer[transformed->W] ALIGNED16;
            vfloat c327d68v = F2V(327.68f);
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
                    vfloat av = LVFU(origblur->a[y][i]);
                    vfloat bv = LVFU(origblur->b[y][i]);
                    STVF(atan2Buffer[i], xatan2f(bv, av));
                    STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    atan2Buffer[i] = xatan2f(origblur->b[y][i], origblur->a[y][i]);
                    sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
                }

#endif

                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int begx = int (lp.xc - lp.lxL);
                    int begy = int (lp.yc - lp.lyT);

                    int zone = 0;
                    float localFactor = 1.f;

                    if (lp.shapmet == 0) {
                        calcTransition(lox, loy, ach, lp, zone, localFactor);
                    } else if (lp.shapmet == 1) {
                        calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                    }


                    if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                        //        transformed->L[y][x] = original->L[y][x];
                        continue;
                    }

#ifdef __SSE2__
                    float rhue = atan2Buffer[x];
#else
                    float rhue = xatan2f(origblur->b[y][x], origblur->a[y][x]);
#endif

                    float rL = origblur->L[y][x] / 327.68f;
                    float csob = 0.f;
                    float rs = 0.f;

                    if (senstype == 1  || senstype == 0) {
                        csob = (blend2[loy - begy][lox - begx]) / 100.f ;

                        if (csob > 60.f) {
                            csob = 60.f;
                        }

                        csob = log(1.f + csob + 0.001f);
                        rs = sobelref / csob;
                    }

                    float dE = 0.f;
                    float rsob = 0.f;

                    if (lp.struexp > 0.f && rs > 0.f && senstype == 1) {
                        rsob =  1.1f * lp.struexp * rs;
                    }

                    if (lp.struco > 0.f && rs > 0.f && senstype == 0) {
                        rsob =  1.1f * lp.struco * rs;
                    }

                    if (usemask  || usemaskcol) {
                        dE = rsob + sqrt(SQR(refa - origblurmask->a[y][x] / 327.68f) + SQR(refb - origblurmask->b[y][x] / 327.68f) + SQR(lumaref - origblurmask->L[y][x] / 327.68f));
                    } else {
                        dE = rsob + sqrt(SQR(refa - origblur->a[y][x] / 327.68f) + SQR(refb - origblur->b[y][x] / 327.68f) + SQR(lumaref - rL));
                    }

                    float cli = 0.f;
                    float clc = 0.f;
                    float cla = 0.f;
                    float clb = 0.f;
                    float hhro = 0.f;

                    if (HHutili) {
                        hhro = bufhh[loy - begy][lox - begx];
                    }

                    cli = (buflight[loy - begy][lox - begx]);
                    clc = (bufchro[loy - begy][lox - begx]);

                    if (senstype == 1  || senstype == 0) {
                        cla = buf_a_cat[loy - begy][lox - begx];
                        clb = buf_b_cat[loy - begy][lox - begx];
                        //   printf("a=%f b=%f", cla, clb);
                    }

                    float reducdE = 0.f;
                    float mindE = 2.f + 0.05f * varsens;//between 2 and 7
                    float maxdE = 5.f + 1.5f * varsens; // between 5 and 150, we can chnage this values

                    float ar = 1.f / (mindE - maxdE);

                    float br = - ar * maxdE;

                    if (dE > maxdE) {
                        reducdE = 0.f;
                    }

                    if (dE > mindE && dE <= maxdE) {
                        reducdE = ar * dE + br;
                    }

                    if (dE <= mindE) {
                        reducdE = 1.f;
                    }


                    float realstrdE = reducdE * cli;
                    float realstradE = reducdE * cla;
                    float realstrbdE = reducdE * clb;
                    float realstrchdE = reducdE * clc;
                    float realhhdE = reducdE * hhro;

                    float addh = 0.f;
//                    float chromhr = sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x]));
//                    float2 sincosval;
//                    sincosval.y = 1.f;
//                    sincosval.x = 0.0f;

                    if (rL > 0.1f) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9



                        switch (zone) {
                            case 0: { // outside selection and outside transition zone => no effect, keep original values
                                transformed->L[y][x] = original->L[y][x];
                                transformed->a[y][x] = original->a[y][x];
                                transformed->b[y][x] = original->b[y][x];

                                break;
                            }

                            case 1: { // inside transition zone
                                float factorx = localFactor;
                                float diflc = 0.f;
                                float newhr = 0.f;

                                if (senstype == 4  || senstype == 6 || senstype == 2 || senstype == 3) {//all except color and light (TODO) and exposure
                                    float lightc = bufexporig->L[loy - begy][lox - begx];

                                    float fli = ((100.f + realstrdE) / 100.f);
                                    float diflc = lightc * fli - original->L[y][x];
                                    diflc *= factorx;
                                    transformed->L[y][x] = CLIP(original->L[y][x] + diflc);
                                } else if (senstype == 1 || senstype == 0) {//exposure the TODO color and light
                                    transformed->L[y][x] = CLIP(original->L[y][x] + 328.f * factorx * realstrdE);//kch fach
                                    diflc = 328.f * factorx * realstrdE;
                                }

                                if (HHutili) {
                                    addh = 0.01f * realhhdE * factorx;
                                    newhr = rhue + addh;

                                    if (newhr > rtengine::RT_PI) {
                                        newhr -= 2 * rtengine::RT_PI;
                                    } else if (newhr < -rtengine::RT_PI) {
                                        newhr  += 2 * rtengine::RT_PI;
                                    }
                                }

                                if (senstype == 7) {
                                    float difab = bufexporig->L[loy - begy][lox - begx] - sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x]));
                                    float difa = difab * cos(rhue);
                                    float difb = difab * sin(rhue);
                                    difa *= factorx * (100.f + realstrchdE) / 100.f;
                                    difb *= factorx * (100.f + realstrchdE) / 100.f;
                                    transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                    transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                                } /* else if (senstype == 0 && HHutili) {

                                float fac = (100.f + factorx * realstrchdE) / 100.f;

                                sincosval = xsincosf(newhr);
                                transformed->a[y][x] = CLIPC(chromhr * sincosval.y * fac) ;
                                transformed->b[y][x] = CLIPC(chromhr * sincosval.x * fac);

                            } */else {

                                    float flia = 1.f;
                                    float flib = 1.f;
                                    float chra = bufexporig->a[loy - begy][lox - begx];
                                    float chrb = bufexporig->b[loy - begy][lox - begx];

                                    if (senstype == 4  || senstype == 6 || senstype == 2 || senstype == 3 || senstype == 0) {
                                        flia = flib = ((100.f + realstrchdE) / 100.f);
                                    } else if (senstype == 1) {
                                        flia = (100.f + realstradE) / 100.f;
                                        flib = (100.f + realstrbdE) / 100.f;
                                    }

                                    float difa = chra * flia - original->a[y][x];
                                    float difb = chrb * flib - original->b[y][x];
                                    difa *= factorx;
                                    difb *= factorx;

                                    transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                    transformed->b[y][x] = CLIPC(original->b[y][x] + difb);

                                    if (expshow  || colshow) {
                                        transformed->L[y][x] = CLIP(12000.f + diflc);
                                        transformed->a[y][x] = CLIPC(difa);
                                        transformed->b[y][x] = CLIPC(difb);
                                    }
                                }

                                break;

                            }

                            case 2: { // inside selection => full effect, no transition
                                float diflc = 0.f;
                                float newhr = 0.f;

                                if (senstype == 4  || senstype == 6  || senstype == 2 || senstype == 3) {//retinex & cbdl
                                    float lightc = bufexporig->L[loy - begy][lox - begx];
                                    float fli = ((100.f + realstrdE) / 100.f);
                                    float diflc = lightc * fli - original->L[y][x];
                                    transformed->L[y][x] = CLIP(original->L[y][x] + diflc);
                                } else if (senstype == 1 || senstype == 0) {
                                    transformed->L[y][x] = CLIP(original->L[y][x] + 328.f * realstrdE);//kch fach
                                    diflc = 328.f * realstrdE;
                                }

                                if (HHutili) {
                                    addh = 0.01f * realhhdE;
                                    newhr = rhue + addh;

                                    if (newhr > rtengine::RT_PI) {
                                        newhr -= 2 * rtengine::RT_PI;
                                    } else if (newhr < -rtengine::RT_PI) {
                                        newhr  += 2 * rtengine::RT_PI;
                                    }
                                }

                                if (senstype == 7) {//cbdl chroma
                                    float difab = bufexporig->L[loy - begy][lox - begx] - sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x]));
                                    float difa = difab * cos(rhue);
                                    float difb = difab * sin(rhue);
                                    difa *= (100.f + realstrchdE) / 100.f;
                                    difb *= (100.f + realstrchdE) / 100.f;
                                    transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                    transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                                }/* else if (senstype == 0 && HHutili) {

                                float fac = (100.f + realstrchdE) / 100.f;

                                sincosval = xsincosf(newhr);
                                transformed->a[y][x] = CLIPC(chromhr * sincosval.y * fac) ;
                                transformed->b[y][x] = CLIPC(chromhr * sincosval.x * fac);
                          */else {
                                    float flia = 1.f;
                                    float flib = 1.f;
                                    float chra = bufexporig->a[loy - begy][lox - begx];
                                    float chrb = bufexporig->b[loy - begy][lox - begx];

                                    if (senstype == 4  || senstype == 6 || senstype == 2 || senstype == 3 || senstype == 0) {
                                        flia = flib = (100.f + realstrchdE) / 100.f;
                                    } else if (senstype == 1) {
                                        flia = (100.f + realstradE) / 100.f;
                                        flib = (100.f + realstrbdE) / 100.f;
                                    }

                                    float difa = chra * flia - original->a[y][x];
                                    float difb = chrb * flib - original->b[y][x];

                                    transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                    transformed->b[y][x] = CLIPC(original->b[y][x] + difb);

                                    if (expshow  || colshow) {
                                        transformed->L[y][x] = CLIP(12000.f + diflc);
                                        transformed->a[y][x] = CLIPC(difa);
                                        transformed->b[y][x] = CLIPC(difb);
                                    }
                                }
                            }
                        }


                    }
                }

            }

            bool execmedian99 = false;

            if (execmedian99)
                //I tested here median to see if action on artifacts...when color differences due to WB or black... or mixed color or ??
                //small action with 9x9 3 times
                //warm cool is hugely better
            {
                float** tmL;
                int wid = transformed->W;
                int hei = transformed->H;
                tmL = new float*[hei];

                for (int i = 0; i < hei; ++i) {
                    tmL[i] = new float[wid];
                }

                Median medianTypeL = Median::TYPE_9X9;
                Median medianTypeAB = Median::TYPE_9X9;


                Median_Denoise(transformed->L, transformed->L, transformed->W, transformed->H, medianTypeL, 3, multiThread, tmL);
                Median_Denoise(transformed->a, transformed->a, transformed->W, transformed->H, medianTypeAB, 3, multiThread, tmL);
                Median_Denoise(transformed->b, transformed->b, transformed->W, transformed->H, medianTypeAB, 3, multiThread, tmL);

                for (int i = 0; i < hei; ++i) {
                    delete[] tmL[i];
                }

                delete[] tmL;
            }

        }
        delete origblur;

        if ((lp.showmaskcolmet == 2 || lp.enaColorMask) && senstype == 1)
        {
            delete origblurmask;
        }

    }
}

void ImProcFunctions::InverseColorLight_Local(const struct local_params & lp, LUTf & lightCurveloc, LabImage * original, LabImage * transformed, int cx, int cy, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, int sk)
{
    // BENCHFUN
    float ach = (float)lp.trans / 100.f;
    const float facc = (100.f + lp.chro) / 100.f; //chroma factor transition

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
    constexpr float delhu = 0.05f; //between 0.05 and 0.2 ==> minima for scope

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
    ImProcFunctions::secondeg_begin(reducac, vi, aa, bb); //parabolic
    ImProcFunctions::secondeg_end(reducac, vinf, aaa, bbb, ccc); //parabolic
    float vinf2 = (lumaref + modlum) / 100.f;
    float vi2 = (lumaref - modlum) / 100.f;
    float aaaa, bbbb, cccc, aO, bO;
    ImProcFunctions::secondeg_end(reducac2, vinf2, aaaa, bbbb, cccc); //parabolic
    ImProcFunctions::secondeg_begin(reducac2, vi2, aO, bO); //parabolic

    int GW = transformed->W;
    int GH = transformed->H;

    LabImage *origblur = nullptr;

    origblur = new LabImage(GW, GH);

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
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < transformed->H; y++) {
            const int loy = cy + y;

#ifdef __SSE2__
            int i = 0;

            for (; i < transformed->W - 3; i += 4) {
                vfloat av = LVFU(origblur->a[y][i]);
                vfloat bv = LVFU(origblur->b[y][i]);
                STVF(atan2Buffer[i], xatan2f(bv, av));
                STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
            }

            for (; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f(origblur->b[y][i], origblur->a[y][i]);
                sqrtBuffer[i] = sqrt(SQR(origblur->b[y][i]) + SQR(origblur->a[y][i])) / 327.68f;
            }

#endif


            for (int x = 0; x < transformed->W; x++) {
                const int lox = cx + x;
                int zone = 0;

                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else if (lp.shapmet == 1) {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);//rect not good
                }



#ifdef __SSE2__
                float rhue = atan2Buffer[x];
                float rchro = sqrtBuffer[x];
#else

                float rhue = xatan2f(origblur->b[y][x], origblur->a[y][x]);

                float rchro = sqrt(SQR(origblur->b[y][x]) + SQR(origblur->a[y][x])) / 327.68f;
#endif

                float rL = origblur->L[y][x] / 327.68f;
                float rLL = origblur->L[y][x] / 327.68f;

                if (fabs(origblur->b[y][x]) < 0.01f) {
                    origblur->b[y][x] = 0.01f;
                }


                float realchro = 1.f;
                //evaluate delta Hue and delta Chro
                float deltachro = fabs(rchro - chromaref);

                float deltahue = fabs(rhue - hueref);

                if (deltahue > rtengine::RT_PI) {
                    deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                }

                //pseudo deltaE
                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280
                float deltaL = fabs(lumaref - rL);  //between 0 and 100

                float kch = 1.f;
                float khu = 0.f;
                float fach = 1.f;
                float falu = 1.f;

                //kch acts on luma
                if (deltachro < 160.f * SQR(lp.sens / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR(lp.sens / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;
                }

                if (lp.sens < 40.f) {
                    kch = pow(kch, pa * lp.sens + pb);    //increase under 40
                }

                bool kzon = false;

                if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                    if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;
                        khu  = apl * rhue + bpl;

                    } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;
                        khu = amo * rhue + bmo;

                    } else {
                        realchro = lp.chro;
                        khu = 1.f;
                    }

                    kzon = true;
                } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
                    if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;
                        khu  = apl * rhue + bpl;

                    } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;
                        khu = amo * rhue + bmo;

                    } else {
                        realchro = lp.chro;
                        khu = 1.f;
                    }

                    kzon = true;
                }

                if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                    if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;
                        khu  = apl * rhue + bpl;

                    } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;
                        khu = amo * rhue + bmo;

                    } else {
                        realchro = lp.chro;
                        khu = 1.f;
                    }

                    kzon = true;
                } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus)) {
                    if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;
                        khu  = apl * rhue + bpl;

                    } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;
                        khu = amo * rhue + bmo;

                    } else {
                        realchro = lp.chro;
                        khu = 1.f;
                    }

                    kzon = true;
                }

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
                                if ((rhue > -2.79f && rhue < -1.11f) && (rchro < 35.f * multchro)) {
                                    realchro *= 0.9f;
                                } else {
                                    realchro = 1.f;

                                }
                            }
                        } else {
                            realchro = lp.chro;

                        }

                        if (lp.sens < 50.f && lp.chro > 0.f) {
                            if (hueref > -0.1f && hueref < 1.6f) { // detect skin
                                if (chromaref > 0.f && chromaref < 55.f * multchroskin) { // detect skin
                                    if ((rhue > -0.09f && rhue < 1.59f) && (rchro < 55.f * multchroskin)) {
                                        realchro *= 0.9f;

                                    } else {
                                        realchro = 1.f;

                                    }
                                }
                            } else {
                                realchro = lp.chro;

                            }
                        }
                    }

                }

                float kLinf = rLL / (100.f);
                float kLsup = kLinf;

                float kdiff = 1.f;

                if (kzon) { ///rhue < hueplus && rhue > huemoins

                    if ((rLL > (lumaref - modlum) && rLL < (lumaref + modlum))) {
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

                    if ((rLL > (lumaref - modlum) && rLL < (lumaref + modlum))) {
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
                    falL *= pow(rchro / kcr, lp.iterat / 10.f);
                }


                float th_r = 0.01f;

                if (rL > th_r) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9

                    switch (zone) {
                        case 2: { // outside selection and outside transition zone => no effect, keep original values
                            transformed->L[y][x] = original->L[y][x];
                            transformed->a[y][x] = original->a[y][x];
                            transformed->b[y][x] = original->b[y][x];
                            break;
                        }

                        case 1: { // inside transition zone
                            float lumnew = original->L[y][x];

                            if (lp.sens < 75.f) {
                                float lightcont;

                                if ((lp.ligh != 0.f || lp.cont != 0)) {
                                    calclight(lumnew, lp.ligh, lumnew, lightCurveloc);  //replace L-curve
                                    lightcont = lumnew;

                                } else {
                                    lightcont = lumnew;
                                }

                                float factorx = 1.f - localFactor;
                                float fli = 1.f;
                                float flisl = 1.f;
                                float flicur = 1.f;
                                float fac = flicur * (100.f + factorx * realchro * falu * falL) / 100.f;  //chroma factor transition
                                float diflc = lightcont * fli * flisl - original->L[y][x];
                                kdiff *= fach * kch;
                                diflc *= kdiff ;

                                diflc *= factorx; //transition lightness
                                transformed->L[y][x] = CLIP(1.f * (original->L[y][x] + diflc));

                                transformed->a[y][x] = CLIPC(original->a[y][x] * fac) ;
                                transformed->b[y][x] = CLIPC(original->b[y][x] * fac);
                            } else {
                                float factorx = 1.f - localFactor;
                                float fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition
                                float lumnew = original->L[y][x];

                                if ((lp.ligh != 0.f || lp.cont != 0)) {
                                    calclight(original->L[y][x], lp.ligh, lumnew, lightCurveloc);
                                }

                                float lightcont = lumnew ; //apply lightness

                                float diflc = lightcont - original->L[y][x];
                                diflc *= factorx;
                                transformed->L[y][x] = CLIP(original->L[y][x] + diflc);
                                transformed->a[y][x] = CLIPC(original->a[y][x] * fac);
                                transformed->b[y][x] = CLIPC(original->b[y][x] * fac);


                            }

                            break;
                        }

                        case 0: { // inside selection => full effect, no transition
                            float lumnew = original->L[y][x];

                            if (lp.sens < 75.f) {

                                float lightcont;

                                if ((lp.ligh != 0.f || lp.cont != 0)) {
                                    calclight(lumnew, lp.ligh, lumnew, lightCurveloc);  //replace L-curve
                                    lightcont = lumnew;

                                } else {
                                    lightcont = lumnew;
                                }

                                float fli = 1.f;
                                float flisl = 1.f;
                                float flicur = 1.f;
                                float fac = flicur * (100.f + realchro * falu * falL) / 100.f; //chroma factor transition7

                                float diflc = lightcont * fli * flisl - original->L[y][x];

                                kdiff *= fach * kch;
                                diflc *= kdiff ;
                                transformed->L[y][x] = CLIP(1.f * (original->L[y][x] + diflc));

                                transformed->a[y][x] = CLIPC(original->a[y][x] * fac) ;
                                transformed->b[y][x] = CLIPC(original->b[y][x] * fac);


                            } else {
                                if ((lp.ligh != 0.f || lp.cont != 0)) {
                                    calclight(original->L[y][x], lp.ligh, lumnew, lightCurveloc);
                                }

                                float lightcont = lumnew ;
                                transformed->L[y][x] = CLIP(lightcont);
                                transformed->a[y][x] = CLIPC(original->a[y][x] * facc);
                                transformed->b[y][x] = CLIPC(original->b[y][x] * facc);

                            }
                        }
                    }
                }

            }
        }
    }
    delete origblur;

}

void ImProcFunctions::calc_ref(int sp, LabImage * original, LabImage * transformed, int cx, int cy, int oW, int oH, int sk, double & huerefblur, double & hueref, double & chromaref, double & lumaref, double & sobelref, LUTu & histogram, float &avg)
{
    if (params->locallab.enabled) {
        //always calculate hueref, chromaref, lumaref  before others operations use in normal mode for all modules exceprt denoise
        struct local_params lp;
        calcLocalParams(sp, oW, oH, params->locallab, lp, 0, 0);
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
//        printf("calc avg=%f \n", avg);
//claculate histogram for area selected
#ifdef _OPENMP
        const int numThreads = min(max(transformed->W * transformed->H / (int)histogram.getSize(), 1), omp_get_max_threads());
        #pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
        {
            LUTu lhist16thrloc(histogram.getSize());
            lhist16thrloc.clear();
            histogram.clear();
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int y = 0; y < transformed->H ; y++) //{
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                        int pos = (int)(original->L[y][x]);
                        lhist16thrloc[pos]++;
                    }
                }


#ifdef _OPENMP
            #pragma omp critical
#endif
            histogram += lhist16thrloc;
        }

// double precision for large summations
        double aveA = 0.;
        double aveB = 0.;
        double aveL = 0.;
        double aveChro = 0.;
        double aveAblur = 0.;
        double aveBblur = 0.;
        float avAblur, avBblur;

        double avesobel = 0.;
// int precision for the counters
        int nab = 0;
        int nso = 0;
        int nsb = 0;
// single precision for the result
        float avA, avB, avL;
        int spotSize = 0.88623f * max(1,  lp.cir / sk);  //18
        //O.88623 = sqrt(PI / 4) ==> sqare equal to circle
        int spotSise2; // = 0.88623f * max (1,  lp.cir / sk); //18

        // very small region, don't use omp here
//      printf("cy=%i cx=%i yc=%f xc=%f circ=%i spot=%i tH=%i tW=%i sk=%i\n", cy, cx, lp.yc, lp.xc, lp.cir, spotSize, transformed->H, transformed->W, sk);
//      printf("ymin=%i ymax=%i\n", max (cy, (int) (lp.yc - spotSize)),min (transformed->H + cy, (int) (lp.yc + spotSize + 1)) );
//      printf("xmin=%i xmax=%i\n", max (cx, (int) (lp.xc - spotSize)),min (transformed->W + cx, (int) (lp.xc + spotSize + 1)) );
        LabImage *sobelL;
        LabImage *deltasobelL;
        LabImage *origsob;
        LabImage *origblur = nullptr;
        LabImage *blurorig = nullptr;

        int spotSi = 1 + 2 * max(1,  lp.cir / sk);

        if (spotSi < 5) {
            spotSi = 5;
        }

        spotSise2 = (spotSi - 1) / 2;

        JaggedArray<float> blend3(spotSi, spotSi);

        origsob = new LabImage(spotSi, spotSi);
        sobelL = new LabImage(spotSi, spotSi);
        deltasobelL = new LabImage(spotSi, spotSi);
        bool isdenoise = false;

        if ((lp.noiself > 0.f || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f) && lp.denoiena) {
            isdenoise = true;
        }

        if (isdenoise) {
            origblur = new LabImage(spotSi, spotSi);
            blurorig = new LabImage(spotSi, spotSi);

            for (int y = max(cy, (int)(lp.yc - spotSise2)); y < min(transformed->H + cy, (int)(lp.yc + spotSise2 + 1)); y++) {
                for (int x = max(cx, (int)(lp.xc - spotSise2)); x < min(transformed->W + cx, (int)(lp.xc + spotSise2 + 1)); x++) {
                    int yb = max(cy, (int)(lp.yc - spotSise2));

                    int xb = max(cx, (int)(lp.xc - spotSise2));

                    int z = y - yb;
                    int u = x - xb;
                    origblur->L[z][u] = original->L[y - cy][x - cx];
                    origblur->a[z][u] = original->a[y - cy][x - cx];
                    origblur->b[z][u] = original->b[y - cy][x - cx];

                }
            }

            float radius = 2.f;
            {
                //No omp
                gaussianBlur(origblur->L, blurorig->L, spotSi, spotSi, radius);
                gaussianBlur(origblur->a, blurorig->a, spotSi, spotSi, radius);
                gaussianBlur(origblur->b, blurorig->b, spotSi, spotSi, radius);

            }

            for (int y = 0; y < spotSi; y++) {
                for (int x = 0; x < spotSi; x++) {
                    aveAblur += blurorig->a[y][x];
                    aveBblur += blurorig->b[y][x];
                    nsb++;

                }
            }
        }

        //ref for luma, chroma, hue
        for (int y = max(cy, (int)(lp.yc - spotSize)); y < min(transformed->H + cy, (int)(lp.yc + spotSize + 1)); y++) {
            for (int x = max(cx, (int)(lp.xc - spotSize)); x < min(transformed->W + cx, (int)(lp.xc + spotSize + 1)); x++) {
                aveL += original->L[y - cy][x - cx];
                aveA += original->a[y - cy][x - cx];
                aveB += original->b[y - cy][x - cx];
                //    aveblend += 100.f * blend2[y - cy][x - cx];
                aveChro += sqrtf(SQR(original->b[y - cy][x - cx]) + SQR(original->a[y - cy][x - cx]));
                nab++;
            }
        }

        //ref for sobel
        bool toto = true;

        if (toto) {
            for (int y = max(cy, (int)(lp.yc - spotSise2)); y < min(transformed->H + cy, (int)(lp.yc + spotSise2 + 1)); y++) {
                for (int x = max(cx, (int)(lp.xc - spotSise2)); x < min(transformed->W + cx, (int)(lp.xc + spotSise2 + 1)); x++) {
                    int yb = max(cy, (int)(lp.yc - spotSise2));

                    int xb = max(cx, (int)(lp.xc - spotSise2));

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
                    // avesobel += blend3[y][x];
                    nbs++;
                }

            sobelref = avesobel / nbs;
            printf("sobelref=%f \n", sobelref);
        }

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

            aveAblur = aveAblur / nsb;
            aveBblur = aveBblur / nsb;
            avAblur = aveAblur / 327.68f;
            avBblur = aveBblur / 327.68f;
        }

        if (isdenoise) {
            huerefblur = xatan2f(avBblur, avAblur);
        } else {
            huerefblur = 0.f;
        }

        //    printf("hueblur=%f hue=%f\n", huerefblur, hueref);
        chromaref = aveChro;
        lumaref = avL;

        printf("Calcref => sp=%i huere=%f chromare=%f lumare=%f \n", sp, hueref, chromaref, lumaref);

        if (isdenoise) {
            delete origblur;
            delete blurorig;
        }

        if (lumaref > 95.f) {//to avoid crash
            lumaref = 95.f;
        }
    }
}

void ImProcFunctions::fftw_denoise(int GW, int GH, int max_numblox_W, int min_numblox_W, float **tmp1, array2D<float> *Lin, int numThreads, const struct local_params & lp, int chrom)
{

    fftwf_plan plan_forward_blox[2];
    fftwf_plan plan_backward_blox[2];

    array2D<float> tilemask_in(TS, TS);
    array2D<float> tilemask_out(TS, TS);

    float *Lbloxtmp  = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
    float *fLbloxtmp = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));

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
    const int border = MAX(2, TS / 16);

    for (int i = 0; i < TS; ++i) {
        float i1 = abs((i > TS / 2 ? i - TS + 1 : i));
        float vmask = (i1 < border ? SQR(sin((rtengine::RT_PI * i1) / (2 * border))) : 1.0f);
        float vmask2 = (i1 < 2 * border ? SQR(sin((rtengine::RT_PI * i1) / (2 * border))) : 1.0f);

        for (int j = 0; j < TS; ++j) {
            float j1 = abs((j > TS / 2 ? j - TS + 1 : j));
            tilemask_in[i][j] = (vmask * (j1 < border ? SQR(sin((rtengine::RT_PI * j1) / (2 * border))) : 1.0f)) + epsilon;
            tilemask_out[i][j] = (vmask2 * (j1 < 2 * border ? SQR(sin((rtengine::RT_PI * j1) / (2 * border))) : 1.0f)) + epsilon;

        }
    }


    float *LbloxArray[numThreads];
    float *fLbloxArray[numThreads];



    const int numblox_W = ceil((static_cast<float>(GW)) / (offset)) + 2 * blkrad;
    const int numblox_H = ceil((static_cast<float>(GH)) / (offset)) + 2 * blkrad;


    //residual between input and denoised L channel
    array2D<float> Ldetail(GW, GH, ARRAY2D_CLEAR_DATA);
    array2D<float> totwt(GW, GH, ARRAY2D_CLEAR_DATA); //weight for combining DCT blocks

    for (int i = 0; i < numThreads; ++i) {
        LbloxArray[i]  = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
        fLbloxArray[i] = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
    }

#ifdef _OPENMP
    int masterThread = omp_get_thread_num();
#endif
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        int subThread = masterThread * 1 + omp_get_thread_num();
#else
        int subThread = 0;
#endif
        float blurbuffer[TS * TS] ALIGNED64;
        float *Lblox = LbloxArray[subThread];
        float *fLblox = fLbloxArray[subThread];
        float pBuf[GW + TS + 2 * blkrad * offset] ALIGNED16;
        float nbrwt[TS * TS] ALIGNED64;
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int vblk = 0; vblk < numblox_H; ++vblk) {

            int top = (vblk - blkrad) * offset;
            float * datarow = pBuf + blkrad * offset;

            for (int i = 0; i < TS; ++i) {
                int row = top + i;
                int rr = row;

                if (row < 0) {
                    rr = MIN(-row, GH - 1);
                } else if (row >= GH) {
                    rr = MAX(0, 2 * GH - 2 - row);
                }

                for (int j = 0; j < GW; ++j) {
                    datarow[j] = ((*Lin)[rr][j] - tmp1[rr][j]);
                }

                for (int j = -blkrad * offset; j < 0; ++j) {
                    datarow[j] = datarow[MIN(-j, GW - 1)];
                }

                for (int j = GW; j < GW + TS + blkrad * offset; ++j) {
                    datarow[j] = datarow[MAX(0, 2 * GW - 2 - j)];
                }//now we have a padded data row

                //now fill this row of the blocks with Lab high pass data
                for (int hblk = 0; hblk < numblox_W; ++hblk) {
                    int left = (hblk - blkrad) * offset;
                    int indx = (hblk) * TS; //index of block in malloc

                    if (top + i >= 0 && top + i < GH) {
                        int j;

                        for (j = 0; j < min((-left), TS); ++j) {
                            Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                        }

                        for (; j < min(TS, GW - left); ++j) {
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

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            //fftwf_print_plan (plan_forward_blox);
            if (numblox_W == max_numblox_W) {
                fftwf_execute_r2r(plan_forward_blox[0], Lblox, fLblox);    // DCT an entire row of tiles
            } else {
                fftwf_execute_r2r(plan_forward_blox[1], Lblox, fLblox);    // DCT an entire row of tiles
            }

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // now process the vblk row of blocks for noise reduction

            float params_Ldetail = 0.f;
            float noisevar_Ldetail = 1.f;

            if (chrom == 0) {
                params_Ldetail = min(float(lp.noiseldetail), 99.9f);    // max out to avoid div by zero when using noisevar_Ldetail as divisor
                noisevar_Ldetail = SQR(static_cast<float>(SQR(100. - params_Ldetail) + 50.*(100. - params_Ldetail)) * TS * 0.5f);
            } else if (chrom == 1) {
                params_Ldetail = min(float(lp.noisechrodetail), 99.9f);
                noisevar_Ldetail = 100.f * pow((static_cast<float>(SQR(100. - params_Ldetail) + 50.*(100. - params_Ldetail)) * TS * 0.5f), 2);//to test ???
            }

            //   float noisevar_Ldetail = SQR(static_cast<float>(SQR(100. - params_Ldetail) + 50.*(100. - params_Ldetail)) * TS * 0.5f);



            for (int hblk = 0; hblk < numblox_W; ++hblk) {
                ImProcFunctions::RGBtile_denoise(fLblox, hblk, noisevar_Ldetail, nbrwt, blurbuffer);
            }//end of horizontal block loop

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            //now perform inverse FT of an entire row of blocks
            if (numblox_W == max_numblox_W) {
                fftwf_execute_r2r(plan_backward_blox[0], fLblox, Lblox);    //for DCT
            } else {
                fftwf_execute_r2r(plan_backward_blox[1], fLblox, Lblox);    //for DCT
            }

            int topproc = (vblk - blkrad) * offset;

            //add row of blocks to output image tile
            ImProcFunctions::RGBoutput_tile_row(Lblox, Ldetail, tilemask_out, GH, GW, topproc);

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        }//end of vertical block loop

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    }
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef _OPENMP

    #pragma omp parallel for
#endif

    for (int i = 0; i < GH; ++i) {
        for (int j = 0; j < GW; ++j) {
            //may want to include masking threshold for large hipass data to preserve edges/detail
            tmp1[i][j] += Ldetail[i][j] / totwt[i][j]; //note that labdn initially stores the denoised hipass data
        }
    }


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

void ImProcFunctions::Lab_Local(int call, int sp, float** shbuffer, LabImage * original, LabImage * transformed, LabImage * reserved, int cx, int cy, int oW, int oH, int sk,
                                const LocretigainCurve & locRETgainCcurve, LUTf & lllocalcurve, bool & locallutili, const LocLHCurve & loclhCurve,  const LocHHCurve & lochhCurve, const LocCCmaskCurve & locccmasCurve, bool & lcmasutili, const  LocLLmaskCurve & locllmasCurve, bool & llmasutili, const  LocHHmaskCurve & lochhmasCurve, bool &lhmasutili, const LocCCmaskexpCurve & locccmasexpCurve, bool &lcmasexputili, const  LocLLmaskexpCurve & locllmasexpCurve, bool &llmasexputili, const  LocHHmaskexpCurve & lochhmasexpCurve, bool & lhmasexputili,
                                bool & LHutili, bool & HHutili, LUTf & cclocalcurve, bool & localskutili, LUTf & sklocalcurve, bool & localexutili, LUTf & exlocalcurve, LUTf & hltonecurveloc, LUTf & shtonecurveloc, LUTf & tonecurveloc, LUTf & lightCurveloc, double & huerefblur, double & hueref, double & chromaref, double & lumaref, double & sobelref, int llColorMask, int llExpMask)
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

        struct local_params lp;
        calcLocalParams(sp, oW, oH, params->locallab, lp, llColorMask, llExpMask);

        const float radius = lp.rad / (sk * 1.4f); //0 to 70 ==> see skip
        int strred = (lp.strucc - 1);

        if (strred > 1) {
            strred = 1;
        }

        float radiussob = strred / (sk * 1.4f);
        double ave = 0.;
        int n = 0;
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
        }

        printf("call= %i sp=%i hueref=%f chromaref=%f lumaref=%f\n", call, sp, hueref, chromaref, lumaref);
//        struct local_contra lco;

// we must here detect : general case, skin, sky,...foliages ???
// delta dhue, luminance and chroma

//sens, sensh, senscb, sensbn, senstm;
        constexpr float ared = (rtengine::RT_PI - 0.05f) / 100.f;

        constexpr float bred = 0.05f;
        float dhue = ared * lp.sens + bred; //delta hue vibr

        float dhuev = ared * lp.sensv + bred; //delta hue vibr

        float dhueex = ared * lp.sensex + bred; //delta hue exp

        float dhueret = ared * lp.sensh + bred; //delta hue retinex

        float dhuebn = ared * lp.sensbn + bred; //delta hue blur

        float dhuetm = ared * lp.senstm + bred; //delta hue tone map

        float dhuesha = ared * lp.senssha + bred; //delta hue sharp

        float dhuelc = ared * lp.senslc + bred; //delta hue local contrast

        float dhuesf = ared * lp.senssf + bred; //delta hue soft light

        float dhuecb = ared * lp.senscb + bred; //delta hue cbdl

        float dhueexclu = ared * lp.sensexclu + bred; //delta hue exclude

        float dhueden = ared * lp.sensden + bred; //delta hue lght chroma

        float moddE = 2.5f;
        float powdE = 6.f;



        if (lp.excmet == 1  && call <= 3) {
            LabImage *deltasobelL = nullptr;
            LabImage *tmpsob = nullptr;
            LabImage *bufsob = nullptr;
            LabImage *bufreserv = nullptr;
            LabImage *bufexclu = nullptr;
            float *origBuffer = nullptr;

            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;
            int begy = lp.yc - lp.lyT;
            int begx = lp.xc - lp.lxL;
            int yEn = lp.yc + lp.ly;
            int xEn = lp.xc + lp.lx;
            bufsob = new LabImage(bfw, bfh);
            bufreserv = new LabImage(bfw, bfh);
            JaggedArray<float> buflight(bfw, bfh);
            JaggedArray<float> bufchro(bfw, bfh);

            float *orig[bfh] ALIGNED16;
            origBuffer = new float[bfh * bfw];

            for (int i = 0; i < bfh; i++) {
                orig[i] = &origBuffer[i * bfw];
            }

            bufexclu = new LabImage(bfw, bfh);


#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int ir = 0; ir < bfh; ir++) //fill with 0
                for (int jr = 0; jr < bfw; jr++) {
                    bufsob->L[ir][jr] = 0.f;
                    bufexclu->L[ir][jr] = 0.f;
                    bufexclu->a[ir][jr] = 0.f;
                    bufexclu->b[ir][jr] = 0.f;
                    buflight[ir][jr] = 0.f;
                    bufchro[ir][jr] = 0.f;
                    bufreserv->L[ir][jr] = 0.f;
                    bufreserv->a[ir][jr] = 0.f;
                    bufreserv->b[ir][jr] = 0.f;
                }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int y = 0; y < transformed->H ; y++) //{
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                        bufreserv->L[loy - begy][lox - begx] = reserved->L[y][x];
                        bufreserv->a[loy - begy][lox - begx] = reserved->a[y][x];
                        bufreserv->b[loy - begy][lox - begx] = reserved->b[y][x];
                        bufexclu->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                        bufexclu->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                        bufexclu->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas

                    }
                }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int y = 0; y < transformed->H ; y++) //{
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                        //  bufsob->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                        bufsob->L[loy - begy][lox - begx] = reserved->L[y][x];//fill square buffer with datas

                    }
                }

            tmpsob = new LabImage(bfw, bfh);
            deltasobelL = new LabImage(bfw, bfh);
            SobelCannyLuma(tmpsob->L, bufsob->L, bfw, bfh, radiussob);
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int ir = 0; ir < bfh; ir++) //fill with 0
                for (int jr = 0; jr < bfw; jr++) {
                    deltasobelL->L[ir][jr] = 1.f;
                }


            //then restore non modified area

            //TODO then use instead of others modifications Color and Light, Blur, etc.
            float hueplus = hueref + dhueexclu;
            float huemoins = hueref - dhueexclu;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhueexclu - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhueexclu + 2.f * rtengine::RT_PI;
            }

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int ir = 0; ir < bfh; ir++)
                for (int jr = 0; jr < bfw; jr++) {
                    float rL;
                    rL = (bufreserv->L[ir][jr] - bufexclu->L[ir][jr]) / 327.68f;
                    buflight[ir][jr] = rL ;


                }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int ir = 0; ir < bfh; ir += 1)
                for (int jr = 0; jr < bfw; jr += 1) {
                    orig[ir][jr] = sqrt(SQR(bufexclu->a[ir][jr]) + SQR(bufexclu->b[ir][jr]));
                }


#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int ir = 0; ir < bfh; ir++)
                for (int jr = 0; jr < bfw; jr++) {
                    float rch;
                    rch = CLIPRET((sqrt((SQR(bufreserv->a[ir][jr]) + SQR(bufreserv->b[ir][jr]))) - orig[ir][jr])) / 327.68f;
                    bufchro[ir][jr] = rch ;
                }

            /*
                        if (lp.sensexclu < 15.f) {
                            JaggedArray<float> blend(bfw, bfh);
                            float contrastf = lp.sensexclu / 100.f;

                            buildBlendMask(bufreserv->L, blend, bfw, bfh, contrastf, 1.f, true);

            #ifdef _OPENMP
                            #pragma omp parallel for
            #endif

                            for (int i = 0; i < bfh; ++i) {
                                for (int j = 0; j < bfw; ++j) {
                                    if (blend[i][j] > 0.05f) {
                                        blend[i][j] = 1.22361f * sqrt(blend[i][j]) - 0.22361f;
                                    }

                                    bufreserv->L[i][j] = intp(blend[i][j], bufexclu->L[i][j], bufreserv->L[i][j]);
                                    bufreserv->a[i][j] = intp(blend[i][j], bufexclu->a[i][j], bufreserv->a[i][j]);
                                    bufreserv->b[i][j] = intp(blend[i][j], bufexclu->b[i][j], bufreserv->b[i][j]);

                                }
                            }
                        }

            #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
            #endif

                        for (int y = 0; y < transformed->H ; y++) //{
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    original->L[y][x] = bufreserv->L[loy - begy][lox - begx];
                                    original->a[y][x] = bufreserv->a[loy - begy][lox - begx];
                                    original->b[y][x] = bufreserv->b[loy - begy][lox - begx];
                                }
                            }
            */
            Exclude_Local(1, deltasobelL->L, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, bufreserv, reserved, cx, cy, sk);


            delete deltasobelL;
            delete tmpsob;


            delete bufexclu;
            delete [] origBuffer;

            delete bufreserv;
            delete bufsob;

        }

//Blur and noise

        if (((radius >= 1.5 * GAUSS_SKIP && lp.rad > 1.) || lp.stren > 0.1)  && lp.blurena) { // radius < GAUSS_SKIP means no gauss, just copy of original image
            LabImage *tmp1 = nullptr;
            LabImage *tmp2 = nullptr;
            LabImage *bufgb = nullptr;
            float *origBuffer = nullptr;
            //       LabImage *deltasobelL = nullptr;
            int GW = transformed->W;
            int GH = transformed->H;
            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;
            JaggedArray<float> buflight(bfw, bfh);
            JaggedArray<float> bufchro(bfw, bfh);

            float *orig[bfh] ALIGNED16;

            if (call <= 3  && lp.blurmet != 1) {
                bufgb = new LabImage(bfw, bfh);

                origBuffer = new float[bfh * bfw];

                for (int i = 0; i < bfh; i++) {
                    orig[i] = &origBuffer[i * bfw];
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
                            bufgb->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufgb->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufgb->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas
                        }
                    }

                tmp1 = new LabImage(bfw, bfh);



                if (lp.blurmet == 2) {
                    tmp2 = new LabImage(transformed->W, transformed->H);
#ifdef _OPENMP
                    #pragma omp parallel
#endif
                    {
                        gaussianBlur(original->L, tmp2->L, GW, GH, radius);
                        gaussianBlur(original->a, tmp2->a, GW, GH, radius);
                        gaussianBlur(original->b, tmp2->b, GW, GH, radius);
                    }

                }



#ifdef _OPENMP
                #pragma omp parallel
#endif

                {
                    gaussianBlur(bufgb->L, tmp1->L, bfw, bfh, radius);
                    gaussianBlur(bufgb->a, tmp1->a, bfw, bfh, radius);
                    gaussianBlur(bufgb->b, tmp1->b, bfw, bfh, radius);
                }


            } else {
                tmp1 = new LabImage(transformed->W, transformed->H);;

#ifdef _OPENMP
                #pragma omp parallel
#endif
                {
                    gaussianBlur(original->L, tmp1->L, GW, GH, radius);
                    gaussianBlur(original->a, tmp1->a, GW, GH, radius);
                    gaussianBlur(original->b, tmp1->b, GW, GH, radius);

                }
            }

            if (lp.stren > 0.1f) {
                if (lp.blurmet <= 1) {

                    float mean = 0.f;//0 best result
                    float variance = lp.stren ; //(double) SQR(lp.stren)/sk;
                    addGaNoise(tmp1, tmp1, mean, variance, sk) ;
                }

            }

            if (lp.blurmet != 1) { //blur and noise (center)
                float hueplus = hueref + dhuebn;
                float huemoins = hueref - dhuebn;

                if (hueplus > rtengine::RT_PI) {
                    hueplus = hueref + dhuebn - 2.f * rtengine::RT_PI;
                }

                if (huemoins < -rtengine::RT_PI) {
                    huemoins = hueref - dhuebn + 2.f * rtengine::RT_PI;
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++)
                    for (int jr = 0; jr < bfw; jr++) {
                        float rL;
                        rL = CLIPRET((tmp1->L[ir][jr] - bufgb->L[ir][jr]) / 328.f);
                        buflight[ir][jr] = rL;

                    }

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int ir = 0; ir < bfh; ir += 1)
                    for (int jr = 0; jr < bfw; jr += 1) {
                        orig[ir][jr] = sqrt(SQR(bufgb->a[ir][jr]) + SQR(bufgb->b[ir][jr]));
                    }


#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++)
                    for (int jr = 0; jr < bfw; jr++) {
                        float rch;
                        rch = CLIPRET((sqrt((SQR(tmp1->a[ir][jr]) + SQR(tmp1->b[ir][jr]))) - orig[ir][jr]) / 328.f);
                        bufchro[ir][jr] = rch;
                    }

                BlurNoise_Local(moddE, powdE, call, tmp1, tmp2, buflight, bufchro, hueplus, huemoins, hueref, dhuebn, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

            } else {

                InverseBlurNoise_Local(lp, original, transformed, tmp1, cx, cy);

            }

            if (call <= 3  && lp.blurmet != 1) {

                delete bufgb;

                delete [] origBuffer;


            }

            delete tmp1;

            if (lp.blurmet == 2) {
                delete tmp2;
            }

        }


        //local impulse
        if ((lp.bilat > 0.f) && lp.denoiena) {
            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;

            LabImage *bufwv = nullptr;
            float hueplus = huerefblur + dhueden;
            float huemoins = huerefblur - dhueden;

            if (hueplus > rtengine::RT_PI) {
                hueplus = huerefblur + dhueden - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = huerefblur - dhueden + 2.f * rtengine::RT_PI;
            }

            if (call == 2) {//simpleprocess
                bufwv = new LabImage(bfw, bfh); //buffer for data in zone limit


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
                            bufwv->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufwv->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufwv->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas
                        }

                    }
            } else {//dcrop.cc

                int GH = transformed->H;
                int GW = transformed->W;

                bufwv = new LabImage(GW, GH);
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int ir = 0; ir < GH; ir++)
                    for (int jr = 0; jr < GW; jr++) {
                        bufwv->L[ir][jr] = original->L[ir][jr];
                        bufwv->a[ir][jr] = original->a[ir][jr];
                        bufwv->b[ir][jr] = original->b[ir][jr];
                    }


            } //end dcrop

            double thr = (float) lp.bilat / 20.0;

            if (bfh > 8 && bfw > 8) {
                ImProcFunctions::impulse_nr(bufwv, thr);
            }

            LabImage tmp1(bufwv->W, bufwv->H);
            //copy bufwv to tmp1 to use same algo for Denoise_local and DeNoise_Local_imp

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int ir = 0; ir < bufwv->H; ir++)
                for (int jr = 0; jr < bufwv->W; jr++) {
                    tmp1.L[ir][jr] = bufwv->L[ir][jr];
                    tmp1.a[ir][jr] = bufwv->a[ir][jr];
                    tmp1.b[ir][jr] = bufwv->b[ir][jr];
                }



            //    DeNoise_Local_imp(call, lp,  levred, hueplus, huemoins, huerefblur, dhueden, original, transformed, bufwv, cx, cy, sk);
            DeNoise_Local(call, lp,  levred, hueplus, huemoins, huerefblur, dhueden, original, transformed, tmp1, cx, cy, sk);

            delete bufwv;
        }

//local denoise
        //all these variables are to prevent use of denoise when non necessary
        // but with qualmet = 2 (default for best quality) we must denoise chroma with little values to prevent artifacts due to variations of Hue
        // but if user select volontary denoise, it is that choice the good (prioritary)
        bool execdenoi = false ;
        bool execcolor = (lp.chro != 0.f || lp.ligh != 0.f || lp.cont != 0); // only if one slider ore more is engaged
        bool execbdl = (lp.mulloc[0] != 1.f || lp.mulloc[1] != 1.f || lp.mulloc[2] != 1.f || lp.mulloc[3] != 1.f || lp.mulloc[4] != 1.f) ;//only if user want cbdl
        execdenoi = noiscfactiv && ((lp.colorena && execcolor) || (lp.tonemapena && lp.strengt != 0.f) || (lp.cbdlena && execbdl) || (lp.sfena && lp.strng > 0.f) || (lp.lcena && lp.lcamount > 0.f) || (lp.sharpena && lp.shrad > 0.42) || (lp.retiena  && lp.str > 0.f)  || (lp.exposena && lp.expcomp != 0.f)  || (lp.expvib  && lp.past != 0.f));

        if (((lp.noiself > 0.f || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f) && lp.denoiena) || execdenoi) {  // sk == 1 ??
            StopWatch Stop1("locallab Denoise called");
            MyMutex::MyLock lock(*fftwMutex);

            if (lp.noisecf >= 0.1f || lp.noisecc >= 0.1f) {
                noiscfactiv = false;
                levred = 7;
            }

            float hueplus = huerefblur + dhueden;
            float huemoins = huerefblur - dhueden;

            if (hueplus > rtengine::RT_PI) {
                hueplus = huerefblur + dhueden - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = huerefblur - dhueden + 2.f * rtengine::RT_PI;
            }

#ifdef _OPENMP
            const int numThreads = omp_get_max_threads();
#else
            const int numThreads = 1;

#endif

            if (call == 1) {
                //printf("OK sk=1\n");
                LabImage tmp1(transformed->W, transformed->H);
                LabImage tmp2(transformed->W, transformed->H);
                tmp2.clear();

                array2D<float> *Lin = nullptr;
                array2D<float> *Ain = nullptr;
                array2D<float> *Bin = nullptr;


                int GW = transformed->W;
                int GH = transformed->H;
                int max_numblox_W = ceil((static_cast<float>(GW)) / (offset)) + 2 * blkrad;
                // calculate min size of numblox_W.
                int min_numblox_W = ceil((static_cast<float>(GW)) / (offset)) + 2 * blkrad;


                for (int ir = 0; ir < GH; ir++)
                    for (int jr = 0; jr < GW; jr++) {
                        tmp1.L[ir][jr] = original->L[ir][jr];
                        tmp1.a[ir][jr] = original->a[ir][jr];
                        tmp1.b[ir][jr] = original->b[ir][jr];
                    }

                int DaubLen = 6;

                int levwavL = levred;
                int skip = 1;

                wavelet_decomposition Ldecomp(tmp1.L[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, DaubLen);
                wavelet_decomposition adecomp(tmp1.a[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, DaubLen);
                wavelet_decomposition bdecomp(tmp1.b[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, DaubLen);

                float madL[8][3];
                int edge = 2;

                if (!Ldecomp.memoryAllocationFailed) {
                    #pragma omp parallel for collapse(2) schedule(dynamic,1)

                    for (int lvl = 0; lvl < levred; lvl++) {
                        for (int dir = 1; dir < 4; dir++) {
                            int Wlvl_L = Ldecomp.level_W(lvl);
                            int Hlvl_L = Ldecomp.level_H(lvl);

                            float ** WavCoeffs_L = Ldecomp.level_coeffs(lvl);

                            madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                        }
                    }

                    float vari[levred];

                    if (levred == 7) {
                        edge = 2;
                        vari[0] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[1] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));

                        vari[3] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[4] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[5] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[6] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    } else if (levred == 4) {
                        edge = 3;
                        vari[0] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[1] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[3] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));

                    }

                    if ((lp.noiself >= 0.1f ||  lp.noiselc >= 0.1f)) {
                        float kr3 = 0.f;
                        float kr4 = 0.f;
                        float kr5 = 0.f;

                        if (lp.noiselc < 30.f) {
                            kr3 = 0.f;
                            kr4 = 0.f;
                            kr5 = 0.f;
                        } else if (lp.noiselc < 50.f) {
                            kr3 = 0.5f;
                            kr4 = 0.3f;
                            kr5 = 0.2f;
                        } else if (lp.noiselc < 70.f) {
                            kr3 = 0.7f;
                            kr4 = 0.5f;
                            kr5 = 0.3f;
                        } else {
                            kr3 = 1.f;
                            kr4 = 1.f;
                            kr5 = 1.f;
                        }

                        vari[0] = max(0.0001f, vari[0]);
                        vari[1] = max(0.0001f, vari[1]);
                        vari[2] = max(0.0001f, vari[2]);
                        vari[3] = max(0.0001f, kr3 * vari[3]);

                        if (levred == 7) {
                            vari[4] = max(0.0001f, kr4 * vari[4]);
                            vari[5] = max(0.0001f, kr5 * vari[5]);
                            vari[6] = max(0.0001f, kr5 * vari[6]);
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
                        #pragma omp parallel for

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


                        if (lp.noiselc < 1.f) {
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
                    noisecfr = lp.noisecf * ((100.f + lp.adjch) / 10.f);
                    noiseccr = lp.noisecc + ((100.f + lp.adjch) / 10.f);
                }

                float noisecfb = lp.noisecf;
                float noiseccb = lp.noisecc;

                if (lp.adjch < 0.f) {
                    noisecfb = lp.noisecf * ((100.f - lp.adjch) / 10.f);
                    noiseccb = lp.noisecc * ((100.f - lp.adjch) / 10.f);
                }


                if (noisecfr < 0.f) {
                    noisecfr = 0.0001f;
                }

                if (noiseccr < 0.f) {
                    noiseccr = 0.0001f;
                }

                if (noisecfb < 0.f) {
                    noisecfb = 0.0001f;
                }

                if (noiseccb < 0.f) {
                    noiseccb = 0.0001f;
                }

                if (!adecomp.memoryAllocationFailed && !bdecomp.memoryAllocationFailed) {

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

                    if ((lp.noisecf >= 0.1f ||  lp.noisecc >= 0.1f  || noiscfactiv)) {
                        float minic = 0.0001f;

                        if (noiscfactiv) {
                            minic = 0.1f;//only for artifact shape detection
                        }

                        float k1 = 0.f;
                        float k2 = 0.f;
                        float k3 = 0.f;

                        if (lp.noisecf < 0.2f) {
                            k1 = 0.f;
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
                        } else if (lp.noisecf < 10.f) {
                            k1 = 0.85f;
                            k2 = 0.7f;
                            k3 = 0.6f;
                        } else if (lp.noisecf < 20.f) {
                            k1 = 0.9f;
                            k2 = 0.8f;
                            k3 = 0.7f;
                        } else if (lp.noisecf < 50.f) {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 0.9f;

                        } else {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 1.f;
                        }

                        variC[0] = max(minic, variC[0]);
                        variC[1] = max(minic, k1 * variC[1]);
                        variC[2] = max(minic, k2 * variC[2]);
                        variC[3] = max(minic, k3 * variC[3]);

                        variCb[0] = max(minic, variCb[0]);
                        variCb[1] = max(minic, k1 * variCb[1]);
                        variCb[2] = max(minic, k2 * variCb[2]);
                        variCb[3] = max(minic, k3 * variCb[3]);

                        if (levred == 7) {
                            float k4 = 0.f;
                            float k5 = 0.f;
                            float k6 = 0.f;

                            if (lp.noisecc == 0.1f) {
                                k4 = 0.f;
                                k5 = 0.0f;
                            } else if (lp.noisecc < 0.2f) {
                                k4 = 0.1f;
                                k5 = 0.0f;
                            } else if (lp.noisecc < 0.5f) {
                                k4 = 0.15f;
                                k5 = 0.0f;
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


                            variC[4] = max(0.0001f, k4 * variC[4]);
                            variC[5] = max(0.0001f, k5 * variC[5]);
                            variCb[4] = max(0.0001f, k4 * variCb[4]);
                            variCb[5] = max(0.0001f, k5 * variCb[5]);

                            if (lp.noisecc < 4.f) {
                                k6 = 0.f;
                            } else if (lp.noisecc < 5.f) {
                                k6 = 0.4f;
                            } else if (lp.noisecc < 6.f) {
                                k6 = 0.7f;
                            } else {
                                k6 = 1.f;
                            }

                            variC[6] = max(0.0001f, k6 * variC[6]);
                            variCb[6] = max(0.0001f, k6 * variCb[6]);

                        }

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
                        #pragma omp parallel for

#endif

                        for (int ir = 0; ir < GH; ir++)
                            for (int jr = 0; jr < GW; jr++) {
                                float cN = sqrt(SQR(tmp1.a[ir][jr]) + SQR(tmp1.b[ir][jr]));

                                if (cN < seuil) {
                                    noisevarchrom[(ir >> 1)*GW2 + (jr >> 1)] =  nvch;
                                } else if (cN < seuil2) {
                                    noisevarchrom[(ir >> 1)*GW2 + (jr >> 1)] = ac * cN + bc;
                                } else {
                                    noisevarchrom[(ir >> 1)*GW2 + (jr >> 1)] =  nvcl;
                                }
                            }


                        float noisevarab_r = 100.f; //SQR(lp.noisecc / 10.0);

                        if (lp.noisecc < 0.1f)  {
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

                if (!Ldecomp.memoryAllocationFailed) {
                    Lin = new array2D<float>(GW, GH);
#ifdef _OPENMP
                    #pragma omp parallel for

#endif

                    for (int i = 0; i < GH; ++i) {
                        for (int j = 0; j < GW; ++j) {
                            (*Lin)[i][j] = tmp1.L[i][j];
                        }
                    }

                    Ldecomp.reconstruct(tmp1.L[0]);
                }

                if (!Ldecomp.memoryAllocationFailed) {
                    if ((lp.noiself >= 0.1f ||  lp.noiselc >= 0.1f)  && levred == 7) {
                        fftw_denoise(GW, GH, max_numblox_W, min_numblox_W, tmp1.L, Lin,  numThreads, lp, 0);
                    }
                }

                if (!adecomp.memoryAllocationFailed) {
                    Ain = new array2D<float>(GW, GH);
#ifdef _OPENMP
                    #pragma omp parallel for

#endif

                    for (int i = 0; i < GH; ++i) {
                        for (int j = 0; j < GW; ++j) {
                            (*Ain)[i][j] = tmp1.a[i][j];
                        }
                    }

                    adecomp.reconstruct(tmp1.a[0]);
                }


                if (!adecomp.memoryAllocationFailed) {
                    if ((lp.noisecf >= 0.1f ||  lp.noisecc >= 0.1f)) {
                        if (lp.noisechrodetail > 1000) { //to avoid all utilisation
                            fftw_denoise(GW, GH, max_numblox_W, min_numblox_W, tmp1.a, Ain,  numThreads, lp, 1);
                        }
                    }




                }


                if (!bdecomp.memoryAllocationFailed) {

                    Bin = new array2D<float>(GW, GH);
#ifdef _OPENMP
                    #pragma omp parallel for

#endif

                    for (int i = 0; i < GH; ++i) {
                        for (int j = 0; j < GW; ++j) {
                            (*Bin)[i][j] = tmp1.b[i][j];
                        }
                    }

                    bdecomp.reconstruct(tmp1.b[0]);
                }


                if (!bdecomp.memoryAllocationFailed) {
                    if ((lp.noisecf >= 0.1f ||  lp.noisecc >= 0.1f)) {
                        if (lp.noisechrodetail > 1000) {//to avoid all utilisation

                            fftw_denoise(GW, GH, max_numblox_W, min_numblox_W, tmp1.b, Bin,  numThreads, lp, 1);
                        }
                    }

                }


                DeNoise_Local(call, lp,  levred, hueplus, huemoins, huerefblur, dhueden, original, transformed, tmp1, cx, cy, sk);

            } else if (call == 2 /* || call == 1 || call == 3 */) { //simpleprocess

                int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                int bfw = int (lp.lx + lp.lxL) + del;
                LabImage bufwv(bfw, bfh);
                bufwv.clear(true);
                array2D<float> *Lin = nullptr;
                //  array2D<float> *Ain = nullptr;
                //  array2D<float> *Bin = nullptr;

                int max_numblox_W = ceil((static_cast<float>(bfw)) / (offset)) + 2 * blkrad;
                // calculate min size of numblox_W.
                int min_numblox_W = ceil((static_cast<float>(bfw)) / (offset)) + 2 * blkrad;
                // these are needed only for creation of the plans and will be freed before entering the parallel loop


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
                wavelet_decomposition Ldecomp(bufwv.L[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, DaubLen);
                wavelet_decomposition adecomp(bufwv.a[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, DaubLen);
                wavelet_decomposition bdecomp(bufwv.b[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, DaubLen);

                float madL[8][3];
                int edge = 2;

                if (!Ldecomp.memoryAllocationFailed) {
                    #pragma omp parallel for collapse(2) schedule(dynamic,1)

                    for (int lvl = 0; lvl < levred; lvl++) {
                        for (int dir = 1; dir < 4; dir++) {
                            int Wlvl_L = Ldecomp.level_W(lvl);
                            int Hlvl_L = Ldecomp.level_H(lvl);

                            float ** WavCoeffs_L = Ldecomp.level_coeffs(lvl);

                            madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                        }
                    }

                    float vari[levred];

                    if (levred == 7) {
                        edge = 2;
                        vari[0] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[1] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));

                        vari[3] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[4] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[5] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[6] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    } else if (levred == 4) {
                        edge = 3;
                        vari[0] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[1] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[3] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));

                    }


                    if ((lp.noiself >= 0.1f ||  lp.noiselc >= 0.1f)) {
                        float kr3 = 0.f;
                        float kr4 = 0.f;
                        float kr5 = 0.f;

                        if (lp.noiselc < 30.f) {
                            kr3 = 0.f;
                            kr4 = 0.f;
                            kr5 = 0.f;
                        } else if (lp.noiselc < 50.f) {
                            kr3 = 0.5f;
                            kr4 = 0.3f;
                            kr5 = 0.2f;
                        } else if (lp.noiselc < 70.f) {
                            kr3 = 0.7f;
                            kr4 = 0.5f;
                            kr5 = 0.3f;
                        } else {
                            kr3 = 1.f;
                            kr4 = 1.f;
                            kr5 = 1.f;
                        }

                        vari[0] = max(0.0001f, vari[0]);
                        vari[1] = max(0.0001f, vari[1]);
                        vari[2] = max(0.0001f, vari[2]);
                        vari[3] = max(0.0001f, kr3 * vari[3]);

                        if (levred == 7) {
                            vari[4] = max(0.0001f, kr4 * vari[4]);
                            vari[5] = max(0.0001f, kr5 * vari[5]);
                            vari[6] = max(0.0001f, kr5 * vari[6]);
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
                        #pragma omp parallel for

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

                        if (lp.noiselc < 1.f) {
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
                    noisecfr = lp.noisecf * ((100.f + lp.adjch) / 10.f);
                    noiseccr = lp.noisecc + ((100.f + lp.adjch) / 10.f);
                }

                float noisecfb = lp.noisecf;
                float noiseccb = lp.noisecc;

                if (lp.adjch < 0.f) {
                    noisecfb = lp.noisecf * ((100.f - lp.adjch) / 10.f);
                    noiseccb = lp.noisecc * ((100.f - lp.adjch) / 10.f);
                }


                if (noisecfr < 0.f) {
                    noisecfr = 0.0001f;
                }

                if (noiseccr < 0.f) {
                    noiseccr = 0.0001f;
                }

                if (noisecfb < 0.f) {
                    noisecfb = 0.0001f;
                }

                if (noiseccb < 0.f) {
                    noiseccb = 0.0001f;
                }


                if (!adecomp.memoryAllocationFailed && !bdecomp.memoryAllocationFailed) {

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


                    if ((lp.noisecf >= 0.1f ||  lp.noisecc >= 0.1f  || noiscfactiv)) {
                        float minic = 0.0001f;

                        if (noiscfactiv) {
                            minic = 0.1f;//only for artifact shape detection
                        }

                        float k1 = 0.f;
                        float k2 = 0.f;
                        float k3 = 0.f;

                        if (lp.noisecf < 0.2f) {
                            k1 = 0.f;
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
                        } else if (lp.noisecf < 10.f) {
                            k1 = 0.85f;
                            k2 = 0.7f;
                            k3 = 0.6f;
                        } else if (lp.noisecf < 20.f) {
                            k1 = 0.9f;
                            k2 = 0.8f;
                            k3 = 0.7f;
                        } else if (lp.noisecf < 50.f) {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 0.9f;

                        } else {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 1.f;
                        }

                        variC[0] = max(minic, variC[0]);
                        variC[1] = max(minic, k1 * variC[1]);
                        variC[2] = max(minic, k2 * variC[2]);
                        variC[3] = max(minic, k3 * variC[3]);

                        variCb[0] = max(minic, variCb[0]);
                        variCb[1] = max(minic, k1 * variCb[1]);
                        variCb[2] = max(minic, k2 * variCb[2]);
                        variCb[3] = max(minic, k3 * variCb[3]);

                        if (levred == 7) {
                            float k4 = 0.f;
                            float k5 = 0.f;
                            float k6 = 0.f;

                            if (lp.noisecc == 0.1f) {
                                k4 = 0.f;
                                k5 = 0.0f;
                            } else if (lp.noisecc < 0.2f) {
                                k4 = 0.1f;
                                k5 = 0.0f;
                            } else if (lp.noisecc < 0.5f) {
                                k4 = 0.15f;
                                k5 = 0.0f;
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


                            variC[4] = max(0.0001f, k4 * variC[4]);
                            variC[5] = max(0.0001f, k5 * variC[5]);
                            variCb[4] = max(0.0001f, k4 * variCb[4]);
                            variCb[5] = max(0.0001f, k5 * variCb[5]);

                            if (lp.noisecc < 4.f) {
                                k6 = 0.f;
                            } else if (lp.noisecc < 5.f) {
                                k6 = 0.4f;
                            } else if (lp.noisecc < 6.f) {
                                k6 = 0.7f;
                            } else {
                                k6 = 1.f;
                            }

                            variC[6] = max(0.0001f, k6 * variC[6]);
                            variCb[6] = max(0.0001f, k6 * variCb[6]);
                        }

                        float* noisevarchrom = new float[bfh * bfw];
                        int bfw2 = (bfw + 1) / 2;
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
                        #pragma omp parallel for

#endif

                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                float cN = sqrt(SQR(bufwv.a[ir][jr]) + SQR(bufwv.b[ir][jr]));

                                if (cN < seuil) {
                                    noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] = nvch;
                                } else if (cN < seuil2) {
                                    noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] = ac * cN + bc;
                                } else {
                                    noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] = nvcl;
                                }
                            }

                        float noisevarab_r = 100.f; //SQR(lp.noisecc / 10.0);


                        if (lp.noisecc < 0.1f)  {
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

                if (!Ldecomp.memoryAllocationFailed) {
                    Lin = new array2D<float>(bfw, bfh);

#ifdef _OPENMP
                    #pragma omp parallel for

#endif

                    for (int i = 0; i < bfh; ++i) {
                        for (int j = 0; j < bfw; ++j) {
                            (*Lin)[i][j] = bufwv.L[i][j];
                        }
                    }

                    Ldecomp.reconstruct(bufwv.L[0]);
                }


                if (!Ldecomp.memoryAllocationFailed) {


                    if ((lp.noiself >= 0.1f ||  lp.noiselc >= 0.1f) && levred == 7) {
                        fftw_denoise(bfw, bfh, max_numblox_W, min_numblox_W, bufwv.L, Lin,  numThreads, lp, 0);
                    }
                }


                if (!adecomp.memoryAllocationFailed) {
                    /*
                    //   Ain = new array2D<float>(bfw, bfh);
                    #ifdef _OPENMP
                    #pragma omp parallel for

                    #endif

                    for (int i = 0; i < bfh; ++i) {
                        for (int j = 0; j < bfw; ++j) {
                            (*Ain)[i][j] = bufwv.a[i][j];
                        }
                    }
                    */
                    adecomp.reconstruct(bufwv.a[0]);
                }

                /*
                                if (!adecomp.memoryAllocationFailed) {


                                    if ((lp.noisecf >= 0.1f ||  lp.noisecc >= 0.1f)) {
                                        //   fftw_denoise(bfw, bfh, max_numblox_W, min_numblox_W, bufwv.a, Ain,  numThreads, lp, 1);

                                    }
                                }
                */

                if (!bdecomp.memoryAllocationFailed) {
                    /*
                    //             Bin = new array2D<float>(bfw, bfh);
                    #ifdef _OPENMP
                    #pragma omp parallel for

                    #endif

                    for (int i = 0; i < bfh; ++i) {
                        for (int j = 0; j < bfw; ++j) {
                            (*Bin)[i][j] = bufwv.b[i][j];
                        }
                    }
                    */
                    bdecomp.reconstruct(bufwv.b[0]);
                }

                /*
                                if (!bdecomp.memoryAllocationFailed) {
                                    if ((lp.noisecf >= 0.1f ||  lp.noisecc >= 0.1f)) {
                                        //    fftw_denoise(bfw, bfh, max_numblox_W, min_numblox_W, bufwv.b, Bin,  numThreads, lp, 1);

                                    }
                                }
                */

                DeNoise_Local(call, lp, levred, hueplus, huemoins, huerefblur, dhueden, original, transformed, bufwv, cx, cy, sk);
            }

        }


//vibrance

        if (lp.expvib && (lp.past != 0.f  || lp.satur != 0.f)) { //interior ellipse renforced lightness and chroma  //locallutili
            //  printf("OK appel vib loc\n");
            float hueplus = hueref + dhuev;
            float huemoins = hueref - dhuev;

            //  printf ("hueplus=%f huemoins=%f dhu=%f\n", hueplus, huemoins, dhuev);

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhuev - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhuev + 2.f * rtengine::RT_PI;
            }

            LabImage *bufexporig = nullptr;
            LabImage *bufexpfin = nullptr;

            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;

            JaggedArray<float> buflight(bfw, bfh);
            JaggedArray<float> bufl_ab(bfw, bfh);


            if (call <= 3) { //simpleprocess, dcrop, improccoordinator


                bufexporig = new LabImage(bfw, bfh); //buffer for data in zone limit
                bufexpfin = new LabImage(bfw, bfh); //buffer for data in zone limit


#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufexporig->L[ir][jr] = 0.f;
                        bufexporig->a[ir][jr] = 0.f;
                        bufexporig->b[ir][jr] = 0.f;
                        bufexpfin->L[ir][jr] = 0.f;
                        bufexpfin->a[ir][jr] = 0.f;
                        bufexpfin->b[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufl_ab[ir][jr] = 0.f;


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

                            bufexporig->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufexporig->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufexporig->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas

                        }
                    }



                ImProcFunctions::vibrancelocal(sp, bfw, bfh, bufexporig, bufexpfin, localskutili, sklocalcurve);



#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {

                            float rL;
                            rL = CLIPRET((bufexpfin->L[loy - begy][lox - begx] - bufexporig->L[loy - begy][lox - begx]) / 328.f);

                            buflight[loy - begy][lox - begx] = rL;


                            float chp;
                            chp = CLIPRET((sqrt(SQR(bufexpfin->a[loy - begy][lox - begx]) + SQR(bufexpfin->b[loy - begy][lox - begx])) - sqrt(SQR(bufexporig->a[loy - begy][lox - begx]) + SQR(bufexporig->b[loy - begy][lox - begx]))) / 250.f);

                            bufl_ab[loy - begy][lox - begx] = chp;

                        }
                    }

                transit_shapedetect(2, bufexporig, nullptr, buflight, bufl_ab, nullptr, nullptr, nullptr, false, hueref, chromaref, lumaref, sobelref, nullptr, lp, original, transformed, cx, cy, sk);
                //call Expo_vibr_Local with first parameter = 2 for vibrance
                //              transit_shapedetect(6, loctemp, nullptr, buflight, bufchrom, nullptr, nullptr, nullptr, false, hueref, chromaref, lumaref, sobelref, nullptr, lp, original, transformed, cx, cy, sk);

            }

            if (call <= 3) {

                delete bufexporig;
                delete bufexpfin;


            }

        }


//Tone mapping

//&& lp.tonemapena
        if (lp.strengt != 0.f  && lp.tonemapena) {
            LabImage *tmp1 = nullptr;

            LabImage *bufgb = nullptr;
            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;
            JaggedArray<float> buflight(bfw, bfh);

            if (call <= 3) { //simpleprocess dcrop improcc

                bufgb = new LabImage(bfw, bfh);

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

                tmp1 = new LabImage(bfw, bfh);
                ImProcFunctions::EPDToneMaplocal(sp, bufgb, tmp1, 5, sk);
            } /*else { //stay here in case of

                tmp = new LabImage (transformed->W, transformed->H);
                tmp->CopyFrom (original);
                tmp1 = new LabImage (transformed->W, transformed->H);
                ImProcFunctions::EPDToneMaplocal (tmp, tmp1, 5 , sk);
                delete tmp;
            }
*/
            float hueplus = hueref + dhuetm;
            float huemoins = hueref - dhuetm;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhuetm - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhuetm + 2.f * rtengine::RT_PI;
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

                        float rL = CLIPRET((tmp1->L[loy - begy][lox - begx] - original->L[y][x]) / 400.f);

                        buflight[loy - begy][lox - begx]  = rL;

                    }
                }


            TM_Local(moddE, powdE, tmp1, buflight, hueplus, huemoins, hueref, dhuetm, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

            if (call <= 3) {
                delete bufgb;


            }

            delete tmp1;


        }

//begin cbdl
        if ((lp.mulloc[0] != 1.f || lp.mulloc[1] != 1.f || lp.mulloc[2] != 1.f || lp.mulloc[3] != 1.f || lp.mulloc[4] != 1.f) && lp.cbdlena) {
            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;
            JaggedArray<float> buflight(bfw, bfh);
            JaggedArray<float> bufchrom(bfw, bfh);
            JaggedArray<float> bufchr(bfw, bfh);
            JaggedArray<float> bufsh(bfw, bfh);
            LabImage *loctemp = nullptr;
            LabImage *loctempch = nullptr;

            float b_l = -5.f;
            float t_l = 25.f;
            float t_r = 120.f;
            float b_r = 170.f;
            double skinprot = 0.;
            int choice = 0;

            // I initialize these variable in case of !
            float hueplus = hueref + dhuecb;
            float huemoins = hueref - dhuecb;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhuecb - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhuecb + 2.f * rtengine::RT_PI;
            }


            if (call <= 3) { //call from simpleprocess dcrop improcc
                loctemp = new LabImage(bfw, bfh);
                loctempch = new LabImage(bfw, bfh);

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufsh[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufchr[ir][jr] = 0.f;
                        bufchrom[ir][jr] = 0.f;
                        loctemp->L[ir][jr] = 0.f;
                        loctemp->a[ir][jr] = 0.f;
                        loctemp->b[ir][jr] = 0.f;
                        loctempch->L[ir][jr] = 0.f;
                        loctempch->a[ir][jr] = 0.f;
                        loctempch->b[ir][jr] = 0.f;
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
                            bufchr[loy - begy][lox - begx] = sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x]));
                            loctemp->L[loy - begy][lox - begx] =  original->L[y][x];
                            loctemp->a[loy - begy][lox - begx] =  original->a[y][x];
                            loctemp->b[loy - begy][lox - begx] =  original->b[y][x];
                            loctempch->L[loy - begy][lox - begx] =  original->L[y][x];
                            loctempch->a[loy - begy][lox - begx] =  original->a[y][x];
                            loctempch->b[loy - begy][lox - begx] =  original->b[y][x];
                        }
                    }

                ImProcFunctions::cbdl_local_temp(bufsh, bufsh, loctemp->L, bfw, bfh, lp.mulloc, 1.f, lp.threshol, skinprot, false,  b_l, t_l, t_r, b_r, choice, sk);


#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            float rL;
                            rL = CLIPRET((loctemp->L[loy - begy][lox - begx] - original->L[y][x]) / 330.f);

                            buflight[loy - begy][lox - begx]  = rL;
                        }
                    }


                transit_shapedetect(6, loctemp, nullptr, buflight, bufchrom, nullptr, nullptr, nullptr, false, hueref, chromaref, lumaref, sobelref, nullptr, lp, original, transformed, cx, cy, sk);

                delete loctemp;

                //chroma CBDL begin here
                if (lp.chromacb > 0.f) {
                    if (lp.chromacb <= 1.f) {
                        lp.chromacb = 1.f;
                    }

                    float multc[5];


                    for (int lv = 0; lv < 5; lv++) {
                        multc[lv] = (lp.chromacb * ((float) lp.mulloc[lv] - 1.f) / 100.f) + 1.f;

                        if (multc[lv] <= 0.f) {
                            multc[lv] = 0.f;
                        }
                    }

                    {
                        ImProcFunctions::cbdl_local_temp(bufchr, bufchr, loctempch->L, bfw, bfh, multc, lp.chromacb, lp.threshol, skinprot, false,  b_l, t_l, t_r, b_r, choice, sk);

                        float rch;

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < transformed->H ; y++) //{
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    rch = CLIPRET((loctempch->L[loy - begy][lox - begx] - sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x]))) / 200.f);
                                    bufchrom[loy - begy][lox - begx]  = rch;
                                }
                            }
                    }

                    transit_shapedetect(7, loctempch, nullptr, buflight, bufchrom, nullptr, nullptr, nullptr, false, hueref, chromaref, lumaref, sobelref, nullptr, lp, original, transformed, cx, cy, sk);

                    delete loctempch;
                }
            }
        }


//end cbdl_Local

// soft light
        if (lp.strng > 0.f && call < 3  && lp.sfena) {
            float hueplus = hueref + dhuesf;
            float huemoins = hueref - dhuesf;


            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhuesf - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhuesf + 2.f * rtengine::RT_PI;
            }

            LabImage *bufexporig = nullptr;
            LabImage *bufexpfin = nullptr;

            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;

            JaggedArray<float> buflight(bfw, bfh);
            JaggedArray<float> bufl_ab(bfw, bfh);


            if (call <= 3) { //simpleprocess, dcrop, improccoordinator


                bufexporig = new LabImage(bfw, bfh); //buffer for data in zone limit
                bufexpfin = new LabImage(bfw, bfh); //buffer for data in zone limit


#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufexporig->L[ir][jr] = 0.f;
                        bufexporig->a[ir][jr] = 0.f;
                        bufexporig->b[ir][jr] = 0.f;
                        bufexpfin->L[ir][jr] = 0.f;
                        bufexpfin->a[ir][jr] = 0.f;
                        bufexpfin->b[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufl_ab[ir][jr] = 0.f;


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

                            bufexporig->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufexporig->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufexporig->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas

                        }
                    }


                ImProcFunctions::softLightloc(bufexporig, bufexpfin, lp.strng);


#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {

                            float rL;
                            rL = CLIPRET((bufexpfin->L[loy - begy][lox - begx] - bufexporig->L[loy - begy][lox - begx]) / 328.f);

                            buflight[loy - begy][lox - begx] = rL;


                            float chp;
                            chp = CLIPRET((sqrt(SQR(bufexpfin->a[loy - begy][lox - begx]) + SQR(bufexpfin->b[loy - begy][lox - begx])) - sqrt(SQR(bufexporig->a[loy - begy][lox - begx]) + SQR(bufexporig->b[loy - begy][lox - begx]))) / 250.f);

                            bufl_ab[loy - begy][lox - begx] = chp;

                        }
                    }

                transit_shapedetect(3, bufexporig, nullptr, buflight, bufl_ab, nullptr, nullptr, nullptr, false, hueref, chromaref, lumaref, sobelref, nullptr, lp, original, transformed, cx, cy, sk);

            }

            if (call <= 3) {

                delete bufexporig;
                delete bufexpfin;


            }


        }

//local contrast
        if (lp.lcamount > 0.f && call < 3  && lp.lcena) { //interior ellipse for sharpening, call = 1 and 2 only with Dcrop and simpleprocess
            int bfh = call == 2 ? int (lp.ly + lp.lyT) + del : original->H; //bfw bfh real size of square zone
            int bfw = call == 2 ? int (lp.lx + lp.lxL) + del : original->W;
            JaggedArray<float> loctemp(bfw, bfh);
            LabImage *bufloca = nullptr;

            if (call == 2) { //call from simpleprocess
                bufloca = new LabImage(bfw, bfh);

                //  JaggedArray<float> hbuffer(bfw, bfh);
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
                            bufloca->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                        }
                    }

                ImProcFunctions::localContrastloc(bufloca, sk, params->locallab.spots.at(sp).lcradius, params->locallab.spots.at(sp).lcamount, params->locallab.spots.at(sp).lcdarkness, params->locallab.spots.at(sp).lightness, loctemp);

            } else { //call from dcrop.cc

                ImProcFunctions::localContrastloc(original, sk, params->locallab.spots.at(sp).lcradius, params->locallab.spots.at(sp).lcamount, params->locallab.spots.at(sp).lcdarkness, params->locallab.spots.at(sp).lightness, loctemp);

            }

            float hueplus = hueref + dhuelc;
            float huemoins = hueref - dhuelc;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhuelc - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhuelc + 2.f * rtengine::RT_PI;
            }

            //sharpen ellipse and transition
            Sharp_Local(call, loctemp, 1, hueplus, huemoins, hueref, dhuesha, chromaref, lp, original, transformed, cx, cy, sk);

            delete bufloca;


        }


        if (!lp.invshar && lp.shrad > 0.42 && call < 3  && lp.sharpena) { //interior ellipse for sharpening, call = 1 and 2 only with Dcrop and simpleprocess
            int bfh = call == 2 ? int (lp.ly + lp.lyT) + del : original->H; //bfw bfh real size of square zone
            int bfw = call == 2 ? int (lp.lx + lp.lxL) + del : original->W;
            JaggedArray<float> loctemp(bfw, bfh);

            if (call == 2) { //call from simpleprocess
                JaggedArray<float> bufsh(bfw, bfh, true);
                JaggedArray<float> hbuffer(bfw, bfh);
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
                ImProcFunctions::deconvsharpeningloc(bufsh, hbuffer, bfw, bfh, loctemp, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius / 100., params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur / 100.);
            } else { //call from dcrop.cc

                ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, bfw, bfh, loctemp, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius / 100., params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur / 100.);

            }

            float hueplus = hueref + dhuesha;
            float huemoins = hueref - dhuesha;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhuesha - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhuesha + 2.f * rtengine::RT_PI;
            }

            //sharpen ellipse and transition
            Sharp_Local(call, loctemp, 0, hueplus, huemoins, hueref, dhuesha, chromaref, lp, original, transformed, cx, cy, sk);

        } else if (lp.invshar && lp.shrad > 0.42 && call < 3 && lp.sharpena) {
            int GW = original->W;
            int GH = original->H;
            JaggedArray<float> loctemp(GW, GH);

            ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, GW, GH, loctemp, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius / 100., params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur / 100.);

            float hueplus = hueref + dhuesha;
            float huemoins = hueref - dhuesha;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhuesha - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhuesha + 2.f * rtengine::RT_PI;
            }

            InverseSharp_Local(loctemp, hueplus, huemoins, hueref, dhuesha, chromaref, lp, original, transformed, cx, cy, sk);
        }

        //      }
//&& lp.retiena
        if (lp.str > 0.f  && lp.retiena) {
            int GW = transformed->W;
            int GH = transformed->H;

            LabImage *bufreti = nullptr;
            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;
            JaggedArray<float> buflight(bfw, bfh);
            JaggedArray<float> bufchro(bfw, bfh);

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
                bufreti = new LabImage(bfw, bfh);

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

                //calc dehaze
                Imagefloat *tmpImage = nullptr;

                if (lp.dehaze > 0) {
                    tmpImage = new Imagefloat(bfw, bfh);
                    lab2rgb(*bufreti, *tmpImage, params->icm.workingProfile);
                    float deha = LIM01(float(0.9f * lp.dehaze + 0.3f * lp.str) / 100.f * 0.9f);
                    float depthcombi = 0.3f * params->locallab.spots.at(sp).neigh + 0.15f * (500.f - params->locallab.spots.at(sp).vart);
                    float depth = -LIM01(depthcombi / 100.f);

                    dehazeloc(tmpImage, deha, depth);

                    rgb2lab(*tmpImage, *bufreti, params->icm.workingProfile);

                    delete tmpImage;
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

                tmpl = new LabImage(Wd, Hd);

            }  else {


#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int ir = 0; ir < Hd; ir += 1)
                    for (int jr = 0; jr < Wd; jr += 1) {
                        orig[ir][jr] = original->L[ir][jr];
                        orig1[ir][jr] = transformed->L[ir][jr];
                    }

                tmpl = new LabImage(transformed->W, transformed->H);


            }

            float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            ImProcFunctions::MSRLocal(sp, orig, tmpl->L, orig1, Wd, Hd, params->locallab, sk, locRETgainCcurve, 0, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int ir = 0; ir < Hd; ir += 1)
                for (int jr = 0; jr < Wd; jr += 1) {
                    tmpl->L[ir][jr] = orig[ir][jr];

                    if (!lp.invret) {
                        float rL;
                        rL = CLIPRET((tmpl->L[ir][jr] - bufreti->L[ir][jr]) / 328.f);
                        buflight[ir][jr] = rL;
                    }
                }


//new shape detection


            if (!lp.invret) {
                transit_shapedetect(4, bufreti, nullptr, buflight, bufchro, nullptr, nullptr, nullptr, false, hueref, chromaref, lumaref, sobelref, nullptr, lp, original, transformed, cx, cy, sk);

            } else {
                InverseReti_Local(lp, original, transformed, tmpl, cx, cy, 0);
            }

            if (params->locallab.spots.at(sp).chrrt > 0) {

                if (!lp.invret && call <= 3) {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int ir = 0; ir < Hd; ir += 1)
                        for (int jr = 0; jr < Wd; jr += 1) {

                            orig[ir][jr] = sqrt(SQR(bufreti->a[ir][jr]) + SQR(bufreti->b[ir][jr]));
                            orig1[ir][jr] = sqrt(SQR(bufreti->a[ir][jr]) + SQR(bufreti->b[ir][jr]));
                        }

                }  else {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int ir = 0; ir < GH; ir += 1)
                        for (int jr = 0; jr < GW; jr += 1) {
                            orig[ir][jr] = sqrt(SQR(original->a[ir][jr]) + SQR(original->b[ir][jr]));
                            orig1[ir][jr] = sqrt(SQR(transformed->a[ir][jr]) + SQR(transformed->b[ir][jr]));
                        }
                }

                ImProcFunctions::MSRLocal(sp, orig, tmpl->L, orig1, Wd, Hd, params->locallab, sk, locRETgainCcurve, 1, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);

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
                                ra = CLIPRET((sqrt(SQR(tmpl->a[ir][jr]) + SQR(tmpl->b[ir][jr])) - Chprov) / 300.f);
                                bufchro[ir][jr] = ra;
                            }

                        }



                }  else {

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


                if (!lp.invret) {
                    transit_shapedetect(5, tmpl, nullptr, buflight, bufchro, nullptr, nullptr, nullptr, false, hueref, chromaref, lumaref, sobelref, nullptr, lp, original, transformed, cx, cy, sk);

                } else {
                    InverseReti_Local(lp, original, transformed, tmpl, cx, cy, 1);
                }

            }

            delete tmpl;
            delete [] origBuffer;
            delete [] origBuffer1;

            if (!lp.invret && call <= 3) {

                delete  bufreti;


            }
        }


        if (lp.exposena && (lp.expcomp != 0.f || lp.war != 0 || lp.showmaskexpmet == 2 || lp.enaExpMask || lp.showmaskexpmet == 3 || lp.showmaskexpmet == 4 || (exlocalcurve  && localexutili))) { //interior ellipse renforced lightness and chroma  //locallutili
            float hueplus = hueref + dhueex;
            float huemoins = hueref - dhueex;


            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhueex - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhueex + 2.f * rtengine::RT_PI;
            }

            LabImage *bufexporig = nullptr;
            LabImage *bufexpfin = nullptr;
            LabImage *bufexptemp = nullptr;
            LabImage *bufcat02fin = nullptr;
            LabImage *bufmaskorigexp = nullptr;
            LabImage *bufmaskblurexp = nullptr;
            LabImage *originalmaskexp = nullptr;

            int bfh = 0.f, bfw = 0.f;
            bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            bfw = int (lp.lx + lp.lxL) + del;
            JaggedArray<float> buflight(bfw, bfh);
            JaggedArray<float> bufl_ab(bfw, bfh);
            JaggedArray<float> buflightcurv(bfw, bfh);
            JaggedArray<float> buf_a_cat(bfw, bfh, true);
            JaggedArray<float> buf_b_cat(bfw, bfh, true);
            JaggedArray<float> blend2(bfw, bfh);


            if (call <= 3) { //simpleprocess, dcrop, improccoordinator


                bufexporig = new LabImage(bfw, bfh);
                bufexpfin = new LabImage(bfw, bfh);
                bufexptemp = new LabImage(bfw, bfh);
                bufcat02fin = new LabImage(bfw, bfh);

                if (lp.showmaskexpmet == 2  || lp.enaExpMask || lp.showmaskexpmet == 3) {
                    int GWm = transformed->W;
                    int GHm = transformed->H;
                    bufmaskorigexp = new LabImage(bfw, bfh);
                    bufmaskblurexp = new LabImage(bfw, bfh);
                    originalmaskexp = new LabImage(GWm, GHm);
                }


#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufexporig->L[ir][jr] = 0.f;
                        bufexporig->a[ir][jr] = 0.f;
                        bufexporig->b[ir][jr] = 0.f;

                        if (lp.showmaskexpmet == 2  || lp.enaExpMask || lp.showmaskexpmet == 3) {
                            bufmaskorigexp->L[ir][jr] = 0.f;
                            bufmaskorigexp->a[ir][jr] = 0.f;
                            bufmaskorigexp->b[ir][jr] = 0.f;
                            bufmaskblurexp->L[ir][jr] = 0.f;
                            bufmaskblurexp->a[ir][jr] = 0.f;
                            bufmaskblurexp->b[ir][jr] = 0.f;
                        }

                        bufexptemp->L[ir][jr] = 0.f;
                        bufexptemp->a[ir][jr] = 0.f;
                        bufexptemp->b[ir][jr] = 0.f;
                        bufexpfin->L[ir][jr] = 0.f;
                        bufexpfin->a[ir][jr] = 0.f;
                        bufexpfin->b[ir][jr] = 0.f;
                        bufcat02fin->L[ir][jr] = 0.f;
                        bufcat02fin->a[ir][jr] = 0.f;
                        bufcat02fin->b[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufl_ab[ir][jr] = 0.f;
                        buflightcurv[ir][jr] = 0.f;
                        buf_a_cat[ir][jr] = 0.f;
                        buf_b_cat[ir][jr] = 0.f;
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
                            bufexporig->L[loy - begy][lox - begx] = original->L[y][x];
                        }
                    }

                const float radius = 3.f / (sk * 1.4f);
                int spotSi = 1 + 2 * max(1,  lp.cir / sk);

                if (spotSi < 5) {
                    spotSi = 5;
                }

                if (bfw > 2 * spotSi && bfh > 2 * spotSi  && lp.struexp > 0.f) {
                    SobelCannyLuma(blend2, bufexporig->L, bfw, bfh, radius);
                    array2D<float> ble(bfw, bfh);
                    array2D<float> guid(bfw, bfh);
#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int ir = 0; ir < bfh; ir++)
                        for (int jr = 0; jr < bfw; jr++) {
                            ble[ir][jr] = blend2[ir][jr] / 32768.f;
                            guid[ir][jr] = bufexporig->L[ir][jr] / 32768.f;
                        }

                    float blur = 25 / sk * (10.f + 1.2f * lp.struexp);

                    rtengine::guidedFilter(guid, ble, ble, blur, 0.001, multiThread);

#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int ir = 0; ir < bfh; ir++)
                        for (int jr = 0; jr < bfw; jr++) {
                            blend2[ir][jr] = ble[ir][jr] * 32768.f;
                        }

                    bool execmedian = true;
                    int passes = 1;

                    if (execmedian) {
                        float** tmL;
                        int wid = bfw;
                        int hei = bfh;
                        tmL = new float*[hei];

                        for (int i = 0; i < hei; ++i) {
                            tmL[i] = new float[wid];
                        }

                        Median medianTypeL = Median::TYPE_3X3_STRONG;
                        Median_Denoise(blend2, blend2, wid, hei, medianTypeL, passes, multiThread, tmL);

                        for (int i = 0; i < hei; ++i) {
                            delete[] tmL[i];
                        }

                        delete[] tmL;
                    }

                    if (lp.showmaskexpmet == 4) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < transformed->H ; y++)
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;
                                int zone = 0;
                                float localFactor = 1.f;
                                const float achm = (float)lp.trans / 100.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else if (lp.shapmet == 1) {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    if (zone > 0) {
                                        transformed->L[y][x] = blend2[loy - begy][lox - begx];
                                        transformed->a[y][x] = 0.f;
                                        transformed->b[y][x] = 0.f;
                                    }
                                }
                            }

                        return;
                    }


                }


#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            if (lp.showmaskexpmet == 2  || lp.enaExpMask || lp.showmaskexpmet == 3) {
                                bufmaskorigexp->L[loy - begy][lox - begx] = original->L[y][x];
                                bufmaskorigexp->a[loy - begy][lox - begx] = original->a[y][x];
                                bufmaskorigexp->b[loy - begy][lox - begx] = original->b[y][x];
                                bufmaskblurexp->L[loy - begy][lox - begx] = original->L[y][x];
                                bufmaskblurexp->a[loy - begy][lox - begx] = original->a[y][x];
                                bufmaskblurexp->b[loy - begy][lox - begx] = original->b[y][x];
                            }

                            bufexporig->L[loy - begy][lox - begx] = original->L[y][x];
                            bufexporig->a[loy - begy][lox - begx] = original->a[y][x];
                            bufexporig->b[loy - begy][lox - begx] = original->b[y][x];

                            float valLLexp = 0.f;
                            float valCC = 0.f;
                            float valHH = 0.f;
                            float kmaskLexp = 0;
                            float kmaskCa = 0;
                            float kmaskCb = 0;

                            float kmaskHL = 0;
                            float kmaskHa = 0;
                            float kmaskHb = 0;


                            if (lp.showmaskexpmet == 2  || lp.enaExpMask || lp.showmaskexpmet == 3) {

                                if (locllmasexpCurve  && llmasexputili) {
                                    float ligh = (bufexporig->L[loy - begy][lox - begx]) / 32768.f;
                                    valLLexp = (float)(locllmasexpCurve[500.f * ligh]);
                                    valLLexp = 1.f - valLLexp;
                                    kmaskLexp = 32768.f * valLLexp;
                                }

                                if (locccmasexpCurve && lcmasexputili) {
                                    float chromask = 0.0001f + (sqrt(SQR(bufexporig->a[loy - begy][lox - begx]) + SQR(bufexporig->b[loy - begy][lox - begx])));
                                    float chromaskr = chromask / 45000.f;
                                    valCC = float (locccmasexpCurve[500.f *  chromaskr]);
                                    valCC = 1.f - valCC;
                                    kmaskCa = valCC;
                                    kmaskCb = valCC;
                                }


                                if (lochhmasexpCurve && lhmasexputili) {
                                    float huema = xatan2f(bufexporig->b[loy - begy][lox - begx], bufexporig->a[loy - begy][lox - begx]);
                                    float h = Color::huelab_to_huehsv2(huema);
                                    h += 1.f / 6.f;

                                    if (h > 1.f) {
                                        h -= 1.f;
                                    }

                                    valHH = float (lochhmasexpCurve[500.f *  h]);
                                    valHH = 1.f - valHH;
                                    kmaskHa = valHH;
                                    kmaskHb = valHH;
                                    kmaskHL = 32768.f * valHH;
                                }

                                bufmaskblurexp->L[loy - begy][lox - begx] = CLIPLOC(kmaskLexp + kmaskHL);
                                bufmaskblurexp->a[loy - begy][lox - begx] = CLIPC(kmaskCa + kmaskHa);
                                bufmaskblurexp->b[loy - begy][lox - begx] = CLIPC(kmaskCb + kmaskHb);


                            }

                        }
                    }

                float radiusb = 3.f / sk;

                if (lp.showmaskexpmet == 2 || lp.enaExpMask || lp.showmaskexpmet == 3) {

#ifdef _OPENMP
                    #pragma omp parallel
#endif
                    {
                        gaussianBlur(bufmaskblurexp->L, bufmaskorigexp->L, bfw, bfh, radiusb);
                        gaussianBlur(bufmaskblurexp->a, bufmaskorigexp->a, bfw, bfh, radiusb);
                        gaussianBlur(bufmaskblurexp->b, bufmaskorigexp->b, bfw, bfh, radiusb);
                    }

                    delete bufmaskblurexp;


                    if (lp.showmaskexpmet != 3 || lp.enaExpMask) {

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < transformed->H ; y++) //{
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;
                                int zone = 0;

                                float localFactor = 1.f;
                                const float achm = (float)lp.trans / 100.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else if (lp.shapmet == 1) {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    if (zone > 0) {
                                        bufexporig->L[loy - begy][lox - begx] += (lp.blendmaexp * bufmaskorigexp->L[loy - begy][lox - begx]);
                                        bufexporig->a[loy - begy][lox - begx] *= (1.f + lp.blendmaexp * bufmaskorigexp->a[loy - begy][lox - begx]);
                                        bufexporig->b[loy - begy][lox - begx] *= (1.f + lp.blendmaexp * bufmaskorigexp->b[loy - begy][lox - begx]);

                                        bufexporig->L[loy - begy][lox - begx] = CLIP(bufexporig->L[loy - begy][lox - begx]);
                                        bufexporig->a[loy - begy][lox - begx] = CLIPC(bufexporig->a[loy - begy][lox - begx]);
                                        bufexporig->b[loy - begy][lox - begx] = CLIPC(bufexporig->b[loy - begy][lox - begx]);

                                        originalmaskexp->L[y][x] = CLIP(bufexporig->L[loy - begy][lox - begx] -  bufmaskorigexp->L[loy - begy][lox - begx]);
                                        originalmaskexp->a[y][x] = CLIPC(bufexporig->a[loy - begy][lox - begx] * (1.f - bufmaskorigexp->a[loy - begy][lox - begx]));
                                        originalmaskexp->b[y][x] = CLIPC(bufexporig->b[loy - begy][lox - begx] * (1.f - bufmaskorigexp->b[loy - begy][lox - begx]));

                                        switch (zone) {

                                            case 1: {
                                                original->L[y][x] += (lp.blendmaexp * localFactor * bufmaskorigexp->L[loy - begy][lox - begx]);
                                                original->a[y][x] *= (1.f + lp.blendmaexp * localFactor * bufmaskorigexp->a[loy - begy][lox - begx]);
                                                original->b[y][x] *= (1.f + lp.blendmaexp * localFactor * bufmaskorigexp->b[loy - begy][lox - begx]);
                                                original->L[y][x] = CLIP(original->L[y][x]);
                                                original->a[y][x] = CLIPC(original->a[y][x]);
                                                original->b[y][x] = CLIPC(original->b[y][x]);
                                                break;
                                            }

                                            case 2: {

                                                original->L[y][x] = bufexporig->L[loy - begy][lox - begx];
                                                original->a[y][x] = bufexporig->a[loy - begy][lox - begx];
                                                original->b[y][x] = bufexporig->b[loy - begy][lox - begx];

                                            }

                                        }
                                    }


                                }
                            }

                        delete bufmaskorigexp;

                    } else if (lp.showmaskexpmet == 3) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < transformed->H ; y++) //{
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;
                                int zone = 0;
                                float localFactor = 1.f;
                                const float achm = (float)lp.trans / 100.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else if (lp.shapmet == 1) {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    if (zone > 0) {
                                        transformed->L[y][x] = CLIPLOC(bufmaskorigexp->L[loy - begy][lox - begx]);
                                        transformed->a[y][x] = bufexporig->a[loy - begy][lox - begx] * (bufmaskorigexp->a[loy - begy][lox - begx]);
                                        transformed->b[y][x] = bufexporig->b[loy - begy][lox - begx] * (bufmaskorigexp->b[loy - begy][lox - begx]);
                                    }
                                }
                            }

                        delete bufmaskorigexp;
                        delete bufexporig;

                        return;
                    }
                }


                if (lp.showmaskexpmet == 0 || lp.showmaskexpmet == 1  || lp.showmaskexpmet == 2 || lp.enaExpMask) {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int y = 0; y < transformed->H ; y++) //{
                        for (int x = 0; x < transformed->W; x++) {
                            int lox = cx + x;
                            int loy = cy + y;

                            if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {

                                bufexptemp->L[loy - begy][lox - begx] = original->L[y][x];
                                bufexptemp->a[loy - begy][lox - begx] = original->a[y][x];
                                bufexptemp->b[loy - begy][lox - begx] = original->b[y][x];
                                bufexpfin->L[loy - begy][lox - begx] = original->L[y][x];
                                bufexpfin->a[loy - begy][lox - begx] = original->a[y][x];
                                bufexpfin->b[loy - begy][lox - begx] = original->b[y][x];
                            }
                        }

//to do Modulate bufexporig and bufexptemp with blend L, H, C and masks


                    if (exlocalcurve  && localexutili) {// L=f(L) curve enhanced
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < transformed->H ; y++) //{
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {

                                    float lighn = bufexporig->L[loy - begy][lox - begx];

                                    float lh;
                                    lh = 0.5f * exlocalcurve[2.f * lighn]; // / ((lighn) / 1.9f) / 3.61f; //lh between 0 and 0 50 or more
                                    bufexptemp->L[loy - begy][lox - begx] = lh;
                                }
                            }

                        if (lp.expcomp == 0.f) {
                            lp.expcomp = 0.1f;    // to enabled
                        }

                        ImProcFunctions::exlabLocal(lp, bfh, bfw, bufexptemp, bufexpfin, hltonecurveloc, shtonecurveloc, tonecurveloc);


                    } else {

                        ImProcFunctions::exlabLocal(lp, bfh, bfw, bufexporig, bufexpfin, hltonecurveloc, shtonecurveloc, tonecurveloc);
                    }

                    //cat02
                    if (params->locallab.spots.at(sp).warm != 0) {
                        ImProcFunctions::ciecamloc_02float(sp, bufexpfin, bufcat02fin);
                    } else {
#ifdef _OPENMP
                        #pragma omp parallel for
#endif

                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                bufcat02fin->L[ir][jr] = bufexpfin->L[ir][jr];
                                bufcat02fin->a[ir][jr] = bufexpfin->a[ir][jr];
                                bufcat02fin->b[ir][jr] = bufexpfin->b[ir][jr];
                            }
                    }



#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int y = 0; y < transformed->H ; y++) //{
                        for (int x = 0; x < transformed->W; x++) {
                            int lox = cx + x;
                            int loy = cy + y;

                            if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {

                                float rL;
                                rL = CLIPRET((bufcat02fin->L[loy - begy][lox - begx] - bufexporig->L[loy - begy][lox - begx]) / 328.f);

                                buflight[loy - begy][lox - begx] = rL;
                                float rA;
                                rA = CLIPRET((bufcat02fin->a[loy - begy][lox - begx] - bufexporig->a[loy - begy][lox - begx]) / 328.f);
                                buf_a_cat[loy - begy][lox - begx] = rA;


                                float rB;
                                rB = CLIPRET((bufcat02fin->b[loy - begy][lox - begx] - bufexporig->b[loy - begy][lox - begx]) / 328.f);
                                buf_b_cat[loy - begy][lox - begx] = rB;


                            }
                        }
                }

                transit_shapedetect(1, bufexporig, originalmaskexp, buflight, bufl_ab, buf_a_cat, buf_b_cat, nullptr, false, hueref, chromaref, lumaref, sobelref, blend2, lp, original, transformed, cx, cy, sk);
                //view mask
            }

            if (call <= 3) {

                delete bufexporig;
                delete bufexpfin;
                delete bufexptemp;
                delete bufcat02fin;

                if (lp.showmaskexpmet == 2 || lp.enaExpMask || lp.showmaskexpmet == 3) {
                    delete originalmaskexp;
                }
            }

        }


//local color and light

        if (!lp.inv  && (lp.chro != 0 || lp.ligh != 0.f || lp.cont != 0 || lp.qualcurvemet != 0 || lp.showmaskcolmet == 2 || lp.enaColorMask || lp.showmaskcolmet == 3  || lp.showmaskcolmet == 4) && lp.colorena) { // || lllocalcurve)) { //interior ellipse renforced lightness and chroma  //locallutili
//        if (!lp.inv  && (lp.chro != 0 || lp.ligh != 0.f || lp.cont != 0  || lp.showmaskcolmet == 2 || lp.enaColorMask || lp.showmaskcolmet == 3  || lp.showmaskcolmet == 4) && lp.colorena) { // || lllocalcurve)) { //interior ellipse renforced lightness and chroma  //locallutili
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

//provisory fixe for contrast slider does not work at all
            //  lp.cont = 0;

            LabImage *bufcolorig = nullptr;
            LabImage *bufmaskorigcol = nullptr;
            LabImage *bufmaskblurcol = nullptr;
            LabImage *originalmaskcol = nullptr;

            float chprosl = 1.f;
            float chprocu = 1.f;
            //   float cligh = 1.f;
            //   float clighL = 1.f;

            int bfh = 0.f, bfw = 0.f;
            bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            bfw = int (lp.lx + lp.lxL) + del;
            JaggedArray<float> buflight(bfw, bfh);
            JaggedArray<float> bufchro(bfw, bfh);
            JaggedArray<float> buflightslid(bfw, bfh);
            JaggedArray<float> bufchroslid(bfw, bfh);
            JaggedArray<float> bufhh(bfw, bfh);
            JaggedArray<float> blend2(bfw, bfh);
            JaggedArray<float> buforigchro(bfw, bfh);
            JaggedArray<float> buf_a(bfw, bfh);
            JaggedArray<float> buf_b(bfw, bfh);

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
                bufcolorig = new LabImage(bfw, bfh);

                if (lp.showmaskcolmet == 2  || lp.enaColorMask || lp.showmaskcolmet == 3) {
                    bufmaskorigcol = new LabImage(bfw, bfh);
                    bufmaskblurcol = new LabImage(bfw, bfh);
                    int GWm = transformed->W;
                    int GHm = transformed->H;
                    originalmaskcol = new LabImage(GWm, GHm);
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufcolorig->L[ir][jr] = 0.f;
                        bufcolorig->a[ir][jr] = 0.f;
                        bufcolorig->b[ir][jr] = 0.f;

                        if (lp.showmaskcolmet == 2  || lp.enaColorMask || lp.showmaskcolmet == 3) {
                            bufmaskorigcol->L[ir][jr] = 0.f;
                            bufmaskorigcol->a[ir][jr] = 0.f;
                            bufmaskorigcol->b[ir][jr] = 0.f;
                            bufmaskblurcol->L[ir][jr] = 0.f;
                            bufmaskblurcol->a[ir][jr] = 0.f;
                            bufmaskblurcol->b[ir][jr] = 0.f;
                        }

                        bufchro[ir][jr] = 0.f;
                        buf_a[ir][jr] = 0.f;
                        buf_b[ir][jr] = 0.f;
                        bufchroslid[ir][jr] = 0.f;
                        buflightslid[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufhh[ir][jr] = 0.f;
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
                            bufcolorig->L[loy - begy][lox - begx] = original->L[y][x];
                            //    buforigchro[loy - begy][lox - begx] = sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x]));
                        }
                    }

                const float radius = 3.f / (sk * 1.4f);
                int spotSi = 1 + 2 * max(1,  lp.cir / sk);

                if (spotSi < 5) {
                    spotSi = 5;
                }

                if (bfw > 2 * spotSi && bfh > 2 * spotSi  && lp.struco > 0.f) {
                    SobelCannyLuma(blend2, bufcolorig->L, bfw, bfh, radius);
                    array2D<float> ble(bfw, bfh);
                    array2D<float> blec(bfw, bfh);
                    array2D<float> guid(bfw, bfh);
#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int ir = 0; ir < bfh; ir++)
                        for (int jr = 0; jr < bfw; jr++) {
                            ble[ir][jr] = blend2[ir][jr] / 32768.f;
                            guid[ir][jr] = bufcolorig->L[ir][jr] / 32768.f;
                        }

                    float blur = 25 / sk * (10.f + 1.2f * lp.struco);
                    printf("Blur=%f \n", blur);

                    rtengine::guidedFilter(guid, ble, ble, blur, 0.001, multiThread);
                    //    rtengine::guidedFilter(guid, blec, blec, blur, 0.001, multiThread);

#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int ir = 0; ir < bfh; ir++)
                        for (int jr = 0; jr < bfw; jr++) {
                            blend2[ir][jr] = ble[ir][jr] * 32768.f;
                        }

                    bool execmedian = true;
                    int passes = 1;

                    if (execmedian) {
                        float** tmL;
                        int wid = bfw;
                        int hei = bfh;
                        tmL = new float*[hei];

                        for (int i = 0; i < hei; ++i) {
                            tmL[i] = new float[wid];
                        }

                        Median medianTypeL = Median::TYPE_3X3_STRONG;
                        Median_Denoise(blend2, blend2, wid, hei, medianTypeL, passes, multiThread, tmL);

                        for (int i = 0; i < hei; ++i) {
                            delete[] tmL[i];
                        }

                        delete[] tmL;
                    }

                    if (lp.showmaskcolmet == 4) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < transformed->H ; y++) //{
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;
                                int zone = 0;
                                float localFactor = 1.f;
                                const float achm = (float)lp.trans / 100.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else if (lp.shapmet == 1) {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    if (zone > 0) {
                                        transformed->L[y][x] = blend2[loy - begy][lox - begx];
                                        transformed->a[y][x] = 0.f;
                                        transformed->b[y][x] = 0.f;
                                    }
                                }
                            }

                        return;
                    }


                }

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            if (lp.showmaskcolmet == 2  || lp.enaColorMask || lp.showmaskcolmet == 3) {
                                bufmaskorigcol->L[loy - begy][lox - begx] = original->L[y][x];
                                bufmaskorigcol->a[loy - begy][lox - begx] = original->a[y][x];
                                bufmaskorigcol->b[loy - begy][lox - begx] = original->b[y][x];
                                bufmaskblurcol->L[loy - begy][lox - begx] = original->L[y][x];
                                bufmaskblurcol->a[loy - begy][lox - begx] = original->a[y][x];
                                bufmaskblurcol->b[loy - begy][lox - begx] = original->b[y][x];
                            }

                            bufcolorig->L[loy - begy][lox - begx] = original->L[y][x];
                            bufcolorig->a[loy - begy][lox - begx] = original->a[y][x];
                            bufcolorig->b[loy - begy][lox - begx] = original->b[y][x];

                            float valLL = 0.f;
                            float valCC = 0.f;
                            float valHH = 0.f;
                            float kmaskL = 0;
                            float kmaskCa = 0;
                            float kmaskCb = 0;

                            float kmaskHL = 0;
                            float kmaskHa = 0;
                            float kmaskHb = 0;

                            if (lp.showmaskcolmet == 2  || lp.enaColorMask || lp.showmaskcolmet == 3) {

                                if (locllmasCurve && llmasutili) {
                                    float ligh = (bufcolorig->L[loy - begy][lox - begx]) / 32768.f;
                                    valLL = (float)(locllmasCurve[500.f * ligh]);
                                    valLL = 1.f - valLL;
                                    kmaskL = 32768.f * valLL;
                                }

                                if (locccmasCurve && lcmasutili) {
                                    float chromask = 0.0001f + (sqrt(SQR(bufcolorig->a[loy - begy][lox - begx]) + SQR(bufcolorig->b[loy - begy][lox - begx])));
                                    float chromaskr = chromask / 45000.f;
                                    valCC = float (locccmasCurve[500.f *  chromaskr]);
                                    valCC = 1.f - valCC;
                                    kmaskCa = valCC;
                                    kmaskCb = valCC;
                                }

                                if (lochhmasCurve && lhmasutili) {
                                    float huema = xatan2f(bufcolorig->b[loy - begy][lox - begx], bufcolorig->a[loy - begy][lox - begx]);
                                    float h = Color::huelab_to_huehsv2(huema);
                                    h += 1.f / 6.f;

                                    if (h > 1.f) {
                                        h -= 1.f;
                                    }

                                    valHH = float (lochhmasCurve[500.f *  h]);
                                    valHH = 1.f - valHH;
                                    kmaskHa = valHH;
                                    kmaskHb = valHH;
                                    kmaskHL = 32768.f * valHH;
                                }

                                bufmaskblurcol->L[loy - begy][lox - begx] = CLIPLOC(kmaskL + kmaskHL);
                                bufmaskblurcol->a[loy - begy][lox - begx] = CLIPC(kmaskCa + kmaskHa);
                                bufmaskblurcol->b[loy - begy][lox - begx] = CLIPC(kmaskCb + kmaskHb);


                            }
                        }
                    }

                float radiusb = 3.f / sk;

                if (lp.showmaskcolmet == 2  || lp.enaColorMask || lp.showmaskcolmet == 3) {
#ifdef _OPENMP
                    #pragma omp parallel
#endif
                    {
                        gaussianBlur(bufmaskblurcol->L, bufmaskorigcol->L, bfw, bfh, radiusb);
                        gaussianBlur(bufmaskblurcol->a, bufmaskorigcol->a, bfw, bfh, radiusb);
                        gaussianBlur(bufmaskblurcol->b, bufmaskorigcol->b, bfw, bfh, radiusb);
                    }
                    delete bufmaskblurcol;

                    if (lp.showmaskcolmet != 3 || lp.enaColorMask) {

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < transformed->H ; y++) //{
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;
                                int zone = 0;

                                float localFactor = 1.f;
                                const float achm = (float)lp.trans / 100.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else if (lp.shapmet == 1) {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    if (zone > 0) {

                                        bufcolorig->L[loy - begy][lox - begx] += (lp.blendmacol * bufmaskorigcol->L[loy - begy][lox - begx]);
                                        bufcolorig->a[loy - begy][lox - begx] *= (1.f + lp.blendmacol * bufmaskorigcol->a[loy - begy][lox - begx]);
                                        bufcolorig->b[loy - begy][lox - begx] *= (1.f + lp.blendmacol * bufmaskorigcol->b[loy - begy][lox - begx]);

                                        bufcolorig->L[loy - begy][lox - begx] = CLIP(bufcolorig->L[loy - begy][lox - begx]);
                                        bufcolorig->a[loy - begy][lox - begx] = CLIPC(bufcolorig->a[loy - begy][lox - begx]);
                                        bufcolorig->b[loy - begy][lox - begx] = CLIPC(bufcolorig->b[loy - begy][lox - begx]);

                                        originalmaskcol->L[y][x] = CLIP(bufcolorig->L[loy - begy][lox - begx] -  bufmaskorigcol->L[loy - begy][lox - begx]);
                                        originalmaskcol->a[y][x] = CLIPC(bufcolorig->a[loy - begy][lox - begx] * (1.f - bufmaskorigcol->a[loy - begy][lox - begx]));
                                        originalmaskcol->b[y][x] = CLIPC(bufcolorig->b[loy - begy][lox - begx] * (1.f - bufmaskorigcol->b[loy - begy][lox - begx]));

                                        switch (zone) {

                                            case 1: {
                                                original->L[y][x] += (lp.blendmacol * localFactor * bufmaskorigcol->L[loy - begy][lox - begx]);
                                                original->a[y][x] *= (1.f + lp.blendmacol * localFactor * bufmaskorigcol->a[loy - begy][lox - begx]);
                                                original->b[y][x] *= (1.f + lp.blendmacol * localFactor * bufmaskorigcol->b[loy - begy][lox - begx]);
                                                original->L[y][x] = CLIP(original->L[y][x]);
                                                original->a[y][x] = CLIPC(original->a[y][x]);
                                                original->b[y][x] = CLIPC(original->b[y][x]);
                                                break;
                                            }

                                            case 2: {

                                                original->L[y][x] = bufcolorig->L[loy - begy][lox - begx];
                                                original->a[y][x] = bufcolorig->a[loy - begy][lox - begx];
                                                original->b[y][x] = bufcolorig->b[loy - begy][lox - begx];

                                            }

                                        }

                                    }
                                }
                            }

                        delete bufmaskorigcol;

                    } else if (lp.showmaskcolmet == 3) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < transformed->H ; y++) //{
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;
                                int zone = 0;
                                float localFactor = 1.f;
                                const float achm = (float)lp.trans / 100.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else if (lp.shapmet == 1) {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    if (zone > 0) {
                                        transformed->L[y][x] = CLIPLOC(bufmaskorigcol->L[loy - begy][lox - begx]);
                                        transformed->a[y][x] = bufcolorig->a[loy - begy][lox - begx] * (bufmaskorigcol->a[loy - begy][lox - begx]);
                                        transformed->b[y][x] = bufcolorig->b[loy - begy][lox - begx] * (bufmaskorigcol->b[loy - begy][lox - begx]);
                                    }
                                }
                            }

                        delete bufmaskorigcol;
                        delete bufcolorig;
                        return;

                    }
                }

                if (lp.showmaskcolmet == 0 || lp.showmaskcolmet == 1 || lp.showmaskcolmet == 2 || lp.enaColorMask) {

                    LabImage *bufcolcalc = nullptr;
                    bufcolcalc = new LabImage(bfw, bfh);

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int y = 0; y < transformed->H ; y++) //{
                        for (int x = 0; x < transformed->W; x++) {
                            int lox = cx + x;
                            int loy = cy + y;

                            if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                bufcolcalc->a[loy - begy][lox - begx] = bufcolorig->a[loy - begy][lox - begx];
                                bufcolcalc->b[loy - begy][lox - begx] = bufcolorig->b[loy - begy][lox - begx];
                                bufcolcalc->L[loy - begy][lox - begx] = bufcolorig->L[loy - begy][lox - begx];

                                if (cclocalcurve  && lp.qualcurvemet != 0) { // C=f(C) curve
                                    float chromat = sqrt(SQR(bufcolcalc->a[loy - begy][lox - begx]) +  SQR(bufcolcalc->b[loy - begy][lox - begx]));

                                    float ch;
                                    float ampli = 25.f;
                                    ch = (cclocalcurve[chromat * adjustr ])  / ((chromat + 0.00001f) * adjustr); //ch between 0 and 0 50 or more
                                    chprocu = CLIPCHRO(ampli * ch - ampli);  //ampli = 25.f arbitrary empirical coefficient between 5 and 50
                                }

                                if (lp.chro != 0.f) {
                                    float ch;
                                    float ampli = 70.f;
                                    ch = (1.f + 0.01f * lp.chro) ; //* (chromat * adjustr)  / ((chromat + 0.00001f) * adjustr); //ch between 0 and 0 50 or more

                                    if (ch <= 1.f) {//convert data curve near values of slider -100 + 100, to be used after to detection shape
                                        chprosl = 99.f * ch - 99.f;
                                    } else {
                                        chprosl = CLIPCHRO(ampli * ch - ampli);  //ampli = 25.f arbitrary empirical coefficient between 5 and 50
                                    }
                                }

                                bufchro[loy - begy][lox - begx] = chprosl + chprocu;

                                if (lochhCurve && HHutili && lp.qualcurvemet != 0) {
                                    float hhforcurv = xatan2f(bufcolcalc->b[loy - begy][lox - begx], bufcolcalc->a[loy - begy][lox - begx]);

                                    float valparam = float ((lochhCurve[500.f * Color::huelab_to_huehsv2(hhforcurv)] - 0.5f));  //get H=f(H)  1.7 optimisation !
                                    float ddhue = CLIPRET(200.f * valparam);
                                    bufhh[loy - begy][lox - begx] = ddhue;//valparamdh; //
                                }


                                if ((lp.ligh != 0.f || lp.cont != 0)) {
                                    float lighLnew = 0.f;
                                    float lighL = bufcolcalc->L[loy - begy][lox - begx];
                                    calclight(lighL, lp.ligh, lighLnew, lightCurveloc);  //replace L-curve
                                    bufcolcalc->L[loy - begy][lox - begx] = lighLnew;
                                }

                                if (lllocalcurve && locallutili  && lp.qualcurvemet != 0) {// L=f(L) curve enhanced
                                    float lh;
                                    float lighn = bufcolcalc->L[loy - begy][lox - begx];
                                    lh = 0.5f * (lllocalcurve[lighn * 2.f]);// / ((lighn + 0.00001f) * 1.9f) ; // / ((lighn) / 1.9f) / 3.61f; //lh between 0 and 0 50 or more
                                    bufcolcalc->L[loy - begy][lox - begx] = lh;
                                }

                                if (loclhCurve  && LHutili  && lp.qualcurvemet != 0) {
                                    float l_r;//Luminance Lab in 0..1
                                    float rhue = xatan2f(bufcolcalc->b[loy - begy][lox - begx], bufcolcalc->a[loy - begy][lox - begx]);
                                    float lighn = bufcolcalc->L[loy - begy][lox - begx];

                                    l_r =  lighn / 32768.f;
                                    {
                                        float khu = 1.9f; //in reserve in case of!

                                        float valparam = float ((loclhCurve[500.f * Color::huelab_to_huehsv2(rhue)] - 0.5f));  //get l_r=f(H)
                                        float valparamneg;
                                        valparamneg = valparam;

                                        if (valparam > 0.f) {
                                            l_r = (1.f - valparam) * l_r + valparam * (1.f - SQR(((SQR(1.f - min(l_r, 1.0f))))));
                                        } else
                                            //for negative
                                        {
                                            l_r *= (1.f + khu * valparamneg);
                                        }
                                    }

                                    bufcolcalc->L[loy - begy][lox - begx] = l_r * 32768.f;

                                }

                                float rL;
                                rL = CLIPRET((bufcolcalc->L[loy - begy][lox - begx] - bufcolorig->L[loy - begy][lox - begx]) / 328.f);
                                buflight[loy - begy][lox - begx] = rL;

                                float rA;
                                rA = CLIPRET((bufcolcalc->a[loy - begy][lox - begx] - bufcolorig->a[loy - begy][lox - begx]) / 328.f);
                                buf_a[loy - begy][lox - begx] = rA;


                                float rB;
                                rB = CLIPRET((bufcolcalc->b[loy - begy][lox - begx] - bufcolorig->b[loy - begy][lox - begx]) / 328.f);
                                buf_b[loy - begy][lox - begx] = rB;


                            }
                        }

                    delete bufcolcalc;
                }

                transit_shapedetect(0, bufcolorig, originalmaskcol, buflight, bufchro, buf_a, buf_b, bufhh, HHutili, hueref, chromaref, lumaref, sobelref, blend2, lp, original, transformed, cx, cy, sk);

                //     ColorLight_Local(moddE, powdE, call, bufcolorig, originalmaskcol, buflight, bufchro, bufchroslid, bufhh, buflightslid, LHutili, HHutili, hueplus, huemoins, hueref, dhue, chromaref, lumaref, sobelref, blend2, lllocalcurve, loclhCurve, lochhCurve, lightCurveloc, lp, original, transformed, cx, cy, sk);

                if (call <= 3) {

                    delete bufcolorig;

                    if (lp.showmaskcolmet == 2  || lp.enaColorMask || lp.showmaskcolmet == 3) {
                        delete originalmaskcol;
                    }

                }
            }
        }
//inverse
        else if (lp.inv  && (lp.chro != 0 || lp.ligh != 0) && lp.colorena) {
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            InverseColorLight_Local(lp, lightCurveloc, original, transformed, cx, cy,  hueplus, huemoins, hueref, dhue, chromaref, lumaref, sk);
        }

// Gamut and Munsell control - very important do not desactivated to avoid crash
        if (params->locallab.spots.at(sp).avoid) {
            TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);
            float wip[3][3] = {
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

                        float Chprov1 = sqrtf(SQR(bb) + SQR(aa));
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
                            HH = xatan2f(bb, aa);
                        }

                        float Chprov1 = sqrtf(SQR(aa) + SQR(bb)) / 327.68f;

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
                        Chprov1 = min(Chprov1, chr);

                        Color::gamutLchonly(sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.92f, neg, more_rgb);
#else
                        Color::pregamutlab(Lprov1, HH, chr);
                        Chprov1 = min(Chprov1, chr);
                        Color::gamutLchonly(sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.92f);
#endif

                        transformed->L[y][x] = Lprov1 * 327.68f;
                        transformed->a[y][x] = 327.68f * Chprov1 * sincosval.y;
                        transformed->b[y][x] = 327.68f * Chprov1 * sincosval.x;

                        if (needHH) {
                            float Lprov2 = original->L[y][x] / 327.68f;
                            float correctionHue = 0.f; // Munsell's correction
                            float correctlum = 0.f;
                            float memChprov = sqrtf(SQR(original->a[y][x]) + SQR(original->b[y][x])) / 327.68f;
                            float Chprov = sqrtf(SQR(transformed->a[y][x]) + SQR(transformed->b[y][x])) / 327.68f;
#ifdef _DEBUG
                            Color::AllMunsellLch(true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum, MunsDebugInfo);
#else
                            Color::AllMunsellLch(true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum);
#endif

                            if (fabs(correctionHue) < 0.015f) {
                                HH += correctlum;    // correct only if correct Munsell chroma very little.
                            }

                            float2 sincosval = xsincosf(HH + correctionHue);

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
            printf("Color::AllMunsellLch (correction performed in %d usec):\n", t2e.etime(t1e));
            //  printf("   Munsell chrominance: MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhue[0],    MunsDebugInfo->maxdhue[1],    MunsDebugInfo->maxdhue[2],    MunsDebugInfo->maxdhue[3],    MunsDebugInfo->depass);
            //  printf("   Munsell luminance  : MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhuelum[0], MunsDebugInfo->maxdhuelum[1], MunsDebugInfo->maxdhuelum[2], MunsDebugInfo->maxdhuelum[3], MunsDebugInfo->depassLum);
        }

        delete MunsDebugInfo;
#endif

    }

}

}
