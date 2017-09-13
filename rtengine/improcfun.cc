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
#include <cmath>
#include <glib.h>
#include <glibmm.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "rtengine.h"
#include "improcfun.h"
#include "curves.h"
#include "mytime.h"
#include "iccstore.h"
#include "impulse_denoise.h"
#include "imagesource.h"
#include "rtthumbnail.h"
#include "utils.h"
#include "iccmatrices.h"
#include "color.h"
#include "calc_distort.h"
#include "rt_math.h"
#include "EdgePreservingDecomposition.h"
#include "improccoordinator.h"
#include "clutstore.h"
#include "ciecam02.h"
//#define BENCHMARK
#include "StopWatch.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/guiutils.h"

#undef CLIPD
#define CLIPD(a) ((a)>0.0f?((a)<1.0f?(a):1.0f):0.0f)

namespace rtengine
{

using namespace procparams;

extern const Settings* settings;

ImProcFunctions::~ImProcFunctions ()
{

    if (monitorTransform) {
        cmsDeleteTransform (monitorTransform);
    }
}

void ImProcFunctions::setScale (double iscale)
{
    scale = iscale;
}

void ImProcFunctions::updateColorProfiles (const Glib::ustring& monitorProfile, RenderingIntent monitorIntent, bool softProof, bool gamutCheck)
{
    // set up monitor transform
    if (monitorTransform) {
        cmsDeleteTransform (monitorTransform);
    }

    monitorTransform = nullptr;

    cmsHPROFILE monitor = nullptr;

    if (!monitorProfile.empty()) {
#if !defined(__APPLE__) // No support for monitor profiles on OS X, all data is sRGB
        monitor = ICCStore::getInstance()->getProfile (monitorProfile);
#else
        monitor = ICCStore::getInstance()->getProfile ("RT_sRGB");
#endif
    }

    if (monitor) {
        MyMutex::MyLock lcmsLock (*lcmsMutex);

        cmsUInt32Number flags;
        cmsHPROFILE iprof  = cmsCreateLab4Profile (nullptr);

        bool softProofCreated = false;

        if (softProof) {
            cmsHPROFILE oprof = nullptr;

            if (!settings->printerProfile.empty()) {
                oprof = ICCStore::getInstance()->getProfile (settings->printerProfile);
            }

            if (oprof) {
                // NOCACHE is for thread safety, NOOPTIMIZE for precision
                flags = cmsFLAGS_SOFTPROOFING | cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

                if (settings->printerBPC) {
                    flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
                }

                if (gamutCheck) {
                    flags |= cmsFLAGS_GAMUTCHECK;
                }

                monitorTransform = cmsCreateProofingTransform (
                                       iprof, TYPE_Lab_FLT,
                                       monitor, TYPE_RGB_8,
                                       oprof,
                                       monitorIntent, settings->printerIntent,
                                       flags
                                   );

                if (monitorTransform) {
                    softProofCreated = true;
                }
            }
        }

        if (!softProofCreated) {
            flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

            if (settings->monitorBPC) {
                flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
            }

            monitorTransform = cmsCreateTransform (iprof, TYPE_Lab_FLT, monitor, TYPE_RGB_8, monitorIntent, flags);
        }

        cmsCloseProfile (iprof);
    }
}

void ImProcFunctions::firstAnalysis (const Imagefloat* const original, const ProcParams &params, LUTu & histogram)
{

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params.icm.working);

    lumimul[0] = wprof[1][0];
    lumimul[1] = wprof[1][1];
    lumimul[2] = wprof[1][2];
    int W = original->getWidth();
    int H = original->getHeight();

    float lumimulf[3] = {static_cast<float> (lumimul[0]), static_cast<float> (lumimul[1]), static_cast<float> (lumimul[2])};

    // calculate histogram of the y channel needed for contrast curve calculation in exposure adjustments
    histogram.clear();

    if (multiThread) {

#ifdef _OPENMP
        const int numThreads = min (max (W * H / (int)histogram.getSize(), 1), omp_get_max_threads());
        #pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
        {
            LUTu hist (histogram.getSize());
            hist.clear();
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int i = 0; i < H; i++) {
                for (int j = 0; j < W; j++) {

                    float r = original->r (i, j);
                    float g = original->g (i, j);
                    float b = original->b (i, j);

                    int y = (lumimulf[0] * r + lumimulf[1] * g + lumimulf[2] * b);
                    hist[y]++;
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            histogram += hist;

        }
#ifdef _OPENMP
        static_cast<void> (numThreads); // to silence cppcheck warning
#endif
    } else {
        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {

                float r = original->r (i, j);
                float g = original->g (i, j);
                float b = original->b (i, j);

                int y = (lumimulf[0] * r + lumimulf[1] * g + lumimulf[2] * b);
                histogram[y]++;
            }
        }
    }
}

// Copyright (c) 2012 Jacques Desmis <jdesmis@gmail.com>
void ImProcFunctions::ciecam_02 (CieImage* ncie, double adap, int begh, int endh, int pW, int pwb, LabImage* lab, const ProcParams* params,
                                 const ColorAppearance & customColCurve1, const ColorAppearance & customColCurve2, const ColorAppearance & customColCurve3,
                                 LUTu & histLCAM, LUTu & histCCAM, LUTf & CAMBrightCurveJ, LUTf & CAMBrightCurveQ, float &mean, int Iterates, int scale, bool execsharp, double &d, double &dj, double &yb, int rtt)
{
    if (params->colorappearance.enabled) {
//int lastskip;
//if(rtt==1) {lastskip=scale;} //not for Rtthumbnail

#ifdef _DEBUG
        MyTime t1e, t2e;
        t1e.set();
#endif
        LUTf dLcurve;
        LUTu hist16JCAM;
        bool jp = false;

        //preparate for histograms CIECAM
        if (pW != 1) { //only with improccoordinator
            dLcurve (65536, 0);
            dLcurve.clear();
            hist16JCAM (65536, 0);
            hist16JCAM.clear();

            for (int i = 0; i < 32768; i++) { //# 32768*1.414  approximation maxi for chroma
                float val = (double)i / 32767.0;
                dLcurve[i] = CLIPD (val);
            }
        }

        LUTf dCcurve;
        LUTu hist16_CCAM;
        bool chropC = false;

        if (pW != 1) { //only with improccoordinator
            dCcurve (65536, 0);
            hist16_CCAM (65536);
            hist16_CCAM.clear();

            for (int i = 0; i < 48000; i++) { //# 32768*1.414  approximation maxi for chroma
                float valc = (double)i / 47999.0;
                dCcurve[i] = CLIPD (valc);
            }
        }

        //end preparate histogram
        int width = lab->W, height = lab->H;
        float minQ = 10000.f;
        float maxQ = -1000.f;
        float w_h;
        float a_w;
        float c_;
        float f_l;
        double Yw;
        Yw = 1.0;
        double Xw, Zw;
        double Xwout, Zwout;
        double Xwsc, Zwsc;

        double f = 0., c = 0., nc = 0., yb = 0., la, xw, yw, zw, f2 = 0., c2 = 0., nc2 = 0., yb2 = 0., la2;
        double fl, n, nbb, ncb, aw;
        double xwd = 0., ywd, zwd = 0.;
        double xws, yws, zws;
        int alg = 0;
        bool algepd = false;
        float sum = 0.f;

        bool ciedata = params->colorappearance.datacie;

        ColorTemp::temp2mulxyz (params->wb.temperature, params->wb.green, params->wb.method, Xw, Zw); //compute white Xw Yw Zw  : white current WB
        ColorTemp::temp2mulxyz (params->colorappearance.tempout, params->colorappearance.greenout, "Custom", Xwout, Zwout);
        ColorTemp::temp2mulxyz (params->colorappearance.tempsc, params->colorappearance.greensc, "Custom", Xwsc, Zwsc);

        //viewing condition for surrsrc
        if (params->colorappearance.surrsrc == "Average") {
            f  = 1.00;
            c  = 0.69;
            nc = 1.00;
        } else if (params->colorappearance.surrsrc == "Dim") {
            f  = 0.9;
            c  = 0.59;
            nc = 0.9;
        } else if (params->colorappearance.surrsrc == "Dark") {
            f  = 0.8;
            c  = 0.525;
            nc = 0.8;
        } else if (params->colorappearance.surrsrc == "ExtremelyDark") {
            f  = 0.8;
            c  = 0.41;
            nc = 0.8;
        }


        //viewing condition for surround
        if (params->colorappearance.surround == "Average") {
            f2 = 1.0, c2 = 0.69, nc2 = 1.0;
        } else if (params->colorappearance.surround == "Dim") {
            f2  = 0.9;
            c2  = 0.59;
            nc2 = 0.9;
        } else if (params->colorappearance.surround == "Dark") {
            f2  = 0.8;
            c2  = 0.525;
            nc2 = 0.8;
        } else if (params->colorappearance.surround == "ExtremelyDark") {
            f2  = 0.8;
            c2  = 0.41;
            nc2 = 0.8;
        }

        /*
                //scene condition for surround
                if (params->colorappearance.surrsource)  {
                    f  = 0.85;    // if user => source image has surround very dark
                    c  = 0.55;
                    nc = 0.85;
                }
        */
        //with which algorithme
        if     (params->colorappearance.algo == "JC") {
            alg = 0;
        } else if (params->colorappearance.algo == "JS") {
            alg = 1;
        } else if (params->colorappearance.algo == "QM")  {
            alg = 2;
            algepd = true;
        } else if (params->colorappearance.algo == "ALL")  {
            alg = 3;
            algepd = true;
        }

        bool needJ = (alg == 0 || alg == 1 || alg == 3);
        bool needQ = (alg == 2 || alg == 3);
        /*
                //settings white point of output device - or illuminant viewing
                if (settings->viewingdevice == 0) {
                    xwd = 96.42;    //5000K
                    ywd = 100.0;
                    zwd = 82.52;
                } else if (settings->viewingdevice == 1) {
                    xwd = 95.68;    //5500
                    ywd = 100.0;
                    zwd = 92.15;
                } else if (settings->viewingdevice == 2) {
                    xwd = 95.24;    //6000
                    ywd = 100.0;
                    zwd = 100.81;
                } else if (settings->viewingdevice == 3)  {
                    xwd = 95.04;    //6500
                    ywd = 100.0;
                    zwd = 108.88;
                } else if (settings->viewingdevice == 4)  {
                    xwd = 109.85;    //tungsten
                    ywd = 100.0;
                    zwd = 35.58;
                } else if (settings->viewingdevice == 5)  {
                    xwd = 99.18;    //fluo F2
                    ywd = 100.0;
                    zwd = 67.39;
                } else if (settings->viewingdevice == 6)  {
                    xwd = 95.04;    //fluo F7
                    ywd = 100.0;
                    zwd = 108.75;
                } else if (settings->viewingdevice == 7)  {
                    xwd = 100.96;    //fluo F11
                    ywd = 100.0;
                    zwd = 64.35;
                }
        */

        xwd = 100. * Xwout;
        zwd = 100. * Zwout;
        ywd = 100. / params->colorappearance.greenout;//approximation to simplify

        xws = 100. * Xwsc;
        zws = 100. * Zwsc;
        yws = 100. / params->colorappearance.greensc;//approximation to simplify

        /*
                //settings mean Luminance Y of output device or viewing
                if (settings->viewingdevicegrey == 0) {
                    yb2 = 5.0;
                } else if (settings->viewingdevicegrey == 1) {
                    yb2 = 10.0;
                } else if (settings->viewingdevicegrey == 2) {
                    yb2 = 15.0;
                } else if (settings->viewingdevicegrey == 3) {
                    yb2 = 18.0;
                } else if (settings->viewingdevicegrey == 4) {
                    yb2 = 23.0;
                } else if (settings->viewingdevicegrey == 5)  {
                    yb2 = 30.0;
                } else if (settings->viewingdevicegrey == 6)  {
                    yb2 = 40.0;
                }
        */
        yb2 = params->colorappearance.ybout;

        //La and la2 = ambiant luminosity scene and viewing
        la = double (params->colorappearance.adapscen);

        if (pwb == 2) {
            if (params->colorappearance.autoadapscen) {
                la = adap;
            }
        }

        la2 = double (params->colorappearance.adaplum);

        // level of adaptation
        double deg = (params->colorappearance.degree) / 100.0;
        double pilot = params->colorappearance.autodegree ? 2.0 : deg;


        const float degout = (params->colorappearance.degreeout) / 100.0;
        const float pilotout = params->colorappearance.autodegreeout ? 2.0 : degout;

        //algoritm's params
        float jli = params->colorappearance.jlight;
        float chr = params->colorappearance.chroma;
        float contra = params->colorappearance.contrast;
        float qbri = params->colorappearance.qbright;
        float schr = params->colorappearance.schroma;
        float mchr = params->colorappearance.mchroma;
        float qcontra = params->colorappearance.qcontrast;
        float hue = params->colorappearance.colorh;
        double rstprotection = 100. - params->colorappearance.rstprotection;

        if (schr > 0.0) {
            schr = schr / 2.0f;    //divide sensibility for saturation
        }

        // extracting datas from 'params' to avoid cache flush (to be confirmed)
        ColorAppearanceParams::eTCModeId curveMode = params->colorappearance.curveMode;
        ColorAppearanceParams::eTCModeId curveMode2 = params->colorappearance.curveMode2;
        bool hasColCurve1 = bool (customColCurve1);
        bool hasColCurve2 = bool (customColCurve2);
        ColorAppearanceParams::eCTCModeId curveMode3 = params->colorappearance.curveMode3;
        bool hasColCurve3 = bool (customColCurve3);


        if (CAMBrightCurveJ.dirty || CAMBrightCurveQ.dirty) {
            LUTu hist16J;
            LUTu hist16Q;

            if (needJ) {
                hist16J (65536);
                hist16J.clear();
            }

            if (needQ) {
                hist16Q (65536);
                hist16Q.clear();
            }

            float koef = 1.0f; //rough correspondence between L and J

            for (int i = 0; i < height; i++)

                //   for (int i=begh; i<endh; i++)
                for (int j = 0; j < width; j++) { //rough correspondence between L and J
                    float currL = lab->L[i][j] / 327.68f;

                    if     (currL > 95.) {
                        koef = 1.f;
                    } else if (currL > 85.) {
                        koef = 0.97f;
                    } else if (currL > 80.) {
                        koef = 0.93f;
                    } else if (currL > 70.) {
                        koef = 0.87f;
                    } else if (currL > 60.) {
                        koef = 0.85f;
                    } else if (currL > 50.) {
                        koef = 0.8f;
                    } else if (currL > 40.) {
                        koef = 0.75f;
                    } else if (currL > 30.) {
                        koef = 0.7f;
                    } else if (currL > 20.) {
                        koef = 0.7f;
                    } else if (currL > 10.) {
                        koef = 0.9f;
                    } else if (currL > 0.) {
                        koef = 1.0f;
                    }

                    if (needJ) {
                        hist16J[CLIP ((int) ((koef * lab->L[i][j])))]++;  //evaluate histogram luminance L # J
                    }

                    if (needQ) {
                        hist16Q[CLIP ((int) (32768.f * sqrt ((koef * (lab->L[i][j])) / 32768.f)))]++;  //for brightness Q : approximation for Q=wh*sqrt(J/100)  J not equal L
                    }

                    sum += koef * lab->L[i][j]; //evaluate mean J to calcualte Yb
                }

            //mean=(sum/((endh-begh)*width))/327.68f;//for Yb  for all image...if one day "pipette" we can adapt Yb for each zone
            mean = (sum / ((height) * width)) / 327.68f; //for Yb  for all image...if one day "pipette" we can adapt Yb for each zone

            //evaluate lightness, contrast
            if (needJ) {
                if (!CAMBrightCurveJ) {
                    CAMBrightCurveJ (65536, 0);
                    CAMBrightCurveJ.dirty = false;
                }

                Ciecam02::curveJ (jli, contra, 1, CAMBrightCurveJ, hist16J);//lightness and contrast J
            }

            if (needQ) {
                if (!CAMBrightCurveQ) {
                    CAMBrightCurveQ (65536, 0);
                    CAMBrightCurveQ.dirty = false;
                }

                Ciecam02::curveJ (qbri, qcontra, 1, CAMBrightCurveQ, hist16Q);//brightness and contrast Q
            }
        }

        //   if (settings->viewinggreySc == 0) { //auto
        if (params->colorappearance.autoybscen  &&  pwb == 2) {//auto

            if     (mean < 15.f) {
                yb = 3.0;
            } else if (mean < 30.f) {
                yb = 5.0;
            } else if (mean < 40.f) {
                yb = 10.0;
            } else if (mean < 45.f) {
                yb = 15.0;
            } else if (mean < 50.f) {
                yb = 18.0;
            } else if (mean < 55.f) {
                yb = 23.0;
            } else if (mean < 60.f) {
                yb = 30.0;
            } else if (mean < 70.f) {
                yb = 40.0;
            } else if (mean < 80.f) {
                yb = 60.0;
            } else if (mean < 90.f) {
                yb = 80.0;
            } else {
                yb = 90.0;
            }
        } else {
            yb = params->colorappearance.ybscen;
        }

        if (settings->viewinggreySc == 1) {
            yb = 18.0;
        }

        int gamu = 0;
        bool highlight = params->toneCurve.hrenabled; //Get the value if "highlight reconstruction" is activated

        if (params->colorappearance.gamut) {
            gamu = 1;    //enabled gamut control
        }

        xw = 100.0 * Xw;
        yw = 100.0 * Yw;
        zw = 100.0 * Zw;
        double xw1 = 0., yw1 = 0., zw1 = 0., xw2 = 0., yw2 = 0., zw2 = 0.;

        // settings of WB: scene and viewing
        if (params->colorappearance.wbmodel == "RawT") {
            xw1 = 96.46;    //use RT WB; CAT 02 is used for output device (see prefreneces)
            yw1 = 100.0;
            zw1 = 82.445;
            xw2 = xwd;
            yw2 = ywd;
            zw2 = zwd;
        } else if (params->colorappearance.wbmodel == "RawTCAT02") {
            xw1 = xw;    // Settings RT WB are used for CAT02 => mix , CAT02 is use for output device (screen: D50 D65, projector: lamp, LED) see preferences
            yw1 = yw;
            zw1 = zw;
            xw2 = xwd;
            yw2 = ywd;
            zw2 = zwd;
        } else if (params->colorappearance.wbmodel == "free") {
            xw1 = xws;    // free temp and green
            yw1 = yws;
            zw1 = zws;
            xw2 = xwd;
            yw2 = ywd;
            zw2 = zwd;
        }

        double cz, wh, pfl;
        Ciecam02::initcam1 (gamu, yb, pilot, f, la, xw, yw, zw, n, d, nbb, ncb, cz, aw, wh, pfl, fl, c);
        double nj, nbbj, ncbj, czj, awj, flj;
        Ciecam02::initcam2 (gamu, yb2, pilotout, f2,  la2,  xw2,  yw2,  zw2, nj, dj, nbbj, ncbj, czj, awj, flj);




#ifndef _DEBUG
        #pragma omp parallel default(shared) firstprivate(lab,xw1,xw2,yw1,yw2,zw1,zw2,pilot,jli,chr,yb,la,yb2,la2,fl,nc,f,c, height,width,begh, endh,nc2,f2,c2, alg,algepd, gamu, highlight, rstprotection, pW, scale)
#endif
        {
            //matrix for current working space
            TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.working);
            double wip[3][3] = {
                {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
                {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
                {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
            };

#ifndef _DEBUG
            #pragma omp for schedule(dynamic, 10)
#endif

            for (int i = 0; i < height; i++)

//  for (int i=begh; i<endh; i++)
                for (int j = 0; j < width; j++) {

                    float L = lab->L[i][j];
                    float a = lab->a[i][j];
                    float b = lab->b[i][j];
                    float x1, y1, z1;
                    double x, y, z;
                    double epsil = 0.0001;
                    //convert Lab => XYZ
                    Color::Lab2XYZ (L, a, b, x1, y1, z1);
                    //  double J, C, h, Q, M, s, aw, fl, wh;
                    double J, C, h, Q, M, s;

                    double Jpro, Cpro, hpro, Qpro, Mpro, spro;
                    bool t1L = false;
                    bool t1B = false;
                    bool t2B = false;
                    int c1s = 0;
                    int c1co = 0;
                    //double n,nbb,ncb,pfl,cz,d;
                    x = (double)x1 / 655.35;
                    y = (double)y1 / 655.35;
                    z = (double)z1 / 655.35;
                    //process source==> normal
                    Ciecam02::xyz2jchqms_ciecam02 ( J, C,  h,
                                                    Q,  M,  s, aw, fl, wh,
                                                    x,  y,  z,
                                                    xw1, yw1,  zw1,
                                                    yb,  la,
                                                    f, c,  nc,  pilot, gamu, n, nbb, ncb, pfl, cz, d );
                    Jpro = J;
                    Cpro = C;
                    hpro = h;
                    Qpro = Q;
                    Mpro = M;
                    spro = s;
                    w_h = wh + epsil;
                    a_w = aw;
                    c_ = c;
                    f_l = fl;

                    // we cannot have all algoritms with all chroma curves
                    if (alg == 1)  {
                        // Lightness saturation
                        if (Jpro > 99.9f) {
                            Jpro = 99.9f;
                        }

                        Jpro = (CAMBrightCurveJ[ (float) (Jpro * 327.68)]) / 327.68; //ligthness CIECAM02 + contrast
                        double sres;
                        double Sp = spro / 100.0;
                        double parsat = 1.5; //parsat=1.5 =>saturation  ; 1.8 => chroma ; 2.5 => colorfullness (personal evaluation)

                        if (schr == -100.0) {
                            schr = -99.8;
                        }

                        Ciecam02::curvecolor (schr, Sp, sres, parsat);
                        double coe = pow (fl, 0.25);
                        float dred = 100.f; // in C mode
                        float protect_red = 80.0f; // in C mode
                        dred = 100.0 * sqrt ((dred * coe) / Qpro);
                        protect_red = 100.0 * sqrt ((protect_red * coe) / Qpro);
                        int sk = 0;
                        float ko = 100.f;
                        Color::skinred (Jpro, hpro, sres, Sp, dred, protect_red, sk, rstprotection, ko, spro);
                        Qpro = ( 4.0 / c ) * sqrt ( Jpro / 100.0 ) * ( aw + 4.0 ) ;
                        Cpro = (spro * spro * Qpro) / (10000.0);
                    } else if (alg == 3 || alg == 0  || alg == 2) {
                        double coef = 32760. / wh;

                        if (alg == 3 || alg == 2) {
                            if (Qpro * coef > 32767.0f) {
                                Qpro = (CAMBrightCurveQ[ (float)32767.0f]) / coef;   //brightness and contrast
                            } else {
                                Qpro = (CAMBrightCurveQ[ (float) (Qpro * coef)]) / coef;  //brightness and contrast
                            }
                        }

                        double Mp, sres;
                        double coe = pow (fl, 0.25);
                        Mp = Mpro / 100.0;
                        double parsat = 2.5;

                        if (mchr == -100.0) {
                            mchr = -99.8 ;
                        }

                        if (mchr == 100.0) {
                            mchr = 99.9;
                        }

                        if (alg == 3 || alg == 2) {
                            Ciecam02::curvecolor (mchr, Mp, sres, parsat);
                        } else {
                            Ciecam02::curvecolor (0.0, Mp, sres, parsat);    //colorfullness
                        }

                        float dred = 100.f; //in C mode
                        float protect_red = 80.0f; // in C mode
                        dred *= coe; //in M mode
                        protect_red *= coe; //M mode
                        int sk = 0;
                        float ko = 100.f;
                        Color::skinred (Jpro, hpro, sres, Mp, dred, protect_red, sk, rstprotection, ko, Mpro);
                        Jpro = (100.0 * Qpro * Qpro) / (wh * wh);
                        Cpro = Mpro / coe;
                        spro = 100.0 * sqrt ( Mpro / Qpro );

                        if (alg != 2) {
                            if (Jpro > 99.9f) {
                                Jpro = 99.9f;
                            }

                            Jpro = (CAMBrightCurveJ[ (float) (Jpro * 327.68f)]) / 327.68f; //ligthness CIECAM02 + contrast
                        }

                        double Cp;
                        double Sp = spro / 100.0;
                        parsat = 1.5;

                        if (schr == -100.0) {
                            schr = -99.;
                        }

                        if (schr == 100.0) {
                            schr = 98.;
                        }

                        if (alg == 3) {
                            Ciecam02::curvecolor (schr, Sp, sres, parsat);
                        }   else {
                            Ciecam02::curvecolor (0.0, Sp, sres, parsat);    //saturation
                        }

                        dred = 100.f; // in C mode
                        protect_red = 80.0f; // in C mode
                        dred = 100.0 * sqrt ((dred * coe) / Q);
                        protect_red = 100.0 * sqrt ((protect_red * coe) / Q);
                        sk = 0;
                        Color::skinred (Jpro, hpro, sres, Sp, dred, protect_red, sk, rstprotection, ko, spro);
                        //double Q1;
                        Qpro = ( 4.0 / c ) * sqrt ( Jpro / 100.0 ) * ( aw + 4.0 ) ;
                        Cpro = (spro * spro * Qpro) / (10000.0);
                        Cp = Cpro / 100.0;
                        parsat = 1.8; //parsat=1.5 =>saturation  ; 1.8 => chroma ; 2.5 => colorfullness (personal evaluation : for not)

                        if (chr == -100.0) {
                            chr = -99.8;
                        }

                        if (alg != 2) {
                            Ciecam02::curvecolor (chr, Cp, sres, parsat);
                        } else {
                            Ciecam02::curvecolor (0.0, Cp, sres, parsat);    //chroma
                        }

                        dred = 55.f;
                        protect_red = 30.0f;
                        sk = 1;
                        Color::skinred (Jpro, hpro, sres, Cp, dred, protect_red, sk, rstprotection, ko, Cpro);

                        if (Jpro < 1. && Cpro > 12.) {
                            Cpro = 12.;    //reduce artifacts by "pseudo gamut control CIECAM"
                        } else if (Jpro < 2. && Cpro > 15.) {
                            Cpro = 15.;
                        } else if (Jpro < 4. && Cpro > 30.) {
                            Cpro = 30.;
                        } else if (Jpro < 7. && Cpro > 50.) {
                            Cpro = 50.;
                        }

                        hpro = hpro + hue;

                        if ( hpro < 0.0 ) {
                            hpro += 360.0;    //hue
                        }
                    }

                    if (hasColCurve1) {//curve 1 with Lightness and Brightness
                        if (curveMode == ColorAppearanceParams::TC_MODE_LIGHT) {
                            /*  float Jj=(float) Jpro*327.68;
                            float Jold=Jj;
                            const Lightcurve& userColCurve = static_cast<const Lightcurve&>(customColCurve1);
                              userColCurve.Apply(Jj);
                              Jj=0.7f*(Jj-Jold)+Jold;//divide sensibility
                            */
                            float Jj = (float) Jpro * 327.68f;
                            float Jold = Jj;
                            float Jold100 = (float) Jpro;
                            float redu = 25.f;
                            float reduc = 1.f;
                            const Lightcurve& userColCurveJ1 = static_cast<const Lightcurve&> (customColCurve1);
                            userColCurveJ1.Apply (Jj);

                            if (Jj > Jold) {
                                if (Jj < 65535.f)  {
                                    if (Jold < 327.68f * redu) {
                                        Jj = 0.3f * (Jj - Jold) + Jold;    //divide sensibility
                                    } else        {
                                        reduc = LIM ((100.f - Jold100) / (100.f - redu), 0.f, 1.f);
                                        Jj = 0.3f * reduc * (Jj - Jold) + Jold; //reduct sensibility in highlights
                                    }
                                }
                            } else if (Jj > 10.f) {
                                Jj = 0.8f * (Jj - Jold) + Jold;
                            } else if (Jj >= 0.f) {
                                Jj = 0.90f * (Jj - Jold) + Jold;    // not zero ==>artifacts
                            }


                            Jpro = (double) (Jj / 327.68f);

                            if (Jpro < 1.) {
                                Jpro = 1.;
                            }

                            t1L = true;
                        } else if (curveMode == ColorAppearanceParams::TC_MODE_BRIGHT) {
                            //attention! Brightness curves are open - unlike Lightness or Lab or RGB==> rendering  and algoritms will be different
                            float coef = ((aw + 4.f) * (4.f / c)) / 100.f;
                            float Qq = (float) Qpro * 327.68f * (1.f / coef);
                            float Qold100 = (float) Qpro / coef;

                            float Qold = Qq;
                            float redu = 20.f;
                            float reduc = 1.f;

                            const Brightcurve& userColCurveB1 = static_cast<const Brightcurve&> (customColCurve1);
                            userColCurveB1.Apply (Qq);

                            if (Qq > Qold) {
                                if (Qq < 65535.f)  {
                                    if (Qold < 327.68f * redu) {
                                        Qq = 0.25f * (Qq - Qold) + Qold;    //divide sensibility
                                    } else            {
                                        reduc = LIM ((100.f - Qold100) / (100.f - redu), 0.f, 1.f);
                                        Qq = 0.25f * reduc * (Qq - Qold) + Qold; //reduct sensibility in highlights
                                    }
                                }
                            } else if (Qq > 10.f) {
                                Qq = 0.5f * (Qq - Qold) + Qold;
                            } else if (Qq >= 0.f) {
                                Qq = 0.7f * (Qq - Qold) + Qold;    // not zero ==>artifacts
                            }

                            Qpro = (double) (Qq * (coef) / 327.68f);
                            Jpro = 100.* (Qpro * Qpro) / ((4.0 / c) * (4.0 / c) * (aw + 4.0) * (aw + 4.0));
                            t1B = true;

                            if (Jpro < 1.) {
                                Jpro = 1.;
                            }

                        }
                    }

                    if (hasColCurve2) {//curve 2 with Lightness and Brightness
                        if (curveMode2 == ColorAppearanceParams::TC_MODE_LIGHT) {
                            float Jj = (float) Jpro * 327.68;
                            float Jold = Jj;
                            /*
                            const Lightcurve& userColCurve = static_cast<const Lightcurve&>(customColCurve2);
                                    userColCurve.Apply(Jj);
                                    Jj=0.7f*(Jj-Jold)+Jold;//divide sensibility
                                    */
                            float Jold100 = (float) Jpro;
                            float redu = 25.f;
                            float reduc = 1.f;
                            const Lightcurve& userColCurveJ2 = static_cast<const Lightcurve&> (customColCurve2);
                            userColCurveJ2.Apply (Jj);

                            if (Jj > Jold) {
                                if (Jj < 65535.f)  {
                                    if (Jold < 327.68f * redu) {
                                        Jj = 0.3f * (Jj - Jold) + Jold;    //divide sensibility
                                    } else        {
                                        reduc = LIM ((100.f - Jold100) / (100.f - redu), 0.f, 1.f);
                                        Jj = 0.3f * reduc * (Jj - Jold) + Jold; //reduct sensibility in highlights
                                    }
                                }
                            } else if (Jj > 10.f) {
                                if (!t1L) {
                                    Jj = 0.8f * (Jj - Jold) + Jold;
                                } else {
                                    Jj = 0.4f * (Jj - Jold) + Jold;
                                }
                            } else if (Jj >= 0.f) {
                                if (!t1L) {
                                    Jj = 0.90f * (Jj - Jold) + Jold;    // not zero ==>artifacts
                                } else {
                                    Jj = 0.5f * (Jj - Jold) + Jold;
                                }
                            }

                            Jpro = (double) (Jj / 327.68f);

                            if (Jpro < 1.) {
                                Jpro = 1.;
                            }

                        } else if (curveMode2 == ColorAppearanceParams::TC_MODE_BRIGHT) { //
                            float coef = ((aw + 4.f) * (4.f / c)) / 100.f;
                            float Qq = (float) Qpro * 327.68f * (1.f / coef);
                            float Qold100 = (float) Qpro / coef;

                            float Qold = Qq;
                            float redu = 20.f;
                            float reduc = 1.f;

                            const Brightcurve& userColCurveB2 = static_cast<const Brightcurve&> (customColCurve2);
                            userColCurveB2.Apply (Qq);

                            if (Qq > Qold) {
                                if (Qq < 65535.f)  {
                                    if (Qold < 327.68f * redu) {
                                        Qq = 0.25f * (Qq - Qold) + Qold;    //divide sensibility
                                    } else            {
                                        reduc = LIM ((100.f - Qold100) / (100.f - redu), 0.f, 1.f);
                                        Qq = 0.25f * reduc * (Qq - Qold) + Qold; //reduct sensibility in highlights
                                    }
                                }
                            } else if (Qq > 10.f) {
                                Qq = 0.5f * (Qq - Qold) + Qold;
                            } else if (Qq >= 0.f) {
                                Qq = 0.7f * (Qq - Qold) + Qold;    // not zero ==>artifacts
                            }

                            Qpro = (double) (Qq * (coef) / 327.68f);
                            Jpro = 100.* (Qpro * Qpro) / ((4.0 / c) * (4.0 / c) * (aw + 4.0) * (aw + 4.0));
                            t2B = true;

                            if (t1L) { //to workaround the problem if we modify curve1-lightnees after curve2 brightness(the cat that bites its own tail!) in fact it's another type of curve only for this case
                                coef = 2.f; //adapt Q to J approximation
                                Qq = (float) Qpro * coef;
                                Qold = Qq;
                                const Lightcurve& userColCurveJ1 = static_cast<const Lightcurve&> (customColCurve1);
                                userColCurveJ1.Apply (Qq);
                                Qq = 0.1f * (Qq - Qold) + Qold; //approximative adaptation
                                Qpro = (double) (Qq / coef);
                                Jpro = 100.* (Qpro * Qpro) / ((4.0 / c) * (4.0 / c) * (aw + 4.0) * (aw + 4.0));
                            }

                            if (Jpro < 1.) {
                                Jpro = 1.;
                            }
                        }
                    }

                    if (hasColCurve3) {//curve 3 with chroma saturation colorfullness
                        if (curveMode3 == ColorAppearanceParams::TC_MODE_CHROMA) {
                            double parsat = 0.8; //0.68;
                            double coef = 327.68 / parsat;
                            float Cc = (float) Cpro * coef;
                            float Ccold = Cc;
                            const Chromacurve& userColCurve = static_cast<const Chromacurve&> (customColCurve3);
                            userColCurve.Apply (Cc);
                            float dred = 55.f;
                            float protect_red = 30.0f;
                            float sk = 1;
                            float ko = 1.f / coef;
                            Color::skinred (Jpro, hpro, Cc, Ccold, dred, protect_red, sk, rstprotection, ko, Cpro);

                            if (Jpro < 1. && Cpro > 12.) {
                                Cpro = 12.;    //reduce artifacts by "pseudo gamut control CIECAM"
                            } else if (Jpro < 2. && Cpro > 15.) {
                                Cpro = 15.;
                            } else if (Jpro < 4. && Cpro > 30.) {
                                Cpro = 30.;
                            } else if (Jpro < 7. && Cpro > 50.) {
                                Cpro = 50.;
                            }

                            //  Cpro=Cc/coef;
                        } else if (curveMode3 == ColorAppearanceParams::TC_MODE_SATUR) { //
                            double parsat = 0.8; //0.6
                            double coef = 327.68 / parsat;
                            float Ss = (float) spro * coef;
                            float Sold = Ss;
                            const Saturcurve& userColCurve = static_cast<const Saturcurve&> (customColCurve3);
                            userColCurve.Apply (Ss);
                            Ss = 0.6f * (Ss - Sold) + Sold; //divide sensibility saturation
                            double coe = pow (fl, 0.25);
                            float dred = 100.f; // in C mode
                            float protect_red = 80.0f; // in C mode
                            dred = 100.0 * sqrt ((dred * coe) / Qpro);
                            protect_red = 100.0 * sqrt ((protect_red * coe) / Qpro);
                            int sk = 0;
                            float ko = 1.f / coef;
                            Color::skinred (Jpro, hpro, Ss, Sold, dred, protect_red, sk, rstprotection, ko, spro);
                            Qpro = ( 4.0 / c ) * sqrt ( Jpro / 100.0 ) * ( aw + 4.0 ) ;
                            Cpro = (spro * spro * Qpro) / (10000.0);
                            c1s = 1;

                        } else if (curveMode3 == ColorAppearanceParams::TC_MODE_COLORF) { //
                            double parsat = 0.8; //0.68;
                            double coef = 327.68 / parsat;
                            float Mm = (float) Mpro * coef;
                            float Mold = Mm;
                            const Colorfcurve& userColCurve = static_cast<const Colorfcurve&> (customColCurve3);
                            userColCurve.Apply (Mm);
                            double coe = pow (fl, 0.25);
                            float dred = 100.f; //in C mode
                            float protect_red = 80.0f; // in C mode
                            dred *= coe; //in M mode
                            protect_red *= coe;
                            int sk = 0;
                            float ko = 1.f / coef;
                            Color::skinred (Jpro, hpro, Mm, Mold, dred, protect_red, sk, rstprotection, ko, Mpro);
                            Cpro = Mpro / coe;

                            if (Jpro < 1. && Mpro > 12.*coe) {
                                Mpro = 12.*coe;    //reduce artifacts by "pseudo gamut control CIECAM"
                            } else if (Jpro < 2. && Mpro > 15.*coe) {
                                Mpro = 15.*coe;
                            } else if (Jpro < 4. && Mpro > 30.*coe) {
                                Mpro = 30.*coe;
                            } else if (Jpro < 7. && Mpro > 50.*coe) {
                                Mpro = 50.*coe;
                            }


                            c1co = 1;
                        }
                    }

                    //to retrieve the correct values of variables
                    if (t2B && t1B) {
                        Jpro = (100.0 * Qpro * Qpro) / (wh * wh);    // for brightness curve
                    }

                    if (c1s == 1) {
                        Qpro = ( 4.0 / c ) * sqrt ( Jpro / 100.0 ) * ( aw + 4.0 ) ; //for saturation curve
                        Cpro = (spro * spro * Qpro) / (10000.0);
                    }

                    if (c1co == 1) {
                        double coe = pow (fl, 0.25);   // for colorfullness curve
                        Cpro = Mpro / coe;
                    }

                    //retrieve values C,J...s
                    C = Cpro;
                    J = Jpro;
                    Q = Qpro;
                    M = Mpro;
                    h = hpro;
                    s = spro;

                    if (params->colorappearance.tonecie  || settings->autocielab) { //use pointer for tonemapping with CIECAM and also sharpening , defringe, contrast detail
                        //  if(params->colorappearance.tonecie  || params->colorappearance.sharpcie){//use pointer for tonemapping with CIECAM and also sharpening , defringe, contrast detail
                        float Qred = ( 4.0 / c)  * ( aw + 4.0 ); //estimate Q max if J=100.0
                        ncie->Q_p[i][j] = (float)Q + epsil; //epsil to avoid Q=0
                        ncie->M_p[i][j] = (float)M + epsil;
                        ncie->J_p[i][j] = (float)J + epsil;
                        ncie->h_p[i][j] = (float)h;
                        ncie->C_p[i][j] = (float)C + epsil;
                        //      ncie->s_p[i][j]=(float)s;
                        ncie->sh_p[i][j] = (float) 32768.* (( 4.0 / c ) * sqrt ( J / 100.0 ) * ( aw + 4.0 )) / Qred ;

                        //      ncie->ch_p[i][j]=(float) 327.68*C;
                        if (ncie->Q_p[i][j] < minQ) {
                            minQ = ncie->Q_p[i][j];    //minima
                        }

                        if (ncie->Q_p[i][j] > maxQ) {
                            maxQ = ncie->Q_p[i][j];    //maxima
                        }
                    }

                    if (!params->colorappearance.tonecie  || !settings->autocielab  || !params->epd.enabled ) {

//          if(!params->epd.enabled || !params->colorappearance.tonecie  || !settings->autocielab){
                        //  if(!params->epd.enabled || !params->colorappearance.tonecie  || !params->colorappearance.sharpcie){
                        double brli = 327.;
                        double chsacol = 327.;
                        int libr = 0;
                        int colch = 0;

                        if (curveMode == ColorAppearanceParams::TC_MODE_BRIGHT) {
                            brli = 70.0;
                            libr = 1;
                        } else if (curveMode == ColorAppearanceParams::TC_MODE_LIGHT) {
                            brli = 327.;
                            libr = 0;
                        }

                        if (curveMode3 == ColorAppearanceParams::TC_MODE_CHROMA) {
                            chsacol = 400.;//327.;
                            colch = 0;
                        } else if (curveMode3 == ColorAppearanceParams::TC_MODE_SATUR) {
                            chsacol = 450.0;
                            colch = 1;
                        } else if (curveMode3 == ColorAppearanceParams::TC_MODE_COLORF) {
                            chsacol = 400.;//327.0;
                            colch = 2;
                        }

                        if (ciedata) {
                            // Data for J Q M s and C histograms
                            //update histogram
                            jp = true;

                            if (pW != 1) { //only with improccoordinator
                                int posl;

                                if (libr == 1) {
                                    posl = CLIP ((int) (Q * brli));  //40.0 to 100.0 approximative factor for Q  - 327 for J
                                } else { /*if(libr == 0)*/
                                    posl = CLIP ((int) (J * brli));  //327 for J
                                }

                                hist16JCAM[posl]++;
                            }

                            chropC = true;

                            if (pW != 1) { //only with improccoordinator
                                int posc;

                                if (colch == 0) {
                                    posc = CLIP ((int) (C * chsacol));  //450.0 approximative factor for s    320 for M
                                } else if (colch == 1) {
                                    posc = CLIP ((int) (s * chsacol));
                                } else { /*if(colch == 2)*/
                                    posc = CLIP ((int) (M * chsacol));
                                }

                                hist16_CCAM[posc]++;
                            }
                        }

                        double xx, yy, zz;
                        //double nj, nbbj, ncbj, flj, czj, dj, awj;
                        //process normal==> viewing
                        Ciecam02::jch2xyz_ciecam02 ( xx, yy, zz,
                                                     J,  C, h,
                                                     xw2, yw2,  zw2,
                                                     yb2, la2,
                                                     f2,  c2, nc2, gamu, nj, nbbj, ncbj, flj, czj, dj, awj);
                        x = (float)xx * 655.35;
                        y = (float)yy * 655.35;
                        z = (float)zz * 655.35;
                        float Ll, aa, bb;
                        //convert xyz=>lab
                        Color::XYZ2Lab (x,  y,  z, Ll, aa, bb);
                        lab->L[i][j] = Ll;
                        lab->a[i][j] = aa;
                        lab->b[i][j] = bb;

                        // gamut control in Lab mode; I must study how to do with cIECAM only
                        if (gamu == 1) {
                            float R, G, B;
                            float HH, Lprov1, Chprov1;
                            Lprov1 = lab->L[i][j] / 327.68f;
                            Chprov1 = sqrt (SQR (lab->a[i][j] / 327.68f) + SQR (lab->b[i][j] / 327.68f));
                            HH = atan2 (lab->b[i][j], lab->a[i][j]);

#ifdef _DEBUG
                            bool neg = false;
                            bool more_rgb = false;
                            //gamut control : Lab values are in gamut
                            Color::gamutLchonly (HH, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                            //gamut control : Lab values are in gamut
                            Color::gamutLchonly (HH, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif

                            lab->L[i][j] = Lprov1 * 327.68f;
                            lab->a[i][j] = 327.68f * Chprov1 * cos (HH);
                            lab->b[i][j] = 327.68f * Chprov1 * sin (HH);

                        }
                    }
                }
        }

        // End of parallelization
        if (!params->epd.enabled || !params->colorappearance.tonecie   || !settings->autocielab) { //normal
//if(!params->epd.enabled || !params->colorappearance.tonecie   || !params->colorappearance.sharpcie){//normal

            if (ciedata) {
                //update histogram J
                if (pW != 1) { //only with improccoordinator
                    for (int i = 0; i < 32768; i++) { //
                        if (jp) {
                            float hval = dLcurve[i];
                            int hi = (int) (255.0 * CLIPD (hval)); //
                            histLCAM[hi] += hist16JCAM[i] ;
                        }
                    }
                }

                if (pW != 1) { //only with improccoordinator
                    for (int i = 0; i < 48000; i++) { //
                        if (chropC) {
                            float hvalc = dCcurve[i];
                            int hic = (int) (255.0 * CLIPD (hvalc)); //
                            histCCAM[hic] += hist16_CCAM[i] ;
                        }
                    }
                }
            }
        }

#ifdef _DEBUG

        if (settings->verbose) {
            t2e.set();
            printf ("CIECAM02 performed in %d usec:\n", t2e.etime (t1e));
            //  printf("minc=%f maxc=%f minj=%f maxj=%f\n",minc,maxc,minj,maxj);
        }

#endif

        if (settings->autocielab) {
//if(params->colorappearance.sharpcie) {

//all this treatments reduce artifacts, but can lead to slightly  different results
            if (params->defringe.enabled) if (execsharp) {
                    ImProcFunctions::defringecam (ncie);    //
                }

//if(params->dirpyrequalizer.enabled) if(execsharp) {
            if (params->dirpyrequalizer.enabled) {
                if (params->dirpyrequalizer.gamutlab  /*&& execsharp*/) {
                    float b_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[0]) / 100.0f;
                    float t_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[1]) / 100.0f;
                    float b_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[2]) / 100.0f;
                    float t_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[3]) / 100.0f;

                    float artifact = (float) settings->artifact_cbdl;

                    if (artifact > 6.f) {
                        artifact = 6.f;
                    }

                    if (artifact < 0.f) {
                        artifact = 1.f;
                    }

                    float chrom = 50.f;
                    {
                        int hotbad = 0;
                        ImProcFunctions::badpixcam (ncie, artifact, 5, 2, b_l, t_l, t_r, b_r, params->dirpyrequalizer.skinprotect, chrom, hotbad);      //enabled remove artifacts for cbDL
                    }
                }
            }

            if (params->colorappearance.badpixsl > 0) if (execsharp) {
                    int mode = params->colorappearance.badpixsl;
                    ImProcFunctions::badpixcam (ncie, 3.4, 5, mode, 0, 0, 0, 0, 0, 0, 1);//for bad pixels CIECAM
                }

            if (params->sharpenMicro.enabled)if (execsharp) {
                    ImProcFunctions::MLmicrocontrastcam (ncie);
                }

            if (params->sharpening.enabled)
                if (execsharp) {
                    float **buffer = lab->L; // We can use the L-buffer from lab as buffer to save some memory
                    ImProcFunctions::sharpeningcam (ncie, buffer); // sharpening adapted to CIECAM
                }

//if(params->dirpyrequalizer.enabled) if(execsharp) {
            if (params->dirpyrequalizer.enabled /*&& (execsharp)*/) {

//  if     (params->dirpyrequalizer.algo=="FI") choice=0;
//  else if(params->dirpyrequalizer.algo=="LA") choice=1;
                if (rtt == 1) {
                    float b_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[0]) / 100.0f;
                    float t_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[1]) / 100.0f;
                    float b_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[2]) / 100.0f;
                    float t_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[3]) / 100.0f;
                    int choice = 0; //not disabled in case of ! always 0
                    dirpyr_equalizercam (ncie, ncie->sh_p, ncie->sh_p, ncie->W, ncie->H, ncie->h_p, ncie->C_p, params->dirpyrequalizer.mult, params->dirpyrequalizer.threshold,  params->dirpyrequalizer.skinprotect, true, params->dirpyrequalizer.gamutlab, b_l, t_l, t_r, b_r, choice, scale);   //contrast by detail adapted to CIECAM
                }
            }

            float Qredi = ( 4.0 / c_)  * ( a_w + 4.0 );
            float co_e = (pow (f_l, 0.25f));

#ifndef _DEBUG
            #pragma omp parallel default(shared) firstprivate(height,width, Qredi,a_w,c_)
#endif
            {
#ifndef _DEBUG
                #pragma omp for schedule(dynamic, 10)
#endif

                for (int i = 0; i < height; i++) // update CieImages with new values after sharpening, defringe, contrast by detail level
                    for (int j = 0; j < width; j++) {
                        float interm = Qredi * ncie->sh_p[i][j] / (32768.f);
                        ncie->J_p[i][j] = 100.0 * interm * interm / ((a_w + 4.) * (a_w + 4.) * (4. / c_) * (4. / c_));
                        ncie->Q_p[i][j] = ( 4.0 / c_)  * ( a_w + 4.0 ) *  sqrt (ncie->J_p[i][j] / 100.f);
                        ncie->M_p[i][j] = ncie->C_p[i][j] * co_e;
                    }
            }
        }

        if ((params->colorappearance.tonecie || (params->colorappearance.tonecie && params->epd.enabled)) || (params->sharpening.enabled && settings->autocielab)
                || (params->dirpyrequalizer.enabled && settings->autocielab) || (params->defringe.enabled && settings->autocielab)  || (params->sharpenMicro.enabled && settings->autocielab)
                ||  (params->colorappearance.badpixsl > 0 && settings->autocielab)) {

            if (params->epd.enabled  && params->colorappearance.tonecie  && algepd) {
                ImProcFunctions::EPDToneMapCIE (ncie, a_w, c_, w_h, width, height, begh, endh, minQ, maxQ, Iterates, scale );
            }

            //EPDToneMapCIE adapted to CIECAM


#ifndef _DEBUG
            #pragma omp parallel default(shared) firstprivate(lab,xw2,yw2,zw2,chr,yb,la2,yb2, height,width,begh, endh, nc2,f2,c2, gamu, highlight,pW)
#endif
            {
                TMatrix wiprofa = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.working);
                double wipa[3][3] = {
                    {wiprofa[0][0], wiprofa[0][1], wiprofa[0][2]},
                    {wiprofa[1][0], wiprofa[1][1], wiprofa[1][2]},
                    {wiprofa[2][0], wiprofa[2][1], wiprofa[2][2]}
                };


#ifndef _DEBUG
                #pragma omp for schedule(dynamic, 10)
#endif

                for (int i = 0; i < height; i++) // update CIECAM with new values after tone-mapping

                    //  for (int i=begh; i<endh; i++)
                    for (int j = 0; j < width; j++) {
                        double xx, yy, zz;
                        float x, y, z;
                        const float eps = 0.0001;
                        float co_e = (pow (f_l, 0.25f)) + eps;

                        //      if(params->epd.enabled) ncie->J_p[i][j]=(100.0* ncie->Q_p[i][j]*ncie->Q_p[i][j])/(w_h*w_h);
                        if (params->epd.enabled) {
                            ncie->J_p[i][j] = (100.0 * ncie->Q_p[i][j] * ncie->Q_p[i][j]) / SQR ((4. / c) * (aw + 4.));
                        }

                        ncie->C_p[i][j] = (ncie->M_p[i][j]) / co_e;
                        //show histogram in CIECAM mode (Q,J, M,s,C)
                        double brli = 327.;
                        double chsacol = 327.;
                        int libr = 0;
                        int colch = 0;

                        if (curveMode == ColorAppearanceParams::TC_MODE_BRIGHT) {
                            brli = 70.0;
                            libr = 1;
                        } else if (curveMode == ColorAppearanceParams::TC_MODE_LIGHT) {
                            brli = 327.;
                            libr = 0;
                        }

                        if (curveMode3 == ColorAppearanceParams::TC_MODE_CHROMA) {
                            chsacol = 400.;//327.;
                            colch = 0;
                        } else if (curveMode3 == ColorAppearanceParams::TC_MODE_SATUR) {
                            chsacol = 450.0;
                            colch = 1;
                        } else if (curveMode3 == ColorAppearanceParams::TC_MODE_COLORF) {
                            chsacol = 400.;//327.0;
                            colch = 2;
                        }

                        if (ciedata) {
                            // Data for J Q M s and C histograms
                            //update histogram
                            jp = true;

                            if (pW != 1) { //only with improccoordinator
                                int posl;

                                if (libr == 1) {
                                    posl = CLIP ((int) (ncie->Q_p[i][j] * brli));  //40.0 to 100.0 approximative factor for Q  - 327 for J
                                } else { /*if(libr == 0)*/
                                    posl = CLIP ((int) (ncie->J_p[i][j] * brli));  //327 for J
                                }

                                hist16JCAM[posl]++;
                            }

                            chropC = true;

                            if (pW != 1) { //only with improccoordinator
                                int posc;

                                if (colch == 0) {
                                    posc = CLIP ((int) (ncie->C_p[i][j] * chsacol));  //450.0 approximative factor for s    320 for M
                                } else if (colch == 1) {
                                    float sa_t = 100.f * sqrt (ncie->C_p[i][j] / ncie->Q_p[i][j]);   //Q_p always > 0
                                    posc = CLIP ((int) (sa_t * chsacol));
                                } else { /*if(colch == 2)*/
                                    posc = CLIP ((int) (ncie->M_p[i][j] * chsacol));
                                }

                                hist16_CCAM[posc]++;
                            }
                        }

                        //end histograms
                        //  double nd, nbbd, ncbd, fld, czd, dd, awd;
                        Ciecam02::jch2xyz_ciecam02 ( xx, yy, zz,
                                                     ncie->J_p[i][j],  ncie->C_p[i][j], ncie->h_p[i][j],
                                                     xw2, yw2,  zw2,
                                                     yb2, la2,
                                                     f2,  c2, nc2, gamu, nj, nbbj, ncbj, flj, czj, dj, awj);
                        x = (float)xx * 655.35;
                        y = (float)yy * 655.35;
                        z = (float)zz * 655.35;
                        float Ll, aa, bb;
                        //convert xyz=>lab
                        Color::XYZ2Lab (x,  y,  z, Ll, aa, bb);
                        lab->L[i][j] = Ll;
                        lab->a[i][j] = aa;
                        lab->b[i][j] = bb;

                        if (gamu == 1) {
                            float R, G, B;
                            float HH, Lprov1, Chprov1;
                            Lprov1 = lab->L[i][j] / 327.68f;
                            Chprov1 = sqrt (SQR (lab->a[i][j] / 327.68f) + SQR (lab->b[i][j] / 327.68f));
                            HH = atan2 (lab->b[i][j], lab->a[i][j]);

#ifdef _DEBUG
                            bool neg = false;
                            bool more_rgb = false;
                            //gamut control : Lab values are in gamut
                            Color::gamutLchonly (HH, Lprov1, Chprov1, R, G, B, wipa, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                            //gamut control : Lab values are in gamut
                            Color::gamutLchonly (HH, Lprov1, Chprov1, R, G, B, wipa, highlight, 0.15f, 0.96f);
#endif

                            lab->L[i][j] = Lprov1 * 327.68f;
                            lab->a[i][j] = 327.68f * Chprov1 * cos (HH);
                            lab->b[i][j] = 327.68f * Chprov1 * sin (HH);
                        }
                    }


            }

            //end parallelization
            //show CIECAM histograms
            if (ciedata) {
                //update histogram J and Q
                if (pW != 1) { //only with improccoordinator
                    for (int i = 0; i < 32768; i++) { //
                        if (jp) {
                            float hval = dLcurve[i];
                            int hi = (int) (255.0 * CLIPD (hval)); //
                            histLCAM[hi] += hist16JCAM[i] ;
                        }
                    }
                }

                //update color histogram M,s,C
                if (pW != 1) { //only with improccoordinator
                    for (int i = 0; i < 48000; i++) { //
                        if (chropC) {
                            float hvalc = dCcurve[i];
                            int hic = (int) (255.0 * CLIPD (hvalc)); //
                            histCCAM[hic] += hist16_CCAM[i] ;
                        }
                    }
                }
            }

        }

    }
}
//end CIECAM


// Copyright (c) 2012 Jacques Desmis <jdesmis@gmail.com>
void ImProcFunctions::ciecam_02float (CieImage* ncie, float adap, int begh, int endh, int pW, int pwb, LabImage* lab, const ProcParams* params,
                                      const ColorAppearance & customColCurve1, const ColorAppearance & customColCurve2, const ColorAppearance & customColCurve3,
                                      LUTu & histLCAM, LUTu & histCCAM, LUTf & CAMBrightCurveJ, LUTf & CAMBrightCurveQ, float &mean, int Iterates, int scale, bool execsharp, float &d, float &dj, float &yb, int rtt)
{
    if (params->colorappearance.enabled) {

#ifdef _DEBUG
        MyTime t1e, t2e;
        t1e.set();
#endif

        //preparate for histograms CIECAM
        LUTu hist16JCAM;
        LUTu hist16_CCAM;

        if (pW != 1 && params->colorappearance.datacie) { //only with improccoordinator
            hist16JCAM (32768);
            hist16JCAM.clear();
            hist16_CCAM (48000);
            hist16_CCAM.clear();
        }

        //end preparate histogram
        int width = lab->W, height = lab->H;
        float minQ = 10000.f;
        float maxQ = -1000.f;
        float Yw;
        Yw = 1.0;
        double Xw, Zw;
        float f = 0.f, nc = 0.f, la, c = 0.f, xw, yw, zw, f2 = 1.f, c2 = 1.f, nc2 = 1.f, yb2;
        float fl, n, nbb, ncb, aw; //d
        float xwd, ywd, zwd, xws, yws, zws;
        int alg = 0;
        bool algepd = false;
        double Xwout, Zwout;
        double Xwsc, Zwsc;

        const bool epdEnabled = params->epd.enabled;
        bool ciedata = (params->colorappearance.datacie && pW != 1) && ! ((params->colorappearance.tonecie && (epdEnabled)) || (params->sharpening.enabled && settings->autocielab && execsharp)
                       || (params->dirpyrequalizer.enabled && settings->autocielab) || (params->defringe.enabled && settings->autocielab)  || (params->sharpenMicro.enabled && settings->autocielab)
                       || (params->impulseDenoise.enabled && settings->autocielab) ||  (params->colorappearance.badpixsl > 0 && settings->autocielab));

        ColorTemp::temp2mulxyz (params->wb.temperature, params->wb.green, params->wb.method, Xw, Zw); //compute white Xw Yw Zw  : white current WB
        ColorTemp::temp2mulxyz (params->colorappearance.tempout, params->colorappearance.greenout, "Custom", Xwout, Zwout);
        ColorTemp::temp2mulxyz (params->colorappearance.tempsc, params->colorappearance.greensc, "Custom", Xwsc, Zwsc);

        //viewing condition for surrsrc
        if (params->colorappearance.surrsrc == "Average") {
            f  = 1.00f;
            c  = 0.69f;
            nc = 1.00f;
        } else if (params->colorappearance.surrsrc == "Dim") {
            f  = 0.9f;
            c  = 0.59f;
            nc = 0.9f;
        } else if (params->colorappearance.surrsrc == "Dark") {
            f  = 0.8f;
            c  = 0.525f;
            nc = 0.8f;
        } else if (params->colorappearance.surrsrc == "ExtremelyDark") {
            f  = 0.8f;
            c  = 0.41f;
            nc = 0.8f;
        }


        //viewing condition for surround
        if (params->colorappearance.surround == "Average") {
            f2 = 1.0f, c2 = 0.69f, nc2 = 1.0f;
        } else if (params->colorappearance.surround == "Dim") {
            f2  = 0.9f;
            c2  = 0.59f;
            nc2 = 0.9f;
        } else if (params->colorappearance.surround == "Dark") {
            f2  = 0.8f;
            c2  = 0.525f;
            nc2 = 0.8f;
        } else if (params->colorappearance.surround == "ExtremelyDark") {
            f2  = 0.8f;
            c2  = 0.41f;
            nc2 = 0.8f;
        }

        /*
                //scene condition for surround
                if (params->colorappearance.surrsource)  {
                    f  = 0.85f;    // if user => source image has surround very dark
                    c  = 0.55f;
                    nc = 0.85f;
                }
        */
        //with which algorithm
        if     (params->colorappearance.algo == "JC") {
            alg = 0;
        } else if (params->colorappearance.algo == "JS") {
            alg = 1;
        } else if (params->colorappearance.algo == "QM")  {
            alg = 2;
            algepd = true;
        } else { /*if(params->colorappearance.algo == "ALL")*/
            alg = 3;
            algepd = true;
        }

        xwd = 100.f * Xwout;
        zwd = 100.f * Zwout;
        ywd = 100.f / params->colorappearance.greenout;//approximation to simplify

        xws = 100.f * Xwsc;
        zws = 100.f * Zwsc;
        yws = 100.f / params->colorappearance.greensc;//approximation to simplify


        /*
                //settings white point of output device - or illuminant viewing
                if (settings->viewingdevice == 0) {
                    xwd = 96.42f;    //5000K
                    ywd = 100.0f;
                    zwd = 82.52f;
                } else if (settings->viewingdevice == 1) {
                    xwd = 95.68f;    //5500
                    ywd = 100.0f;
                    zwd = 92.15f;
                } else if (settings->viewingdevice == 2) {
                    xwd = 95.24f;    //6000
                    ywd = 100.0f;
                    zwd = 100.81f;
                } else if (settings->viewingdevice == 3)  {
                    xwd = 95.04f;    //6500
                    ywd = 100.0f;
                    zwd = 108.88f;
                } else if (settings->viewingdevice == 4)  {
                    xwd = 109.85f;    //tungsten
                    ywd = 100.0f;
                    zwd = 35.58f;
                } else if (settings->viewingdevice == 5)  {
                    xwd = 99.18f;    //fluo F2
                    ywd = 100.0f;
                    zwd = 67.39f;
                } else if (settings->viewingdevice == 6)  {
                    xwd = 95.04f;    //fluo F7
                    ywd = 100.0f;
                    zwd = 108.75f;
                } else {
                    xwd = 100.96f;    //fluo F11
                    ywd = 100.0f;
                    zwd = 64.35f;
                }
        */
        yb2 = params->colorappearance.ybout;
        /*
                //settings mean Luminance Y of output device or viewing
                if (settings->viewingdevicegrey == 0) {
                    yb2 = 5.0f;
                } else if (settings->viewingdevicegrey == 1) {
                    yb2 = 10.0f;
                } else if (settings->viewingdevicegrey == 2) {
                    yb2 = 15.0f;
                } else if (settings->viewingdevicegrey == 3) {
                    yb2 = 18.0f;
                } else if (settings->viewingdevicegrey == 4) {
                    yb2 = 23.0f;
                } else if (settings->viewingdevicegrey == 5)  {
                    yb2 = 30.0f;
                } else {
                    yb2 = 40.0f;
                }
        */
        //La and la2 = ambiant luminosity scene and viewing
        la = float (params->colorappearance.adapscen);

        if (pwb == 2) {
            if (params->colorappearance.autoadapscen) {
                la = adap;
            }
        }

        /*
                if (alg >= 2 && la < 200.f) {
                    la = 200.f;
                }
        */
        const float la2 = float (params->colorappearance.adaplum);

        // level of adaptation
        const float deg = (params->colorappearance.degree) / 100.0f;
        const float pilot = params->colorappearance.autodegree ? 2.0f : deg;

        const float degout = (params->colorappearance.degreeout) / 100.0f;
        const float pilotout = params->colorappearance.autodegreeout ? 2.0f : degout;

        //algoritm's params
        float chr = 0.f;

        if (alg == 0 || alg == 3) {
            chr = params->colorappearance.chroma;

            if (chr == -100.0f) {
                chr = -99.8f;
            }
        }

        float schr = 0.f;

        if (alg == 3 || alg == 1) {
            schr = params->colorappearance.schroma;

            if (schr > 0.0) {
                schr = schr / 2.0f;    //divide sensibility for saturation
            }

            if (alg == 3) {
                if (schr == -100.0f) {
                    schr = -99.f;
                }

                if (schr == 100.0f) {
                    schr = 98.f;
                }
            } else {
                if (schr == -100.0f) {
                    schr = -99.8f;
                }
            }
        }

        float mchr = 0.f;

        if (alg == 3 || alg == 2) {
            mchr = params->colorappearance.mchroma;

            if (mchr == -100.0f) {
                mchr = -99.8f ;
            }

            if (mchr == 100.0f) {
                mchr = 99.9f;
            }
        }

        const float hue = params->colorappearance.colorh;
        const float rstprotection = 100. - params->colorappearance.rstprotection;

        // extracting datas from 'params' to avoid cache flush (to be confirmed)
        const ColorAppearanceParams::eTCModeId curveMode = params->colorappearance.curveMode;
        const bool hasColCurve1 = bool (customColCurve1);
        const bool t1L = hasColCurve1 && curveMode == ColorAppearanceParams::TC_MODE_LIGHT;

        const ColorAppearanceParams::eTCModeId curveMode2 = params->colorappearance.curveMode2;
        const bool hasColCurve2 = bool (customColCurve2);

        const ColorAppearanceParams::eCTCModeId curveMode3 = params->colorappearance.curveMode3;
        const bool hasColCurve3 = bool (customColCurve3);

        bool needJ = (alg == 0 || alg == 1 || alg == 3);
        bool needQ = (alg == 2 || alg == 3);
        LUTu hist16J;
        LUTu hist16Q;

        if ((needJ && CAMBrightCurveJ.dirty) || (needQ && CAMBrightCurveQ.dirty) || (std::isnan (mean) && settings->viewinggreySc != 0)) {

            if (needJ) {
                hist16J (32768);
                hist16J.clear();
            }

            if (needQ) {
                hist16Q (32768);
                hist16Q.clear();
            }

            float sum = 0.f;
            float sumQ = 0.f;

#ifdef _OPENMP
            const int numThreads = min (max (width * height / 65536, 1), omp_get_max_threads());
            #pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
            {
                LUTu hist16Jthr;
                LUTu hist16Qthr;

                if (needJ) {
                    hist16Jthr (hist16J.getSize());
                    hist16Jthr.clear();
                }

                if (needQ) {
                    hist16Qthr (hist16Q.getSize());
                    hist16Qthr.clear();
                }

                #pragma omp for reduction(+:sum,sumQ)

                for (int i = 0; i < height; i++)
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

                        if (needJ) {
                            hist16Jthr[ (int) ((koef * lab->L[i][j]))]++;  //evaluate histogram luminance L # J
                        }

                        //estimation of wh only with La
                        float whestim = 500.f;

                        if (la < 200.f) {
                            whestim = 200.f;
                        } else if (la < 2500.f) {
                            whestim = 350.f;
                        } else {
                            whestim = 500.f;
                        }

                        if (needQ) {
                            hist16Qthr[CLIP ((int) (32768.f * sqrt ((koef * (lab->L[i][j])) / 32768.f)))]++;  //for brightness Q : approximation for Q=wh*sqrt(J/100)  J not equal L
                            //perhaps  needs to introduce whestim ??
                            //   hist16Qthr[ (int) (sqrtf ((koef * (lab->L[i][j])) * 32768.f))]++;  //for brightness Q : approximation for Q=wh*sqrt(J/100)  J not equal L
                        }

                        sum += koef * lab->L[i][j]; //evaluate mean J to calculate Yb
                        sumQ += whestim * sqrt ((koef * (lab->L[i][j])) / 32768.f);
                        //can be used in case of...
                    }

                #pragma omp critical
                {
                    if (needJ) {
                        hist16J += hist16Jthr;
                    }

                    if (needQ) {
                        hist16Q += hist16Qthr;
                    }

                }
                float meanQ;

                if (std::isnan (mean)) {
                    mean = (sum / ((height) * width)) / 327.68f; //for Yb  for all image...if one day "pipette" we can adapt Yb for each zone
                    meanQ = (sumQ / ((height) * width));//in case of

                }
            }


            //evaluate lightness, contrast
        }



        //  if (settings->viewinggreySc == 0) { //auto
        if (params->colorappearance.autoybscen  &&  pwb == 2) {//auto

            if     (mean < 15.f) {
                yb = 3.0f;
            } else if (mean < 30.f) {
                yb = 5.0f;
            } else if (mean < 40.f) {
                yb = 10.0f;
            } else if (mean < 45.f) {
                yb = 15.0f;
            } else if (mean < 50.f) {
                yb = 18.0f;
            } else if (mean < 55.f) {
                yb = 23.0f;
            } else if (mean < 60.f) {
                yb = 30.0f;
            } else if (mean < 70.f) {
                yb = 40.0f;
            } else if (mean < 80.f) {
                yb = 60.0f;
            } else if (mean < 90.f) {
                yb = 80.0f;
            } else {
                yb = 90.0f;
            }

//        } else if (settings->viewinggreySc == 1) {
        } else {
            yb =  (float) params->colorappearance.ybscen;
        }

        const bool highlight = params->toneCurve.hrenabled; //Get the value if "highlight reconstruction" is activated

        const int gamu = (params->colorappearance.gamut) ? 1 : 0;
        xw = 100.0f * Xw;
        yw = 100.0f * Yw;
        zw = 100.0f * Zw;
        float xw1 = 0.f, yw1 = 0.f, zw1 = 0.f, xw2 = 0.f, yw2 = 0.f, zw2 = 0.f;

        // settings of WB: scene and viewing
        if (params->colorappearance.wbmodel == "RawT") {
            xw1 = 96.46f;    //use RT WB; CAT 02 is used for output device (see prefreneces)
            yw1 = 100.0f;
            zw1 = 82.445f;
            xw2 = xwd;
            yw2 = ywd;
            zw2 = zwd;
        } else if (params->colorappearance.wbmodel == "RawTCAT02") {
            xw1 = xw;    // Settings RT WB are used for CAT02 => mix , CAT02 is use for output device (screen: D50 D65, projector: lamp, LED) see preferences
            yw1 = yw;
            zw1 = zw;
            xw2 = xwd;
            yw2 = ywd;
            zw2 = zwd;
        } else if (params->colorappearance.wbmodel == "free") {
            xw1 = xws;    // free temp and green
            yw1 = yws;
            zw1 = zws;
            xw2 = xwd;
            yw2 = ywd;
            zw2 = zwd;
        }


        float cz, wh, pfl;
        Ciecam02::initcam1float (gamu, yb, pilot, f, la, xw, yw, zw, n, d, nbb, ncb, cz, aw, wh, pfl, fl, c);
        //printf ("wh=%f \n", wh);

        const float pow1 = pow_F ( 1.64f - pow_F ( 0.29f, n ), 0.73f );
        float nj, nbbj, ncbj, czj, awj, flj;
        Ciecam02::initcam2float (gamu, yb2, pilotout, f2,  la2,  xw2,  yw2,  zw2, nj, dj, nbbj, ncbj, czj, awj, flj);
        const float reccmcz = 1.f / (c2 * czj);
        const float pow1n = pow_F ( 1.64f - pow_F ( 0.29f, nj ), 0.73f );

        const float epsil = 0.0001f;
        const float w_h = wh + epsil;
        const float coefQ = 32767.f / wh;
        const float a_w = aw;
        const float c_ = c;
        const float f_l = fl;
        const float coe = pow_F (fl, 0.25f);
        const float QproFactor = ( 0.4f / c ) * ( aw + 4.0f ) ;
        const bool LabPassOne = ! ((params->colorappearance.tonecie && (epdEnabled)) || (params->sharpening.enabled && settings->autocielab && execsharp)
                                   || (params->dirpyrequalizer.enabled && settings->autocielab) || (params->defringe.enabled && settings->autocielab)  || (params->sharpenMicro.enabled && settings->autocielab)
                                   || (params->impulseDenoise.enabled && settings->autocielab) ||  (params->colorappearance.badpixsl > 0 && settings->autocielab));
        //printf("coQ=%f\n", coefQ);

        if (needJ) {
            if (!CAMBrightCurveJ) {
                CAMBrightCurveJ (32768, LUT_CLIP_ABOVE);
            }

            if (CAMBrightCurveJ.dirty) {
                Ciecam02::curveJfloat (params->colorappearance.jlight, params->colorappearance.contrast, hist16J, CAMBrightCurveJ);//lightness and contrast J
                CAMBrightCurveJ /= 327.68f;
                CAMBrightCurveJ.dirty = false;
            }
        }

        if (needQ) {
            if (!CAMBrightCurveQ) {
                CAMBrightCurveQ (32768, LUT_CLIP_ABOVE);
            }

            if (CAMBrightCurveQ.dirty) {
                Ciecam02::curveJfloat (params->colorappearance.qbright, params->colorappearance.qcontrast, hist16Q, CAMBrightCurveQ);//brightness and contrast Q
                //  CAMBrightCurveQ /= coefQ;
                CAMBrightCurveQ.dirty = false;
            }
        }


        //matrix for current working space
        TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.working);
        const float wip[3][3] = {
            { (float)wiprof[0][0], (float)wiprof[0][1], (float)wiprof[0][2]},
            { (float)wiprof[1][0], (float)wiprof[1][1], (float)wiprof[1][2]},
            { (float)wiprof[2][0], (float)wiprof[2][1], (float)wiprof[2][2]}
        };

#ifdef __SSE2__
        int bufferLength = ((width + 3) / 4) * 4; // bufferLength has to be a multiple of 4
#endif
#ifndef _DEBUG
        #pragma omp parallel
#endif
        {
            float minQThr = 10000.f;
            float maxQThr = -1000.f;
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

                vfloat c655d35 = F2V (655.35f);

                for (k = 0; k < width - 3; k += 4) {
                    Color::Lab2XYZ (LVFU (lab->L[i][k]), LVFU (lab->a[i][k]), LVFU (lab->b[i][k]), x, y, z);
                    x = x / c655d35;
                    y = y / c655d35;
                    z = z / c655d35;
                    Ciecam02::xyz2jchqms_ciecam02float ( J, C,  h,
                                                         Q,  M,  s, F2V (aw), F2V (fl), F2V (wh),
                                                         x,  y,  z,
                                                         F2V (xw1), F2V (yw1),  F2V (zw1),
                                                         F2V (c),  F2V (nc), F2V (pow1), F2V (nbb), F2V (ncb), F2V (pfl), F2V (cz), F2V (d));
                    STVF (Jbuffer[k], J);
                    STVF (Cbuffer[k], C);
                    STVF (hbuffer[k], h);
                    STVF (Qbuffer[k], Q);
                    STVF (Mbuffer[k], M);
                    STVF (sbuffer[k], s);
                }

                for (; k < width; k++) {
                    float L = lab->L[i][k];
                    float a = lab->a[i][k];
                    float b = lab->b[i][k];
                    float x, y, z;
                    //convert Lab => XYZ
                    Color::Lab2XYZ (L, a, b, x, y, z);
                    x = x / 655.35f;
                    y = y / 655.35f;
                    z = z / 655.35f;
                    float J, C, h, Q, M, s;
                    Ciecam02::xyz2jchqms_ciecam02float ( J, C,  h,
                                                         Q,  M,  s, aw, fl, wh,
                                                         x,  y,  z,
                                                         xw1, yw1,  zw1,
                                                         c,  nc, gamu, pow1, nbb, ncb, pfl, cz, d);
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
                    Color::Lab2XYZ (L, a, b, x1, y1, z1);
                    x = (float)x1 / 655.35f;
                    y = (float)y1 / 655.35f;
                    z = (float)z1 / 655.35f;
                    //process source==> normal
                    Ciecam02::xyz2jchqms_ciecam02float ( J, C,  h,
                                                         Q,  M,  s, aw, fl, wh,
                                                         x,  y,  z,
                                                         xw1, yw1,  zw1,
                                                         c,  nc, gamu, pow1, nbb, ncb, pfl, cz, d);
#endif
                    float Jpro, Cpro, hpro, Qpro, Mpro, spro;
                    Jpro = J;
                    Cpro = C;
                    hpro = h;
                    Qpro = Q;
                    Mpro = M;
                    spro = s;

                    // we cannot have all algorithms with all chroma curves
                    if (alg == 0) {
                        Jpro = CAMBrightCurveJ[Jpro * 327.68f]; //lightness CIECAM02 + contrast
                        Qpro = QproFactor * sqrtf (Jpro);
                        float Cp = (spro * spro * Qpro) / (1000000.f);
                        Cpro = Cp * 100.f;
                        float sres;
                        Ciecam02::curvecolorfloat (chr, Cp, sres, 1.8f);
                        Color::skinredfloat (Jpro, hpro, sres, Cp, 55.f, 30.f, 1, rstprotection, 100.f, Cpro);
                    } else if (alg == 1)  {
                        // Lightness saturation
                        Jpro = CAMBrightCurveJ[Jpro * 327.68f]; //lightness CIECAM02 + contrast
                        float sres;
                        float Sp = spro / 100.0f;
                        float parsat = 1.5f; //parsat=1.5 =>saturation  ; 1.8 => chroma ; 2.5 => colorfullness (personal evaluation)
                        Ciecam02::curvecolorfloat (schr, Sp, sres, parsat);
                        float dred = 100.f; // in C mode
                        float protect_red = 80.0f; // in C mode
                        dred = 100.0f * sqrtf ((dred * coe) / Qpro);
                        protect_red = 100.0f * sqrtf ((protect_red * coe) / Qpro);
                        Color::skinredfloat (Jpro, hpro, sres, Sp, dred, protect_red, 0, rstprotection, 100.f, spro);
                        Qpro = QproFactor * sqrtf (Jpro);
                        Cpro = (spro * spro * Qpro) / (10000.0f);
                    } else if (alg == 2) {
                        //printf("Qp0=%f ", Qpro);

                        Qpro = CAMBrightCurveQ[ (float) (Qpro * coefQ)] / coefQ; //brightness and contrast
                        //printf("Qpaf=%f ", Qpro);

                        float Mp, sres;
                        Mp = Mpro / 100.0f;
                        Ciecam02::curvecolorfloat (mchr, Mp, sres, 2.5f);
                        float dred = 100.f; //in C mode
                        float protect_red = 80.0f; // in C mode
                        dred *= coe; //in M mode
                        protect_red *= coe; //M mode
                        Color::skinredfloat (Jpro, hpro, sres, Mp, dred, protect_red, 0, rstprotection, 100.f, Mpro);
                        Jpro = SQR ((10.f * Qpro) / wh);
                        Cpro = Mpro / coe;
                        Qpro = (Qpro == 0.f ? epsil : Qpro); // avoid division by zero
                        spro = 100.0f * sqrtf ( Mpro / Qpro );
                    } else { /*if(alg == 3) */
                        Qpro = CAMBrightCurveQ[ (float) (Qpro * coefQ)] / coefQ; //brightness and contrast

                        float Mp, sres;
                        Mp = Mpro / 100.0f;
                        Ciecam02::curvecolorfloat (mchr, Mp, sres, 2.5f);
                        float dred = 100.f; //in C mode
                        float protect_red = 80.0f; // in C mode
                        dred *= coe; //in M mode
                        protect_red *= coe; //M mode
                        Color::skinredfloat (Jpro, hpro, sres, Mp, dred, protect_red, 0, rstprotection, 100.f, Mpro);
                        Jpro = SQR ((10.f * Qpro) / wh);
                        Cpro = Mpro / coe;
                        Qpro = (Qpro == 0.f ? epsil : Qpro); // avoid division by zero
                        spro = 100.0f * sqrtf ( Mpro / Qpro );

                        if (Jpro > 99.9f) {
                            Jpro = 99.9f;
                        }

                        Jpro = CAMBrightCurveJ[ (float) (Jpro * 327.68f)]; //lightness CIECAM02 + contrast
                        float Sp = spro / 100.0f;
                        Ciecam02::curvecolorfloat (schr, Sp, sres, 1.5f);
                        dred = 100.f; // in C mode
                        protect_red = 80.0f; // in C mode
                        dred = 100.0f * sqrtf ((dred * coe) / Q);
                        protect_red = 100.0f * sqrtf ((protect_red * coe) / Q);
                        Color::skinredfloat (Jpro, hpro, sres, Sp, dred, protect_red, 0, rstprotection, 100.f, spro);
                        Qpro = QproFactor * sqrtf (Jpro);
                        float Cp = (spro * spro * Qpro) / (1000000.f);
                        Cpro = Cp * 100.f;
                        Ciecam02::curvecolorfloat (chr, Cp, sres, 1.8f);
                        Color::skinredfloat (Jpro, hpro, sres, Cp, 55.f, 30.f, 1, rstprotection, 100.f, Cpro);
// disabled this code, Issue 2690
//              if(Jpro < 1.f && Cpro > 12.f) Cpro=12.f;//reduce artifacts by "pseudo gamut control CIECAM"
//              else if(Jpro < 2.f && Cpro > 15.f) Cpro=15.f;
//              else if(Jpro < 4.f && Cpro > 30.f) Cpro=30.f;
//              else if(Jpro < 7.f && Cpro > 50.f) Cpro=50.f;
                        hpro = hpro + hue;

                        if ( hpro < 0.0f ) {
                            hpro += 360.0f;    //hue
                        }
                    }

                    if (hasColCurve1) {//curve 1 with Lightness and Brightness
                        if (curveMode == ColorAppearanceParams::TC_MODE_LIGHT) {
                            float Jj = (float) Jpro * 327.68f;
                            float Jold = Jj;
                            float Jold100 = (float) Jpro;
                            float redu = 25.f;
                            float reduc = 1.f;
                            const Lightcurve& userColCurveJ1 = static_cast<const Lightcurve&> (customColCurve1);
                            userColCurveJ1.Apply (Jj);

                            if (Jj > Jold) {
                                if (Jj < 65535.f)  {
                                    if (Jold < 327.68f * redu) {
                                        Jj = 0.3f * (Jj - Jold) + Jold;    //divide sensibility
                                    } else        {
                                        reduc = LIM ((100.f - Jold100) / (100.f - redu), 0.f, 1.f);
                                        Jj = 0.3f * reduc * (Jj - Jold) + Jold; //reduct sensibility in highlights
                                    }
                                }
                            } else if (Jj > 10.f) {
                                Jj = 0.8f * (Jj - Jold) + Jold;
                            } else if (Jj >= 0.f) {
                                Jj = 0.90f * (Jj - Jold) + Jold;    // not zero ==>artifacts
                            }

                            Jpro = (float) (Jj / 327.68f);

                            if (Jpro < 1.f) {
                                Jpro = 1.f;
                            }
                        } else if (curveMode == ColorAppearanceParams::TC_MODE_BRIGHT) {
                            //attention! Brightness curves are open - unlike Lightness or Lab or RGB==> rendering  and algoritms will be different
                            float coef = ((aw + 4.f) * (4.f / c)) / 100.f;
                            float Qanc = Qpro;
                            float Qq = (float) Qpro * 327.68f * (1.f / coef);
                            float Qold100 = (float) Qpro / coef;

                            float Qold = Qq;
                            float redu = 20.f;
                            float reduc = 1.f;

                            const Brightcurve& userColCurveB1 = static_cast<const Brightcurve&> (customColCurve1);
                            userColCurveB1.Apply (Qq);

                            if (Qq > Qold) {
                                if (Qq < 65535.f)  {
                                    if (Qold < 327.68f * redu) {
                                        Qq = 0.25f * (Qq - Qold) + Qold;    //divide sensibility
                                    } else            {
                                        reduc = LIM ((100.f - Qold100) / (100.f - redu), 0.f, 1.f);
                                        Qq = 0.25f * reduc * (Qq - Qold) + Qold; //reduct sensibility in highlights
                                    }
                                }
                            } else if (Qq > 10.f) {
                                Qq = 0.5f * (Qq - Qold) + Qold;
                            } else if (Qq >= 0.f) {
                                Qq = 0.7f * (Qq - Qold) + Qold;    // not zero ==>artifacts
                            }

                            Qpro = Qanc * (Qq / Qold);
                            //   Jpro = 100.f * (Qpro * Qpro) / ((4.0f / c) * (4.0f / c) * (aw + 4.0f) * (aw + 4.0f));
                            Jpro = Jpro * SQR (Qq / Qold);

                            if (Jpro < 1.f) {
                                Jpro = 1.f;
                            }
                        }
                    }

                    if (hasColCurve2) {//curve 2 with Lightness and Brightness
                        if (curveMode2 == ColorAppearanceParams::TC_MODE_LIGHT) {
                            float Jj = (float) Jpro * 327.68f;
                            float Jold = Jj;
                            float Jold100 = (float) Jpro;
                            float redu = 25.f;
                            float reduc = 1.f;
                            const Lightcurve& userColCurveJ2 = static_cast<const Lightcurve&> (customColCurve2);
                            userColCurveJ2.Apply (Jj);

                            if (Jj > Jold) {
                                if (Jj < 65535.f)  {
                                    if (Jold < 327.68f * redu) {
                                        Jj = 0.3f * (Jj - Jold) + Jold;    //divide sensibility
                                    } else        {
                                        reduc = LIM ((100.f - Jold100) / (100.f - redu), 0.f, 1.f);
                                        Jj = 0.3f * reduc * (Jj - Jold) + Jold; //reduct sensibility in highlights
                                    }
                                }
                            } else if (Jj > 10.f) {
                                if (!t1L) {
                                    Jj = 0.8f * (Jj - Jold) + Jold;
                                } else {
                                    Jj = 0.4f * (Jj - Jold) + Jold;
                                }
                            } else if (Jj >= 0.f) {
                                if (!t1L) {
                                    Jj = 0.90f * (Jj - Jold) + Jold;    // not zero ==>artifacts
                                } else {
                                    Jj = 0.5f * (Jj - Jold) + Jold;
                                }
                            }

                            Jpro = (float) (Jj / 327.68f);

                            if (Jpro < 1.f) {
                                Jpro = 1.f;
                            }

                        } else if (curveMode2 == ColorAppearanceParams::TC_MODE_BRIGHT) { //
                            float Qanc = Qpro;

                            float coef = ((aw + 4.f) * (4.f / c)) / 100.f;
                            float Qq = (float) Qpro * 327.68f * (1.f / coef);
                            float Qold100 = (float) Qpro / coef;

                            float Qold = Qq;
                            float redu = 20.f;
                            float reduc = 1.f;

                            const Brightcurve& userColCurveB2 = static_cast<const Brightcurve&> (customColCurve2);
                            userColCurveB2.Apply (Qq);

                            if (Qq > Qold) {
                                if (Qq < 65535.f)  {
                                    if (Qold < 327.68f * redu) {
                                        Qq = 0.25f * (Qq - Qold) + Qold;    //divide sensibility
                                    } else            {
                                        reduc = LIM ((100.f - Qold100) / (100.f - redu), 0.f, 1.f);
                                        Qq = 0.25f * reduc * (Qq - Qold) + Qold; //reduct sensibility in highlights
                                    }
                                }
                            } else if (Qq > 10.f) {
                                Qq = 0.5f * (Qq - Qold) + Qold;
                            } else if (Qq >= 0.f) {
                                Qq = 0.7f * (Qq - Qold) + Qold;    // not zero ==>artifacts
                            }

                            //  Qpro = (float) (Qq * (coef) / 327.68f);
                            Qpro = Qanc * (Qq / Qold);
                            Jpro = Jpro * SQR (Qq / Qold);

                            // Jpro = 100.f * (Qpro * Qpro) / ((4.0f / c) * (4.0f / c) * (aw + 4.0f) * (aw + 4.0f));

                            if (t1L) { //to workaround the problem if we modify curve1-lightnees after curve2 brightness(the cat that bites its own tail!) in fact it's another type of curve only for this case
                                coef = 2.f; //adapt Q to J approximation
                                Qq = (float) Qpro * coef;
                                Qold = Qq;
                                const Lightcurve& userColCurveJ1 = static_cast<const Lightcurve&> (customColCurve1);
                                userColCurveJ1.Apply (Qq);
                                Qq = 0.05f * (Qq - Qold) + Qold; //approximative adaptation
                                Qpro = (float) (Qq / coef);
                                Jpro = 100.f * (Qpro * Qpro) / ((4.0f / c) * (4.0f / c) * (aw + 4.0f) * (aw + 4.0f));
                            }

                            if (Jpro < 1.f) {
                                Jpro = 1.f;
                            }
                        }
                    }

                    if (hasColCurve3) {//curve 3 with chroma saturation colorfullness
                        if (curveMode3 == ColorAppearanceParams::TC_MODE_CHROMA) {
                            float parsat = 0.8f; //0.68;
                            float coef = 327.68f / parsat;
                            float Cc = (float) Cpro * coef;
                            float Ccold = Cc;
                            const Chromacurve& userColCurve = static_cast<const Chromacurve&> (customColCurve3);
                            userColCurve.Apply (Cc);
                            float dred = 55.f;
                            float protect_red = 30.0f;
                            int sk = 1;
                            float ko = 1.f / coef;
                            Color::skinredfloat (Jpro, hpro, Cc, Ccold, dred, protect_red, sk, rstprotection, ko, Cpro);
                            /*
                                                        if(Jpro < 1.f && Cpro > 12.f) {
                                                            Cpro = 12.f;    //reduce artifacts by "pseudo gamut control CIECAM"
                                                        } else if(Jpro < 2.f && Cpro > 15.f) {
                                                            Cpro = 15.f;
                                                        } else if(Jpro < 4.f && Cpro > 30.f) {
                                                            Cpro = 30.f;
                                                        } else if(Jpro < 7.f && Cpro > 50.f) {
                                                            Cpro = 50.f;
                                                        }
                            */
                        } else if (curveMode3 == ColorAppearanceParams::TC_MODE_SATUR) { //
                            float parsat = 0.8f; //0.6
                            float coef = 327.68f / parsat;
                            float Ss = (float) spro * coef;
                            float Sold = Ss;
                            const Saturcurve& userColCurve = static_cast<const Saturcurve&> (customColCurve3);
                            userColCurve.Apply (Ss);
                            Ss = 0.6f * (Ss - Sold) + Sold; //divide sensibility saturation
                            float dred = 100.f; // in C mode
                            float protect_red = 80.0f; // in C mode
                            dred = 100.0f * sqrtf ((dred * coe) / Qpro);
                            protect_red = 100.0f * sqrtf ((protect_red * coe) / Qpro);
                            int sk = 0;
                            float ko = 1.f / coef;
                            Color::skinredfloat (Jpro, hpro, Ss, Sold, dred, protect_red, sk, rstprotection, ko, spro);
                            Qpro = ( 4.0f / c ) * sqrtf ( Jpro / 100.0f ) * ( aw + 4.0f ) ;
                            Cpro = (spro * spro * Qpro) / (10000.0f);
                        } else if (curveMode3 == ColorAppearanceParams::TC_MODE_COLORF) { //
                            float parsat = 0.8f; //0.68;
                            float coef = 327.68f / parsat;
                            float Mm = (float) Mpro * coef;
                            float Mold = Mm;
                            const Colorfcurve& userColCurve = static_cast<const Colorfcurve&> (customColCurve3);
                            userColCurve.Apply (Mm);
                            float dred = 100.f; //in C mode
                            float protect_red = 80.0f; // in C mode
                            dred *= coe; //in M mode
                            protect_red *= coe;
                            int sk = 0;
                            float ko = 1.f / coef;
                            Color::skinredfloat (Jpro, hpro, Mm, Mold, dred, protect_red, sk, rstprotection, ko, Mpro);
                            /*
                                                        if(Jpro < 1.f && Mpro > 12.f * coe) {
                                                            Mpro = 12.f * coe;    //reduce artifacts by "pseudo gamut control CIECAM"
                                                        } else if(Jpro < 2.f && Mpro > 15.f * coe) {
                                                            Mpro = 15.f * coe;
                                                        } else if(Jpro < 4.f && Mpro > 30.f * coe) {
                                                            Mpro = 30.f * coe;
                                                        } else if(Jpro < 7.f && Mpro > 50.f * coe) {
                                                            Mpro = 50.f * coe;
                                                        }
                            */
                            Cpro = Mpro / coe;
                        }
                    }

                    //to retrieve the correct values of variables


                    //retrieve values C,J...s
                    C = Cpro;
                    J = Jpro;
                    Q = Qpro;
                    M = Mpro;
                    h = hpro;
                    s = spro;

                    if (params->colorappearance.tonecie  || settings->autocielab) { //use pointer for tonemapping with CIECAM and also sharpening , defringe, contrast detail
                        ncie->Q_p[i][j] = (float)Q + epsil; //epsil to avoid Q=0
                        ncie->M_p[i][j] = (float)M + epsil;
                        ncie->J_p[i][j] = (float)J + epsil;
                        ncie->h_p[i][j] = (float)h;
                        ncie->C_p[i][j] = (float)C + epsil;
                        ncie->sh_p[i][j] = (float) 3276.8f * (sqrtf ( J ) ) ;

                        if (epdEnabled) {
                            if (ncie->Q_p[i][j] < minQThr) {
                                minQThr = ncie->Q_p[i][j];    //minima
                            }

                            if (ncie->Q_p[i][j] > maxQThr) {
                                maxQThr = ncie->Q_p[i][j];    //maxima
                            }
                        }
                    }

                    if (!params->colorappearance.tonecie  || !settings->autocielab || !epdEnabled) {

                        if (ciedata) { //only with improccoordinator
                            // Data for J Q M s and C histograms
                            int posl, posc;
                            float brli = 327.f;
                            float chsacol = 327.f;
                            float libr;
                            float colch;

                            //update histogram
                            if (curveMode == ColorAppearanceParams::TC_MODE_BRIGHT) {
                                brli = 70.0f;
                                libr = Q;     //40.0 to 100.0 approximative factor for Q  - 327 for J
                            } else { /*if(curveMode == ColorAppearanceParams::TC_MODE_LIGHT)*/
                                brli = 327.f;
                                libr = J;    //327 for J
                            }

                            posl = (int) (libr * brli);
                            hist16JCAM[posl]++;

                            if (curveMode3 == ColorAppearanceParams::TC_MODE_CHROMA) {
                                chsacol = 400.f;//327
                                colch = C;    //450.0 approximative factor for s    320 for M
                            } else if (curveMode3 == ColorAppearanceParams::TC_MODE_SATUR) {
                                chsacol = 450.0f;
                                colch = s;
                            } else { /*if(curveMode3 == ColorAppearanceParams::TC_MODE_COLORF)*/
                                chsacol = 400.0f;//327
                                colch = M;
                            }

                            posc = (int) (colch * chsacol);
                            hist16_CCAM[posc]++;

                        }

                        if (LabPassOne) {
#ifdef __SSE2__
                            // write to line buffers
                            Jbuffer[j] = J;
                            Cbuffer[j] = C;
                            hbuffer[j] = h;
#else
                            float xx, yy, zz;
                            //process normal==> viewing

                            Ciecam02::jch2xyz_ciecam02float ( xx, yy, zz,
                                                              J,  C, h,
                                                              xw2, yw2,  zw2,
                                                              f2,  c2, nc2, gamu, pow1n, nbbj, ncbj, flj, czj, dj, awj);
                            float x, y, z;
                            x = xx * 655.35f;
                            y = yy * 655.35f;
                            z = zz * 655.35f;
                            float Ll, aa, bb;
                            //convert xyz=>lab
                            Color::XYZ2Lab (x,  y,  z, Ll, aa, bb);

                            // gamut control in Lab mode; I must study how to do with cIECAM only
                            if (gamu == 1) {
                                float HH, Lprov1, Chprov1;
                                Lprov1 = Ll / 327.68f;
                                Chprov1 = sqrtf (SQR (aa) + SQR (bb)) / 327.68f;
                                HH = xatan2f (bb, aa);
                                float2  sincosval;

                                if (Chprov1 == 0.0f) {
                                    sincosval.y = 1.f;
                                    sincosval.x = 0.0f;
                                } else {
                                    sincosval.y = aa / (Chprov1 * 327.68f);
                                    sincosval.x = bb / (Chprov1 * 327.68f);
                                }


#ifdef _DEBUG
                                bool neg = false;
                                bool more_rgb = false;
                                //gamut control : Lab values are in gamut
                                Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                                //gamut control : Lab values are in gamut
                                Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f);
#endif

                                lab->L[i][j] = Lprov1 * 327.68f;
                                lab->a[i][j] = 327.68f * Chprov1 * sincosval.y;
                                lab->b[i][j] = 327.68f * Chprov1 * sincosval.x;

                            } else {
                                lab->L[i][j] = Ll;
                                lab->a[i][j] = aa;
                                lab->b[i][j] = bb;
                            }

#endif
                        }
                    }
                }

#ifdef __SSE2__
                // process line buffers
                float *xbuffer = Qbuffer;
                float *ybuffer = Mbuffer;
                float *zbuffer = sbuffer;

                for (k = 0; k < bufferLength; k += 4) {
                    Ciecam02::jch2xyz_ciecam02float ( x, y, z,
                                                      LVF (Jbuffer[k]), LVF (Cbuffer[k]), LVF (hbuffer[k]),
                                                      F2V (xw2), F2V (yw2), F2V (zw2),
                                                      F2V (f2),  F2V (nc2), F2V (pow1n), F2V (nbbj), F2V (ncbj), F2V (flj), F2V (dj), F2V (awj), F2V (reccmcz));
                    STVF (xbuffer[k], x * c655d35);
                    STVF (ybuffer[k], y * c655d35);
                    STVF (zbuffer[k], z * c655d35);
                }

                // XYZ2Lab uses a lookup table. The function behind that lut is a cube root.
                // SSE can't beat the speed of that lut, so it doesn't make sense to use SSE
                for (int j = 0; j < width; j++) {
                    float Ll, aa, bb;
                    //convert xyz=>lab
                    Color::XYZ2Lab (xbuffer[j], ybuffer[j], zbuffer[j], Ll, aa, bb);

                    // gamut control in Lab mode; I must study how to do with cIECAM only
                    if (gamu == 1) {
                        float Lprov1, Chprov1;
                        Lprov1 = Ll / 327.68f;
                        Chprov1 = sqrtf (SQR (aa) + SQR (bb)) / 327.68f;
                        float2  sincosval;

                        if (Chprov1 == 0.0f) {
                            sincosval.y = 1.f;
                            sincosval.x = 0.0f;
                        } else {
                            sincosval.y = aa / (Chprov1 * 327.68f);
                            sincosval.x = bb / (Chprov1 * 327.68f);
                        }

#ifdef _DEBUG
                        bool neg = false;
                        bool more_rgb = false;
                        //gamut control : Lab values are in gamut
                        Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                        //gamut control : Lab values are in gamut
                        Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f);
#endif
                        lab->L[i][j] = Lprov1 * 327.68f;
                        lab->a[i][j] = 327.68f * Chprov1 * sincosval.y;
                        lab->b[i][j] = 327.68f * Chprov1 * sincosval.x;
                    } else {
                        lab->L[i][j] = Ll;
                        lab->a[i][j] = aa;
                        lab->b[i][j] = bb;
                    }
                }

#endif
            }

            #pragma omp critical
            {
                if (minQThr < minQ) {
                    minQ = minQThr;
                }

                if (maxQThr > maxQ) {
                    maxQ = maxQThr;
                }
            }
        }

        // End of parallelization
        if (!params->colorappearance.tonecie || !settings->autocielab) { //normal

            if (ciedata) {
                //update histogram J
                hist16JCAM.compressTo (histLCAM);
                //update histogram C
                hist16_CCAM.compressTo (histCCAM);
            }
        }

#ifdef _DEBUG

        if (settings->verbose) {
            t2e.set();
            printf ("CIECAM02 performed in %d usec:\n", t2e.etime (t1e));
            //  printf("minc=%f maxc=%f minj=%f maxj=%f\n",minc,maxc,minj,maxj);
        }

#endif

        if (settings->autocielab) {
            if ((params->colorappearance.tonecie && (epdEnabled)) || (params->sharpening.enabled && settings->autocielab && execsharp)
                    || (params->dirpyrequalizer.enabled && settings->autocielab) || (params->defringe.enabled && settings->autocielab)  || (params->sharpenMicro.enabled && settings->autocielab)
                    || (params->impulseDenoise.enabled && settings->autocielab) ||  (params->colorappearance.badpixsl > 0 && settings->autocielab)) {



//all this treatments reduce artefacts, but can leed to slighty  different results

                if (params->defringe.enabled)
                    if (execsharp) {
                        lab->deleteLab();
                        ImProcFunctions::defringecam (ncie);//defringe adapted to CIECAM
                        lab->reallocLab();
                    }

//if(params->dirpyrequalizer.enabled) if(execsharp) {
                if (params->dirpyrequalizer.enabled)  {
                    if (params->dirpyrequalizer.gamutlab  /*&& execsharp*/) { //remove artifacts by gaussian blur - skin control
                        float b_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[0]) / 100.0f;
                        float t_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[1]) / 100.0f;
                        float b_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[2]) / 100.0f;
                        float t_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[3]) / 100.0f;
                        float artifact = (float) settings->artifact_cbdl;

                        if (artifact > 6.f) {
                            artifact = 6.f;
                        }

                        if (artifact < 0.f) {
                            artifact = 1.f;
                        }

                        int hotbad = 0;
                        float chrom = 50.f;
                        lab->deleteLab();
                        ImProcFunctions::badpixcam (ncie, artifact, 5, 2, b_l, t_l, t_r, b_r, params->dirpyrequalizer.skinprotect, chrom, hotbad);  //enabled remove artifacts for cbDL
                        lab->reallocLab();
                    }
                }

//if(params->colorappearance.badpixsl > 0) { int mode=params->colorappearance.badpixsl;
                if (params->colorappearance.badpixsl > 0) if (execsharp) {
                        int mode = params->colorappearance.badpixsl;
                        lab->deleteLab();
                        ImProcFunctions::badpixcam (ncie, 3.0, 10, mode, 0, 0, 0, 0, 0, 0, 1);//for bad pixels CIECAM
                        lab->reallocLab();
                    }

                if (params->impulseDenoise.enabled) if (execsharp) {
                        float **buffers[3];
                        buffers[0] = lab->L;
                        buffers[1] = lab->a;
                        buffers[2] = lab->b;
                        ImProcFunctions::impulsedenoisecam (ncie, buffers); //impulse adapted to CIECAM
                    }

                if (params->sharpenMicro.enabled)if (execsharp) {
                        ImProcFunctions::MLmicrocontrastcam (ncie);
                    }

                if (params->sharpening.enabled)
                    if (execsharp) {
                        float **buffer = lab->L; // We can use the L-buffer from lab as buffer to save some memory
                        ImProcFunctions::sharpeningcam (ncie, buffer); // sharpening adapted to CIECAM
                    }

//if(params->dirpyrequalizer.enabled) if(execsharp) {
                if (params->dirpyrequalizer.enabled /*&& execsharp*/)  {
//  if(params->dirpyrequalizer.algo=="FI") choice=0;
//  else if(params->dirpyrequalizer.algo=="LA") choice=1;

                    if (rtt == 1) {
                        float b_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[0]) / 100.0f;
                        float t_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[1]) / 100.0f;
                        float b_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[2]) / 100.0f;
                        float t_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[3]) / 100.0f;
                        int choice = 0; // I have not suppress this statement in case of !! always to 0
                        lab->deleteLab();
                        dirpyr_equalizercam (ncie, ncie->sh_p, ncie->sh_p, ncie->W, ncie->H, ncie->h_p, ncie->C_p, params->dirpyrequalizer.mult, params->dirpyrequalizer.threshold, params->dirpyrequalizer.skinprotect,  true, params->dirpyrequalizer.gamutlab, b_l, t_l, t_r, b_r, choice, scale); //contrast by detail adapted to CIECAM
                        lab->reallocLab();
                    }

                    /*
                    if(params->colorappearance.badpixsl > 0) if(execsharp){ int mode=params->colorappearance.badpixsl;
                    printf("BADPIX");
                                                                ImProcFunctions::badpixcam (ncie, 8.0, 10, mode);//for bad pixels
                                                            }
                                                            */
                }

                const float Qredi = ( 4.0f / c_)  * ( a_w + 4.0f );
                const float co_e = (pow_F (f_l, 0.25f));


#ifndef _DEBUG
                #pragma omp parallel
#endif
                {
#ifndef _DEBUG
                    #pragma omp for schedule(dynamic, 10)
#endif

                    for (int i = 0; i < height; i++) // update CieImages with new values after sharpening, defringe, contrast by detail level
                        for (int j = 0; j < width; j++) {
                            float interm = fabsf (ncie->sh_p[i][j] / (32768.f));
                            ncie->J_p[i][j] = 100.0f * SQR (interm);
                            ncie->Q_p[i][j] = interm * Qredi;
                            ncie->M_p[i][j] = ncie->C_p[i][j] * co_e;
                        }
                }
            }
        }

        if ((params->colorappearance.tonecie && (epdEnabled)) || (params->sharpening.enabled && settings->autocielab && execsharp)
                || (params->dirpyrequalizer.enabled && settings->autocielab) || (params->defringe.enabled && settings->autocielab)  || (params->sharpenMicro.enabled && settings->autocielab)
                || (params->impulseDenoise.enabled && settings->autocielab) ||  (params->colorappearance.badpixsl > 0 && settings->autocielab)) {

            ciedata = (params->colorappearance.datacie && pW != 1);

            if (epdEnabled  && params->colorappearance.tonecie && algepd) {
                lab->deleteLab();
                ImProcFunctions::EPDToneMapCIE (ncie, a_w, c_, w_h, width, height, begh, endh, minQ, maxQ, Iterates, scale );
                lab->reallocLab();
            }

            //EPDToneMapCIE adated to CIECAM


            constexpr float eps = 0.0001f;
            const float co_e = (pow_F (f_l, 0.25f)) + eps;

#ifndef _DEBUG
            #pragma omp parallel
#endif
            {
#ifdef __SSE2__
                // one line buffer per channel
                float Jbuffer[bufferLength] ALIGNED16;
                float Cbuffer[bufferLength] ALIGNED16;
                float hbuffer[bufferLength] ALIGNED16;
                float *xbuffer = Jbuffer; // we can use one of the above buffers
                float *ybuffer = Cbuffer; //             "
                float *zbuffer = hbuffer; //             "
#endif

#ifndef _DEBUG
                #pragma omp for schedule(dynamic, 10)
#endif

                for (int i = 0; i < height; i++) { // update CIECAM with new values after tone-mapping
                    for (int j = 0; j < width; j++) {

                        //  if(epdEnabled) ncie->J_p[i][j]=(100.0f* ncie->Q_p[i][j]*ncie->Q_p[i][j])/(w_h*w_h);
                        if (epdEnabled) {
                            ncie->J_p[i][j] = (100.0f * ncie->Q_p[i][j] * ncie->Q_p[i][j]) / SQR ((4.f / c) * (aw + 4.f));
                        }

                        const float ncie_C_p = (ncie->M_p[i][j]) / co_e;

                        //show histogram in CIECAM mode (Q,J, M,s,C)
                        if (ciedata) {
                            // Data for J Q M s and C histograms
                            int posl, posc;
                            float brli = 327.f;
                            float chsacol = 327.f;
                            float libr;
                            float colch;

                            if (curveMode == ColorAppearanceParams::TC_MODE_BRIGHT) {
                                brli = 70.0f;
                                libr = ncie->Q_p[i][j];    //40.0 to 100.0 approximative factor for Q  - 327 for J
                            } else { /*if(curveMode == ColorAppearanceParams::TC_MODE_LIGHT)*/
                                brli = 327.f;
                                libr = ncie->J_p[i][j];    //327 for J
                            }

                            posl = (int) (libr * brli);
                            hist16JCAM[posl]++;

                            if (curveMode3 == ColorAppearanceParams::TC_MODE_CHROMA) {
                                chsacol = 400.f;//327.f;
                                colch = ncie_C_p;
                            } else if (curveMode3 == ColorAppearanceParams::TC_MODE_SATUR) {
                                chsacol = 450.0f;
                                colch = 100.f * sqrtf (ncie_C_p / ncie->Q_p[i][j]);
                            } else { /*if(curveMode3 == ColorAppearanceParams::TC_MODE_COLORF)*/
                                chsacol = 400.f;//327.0f;
                                colch = ncie->M_p[i][j];
                            }

                            posc = (int) (colch * chsacol);
                            hist16_CCAM[posc]++;
                        }

                        //end histograms

#ifdef __SSE2__
                        Jbuffer[j] = ncie->J_p[i][j];
                        Cbuffer[j] = ncie_C_p;
                        hbuffer[j] = ncie->h_p[i][j];
#else
                        float xx, yy, zz;
                        Ciecam02::jch2xyz_ciecam02float ( xx, yy, zz,
                                                          ncie->J_p[i][j],  ncie_C_p, ncie->h_p[i][j],
                                                          xw2, yw2,  zw2,
                                                          f2,  c2, nc2, gamu, pow1n, nbbj, ncbj, flj, czj, dj, awj);
                        float x = (float)xx * 655.35f;
                        float y = (float)yy * 655.35f;
                        float z = (float)zz * 655.35f;
                        float Ll, aa, bb;
                        //convert xyz=>lab
                        Color::XYZ2Lab (x,  y,  z, Ll, aa, bb);

                        if (gamu == 1) {
                            float Lprov1, Chprov1;
                            Lprov1 = Ll / 327.68f;
                            Chprov1 = sqrtf (SQR (aa) + SQR (bb)) / 327.68f;
                            float2  sincosval;

                            if (Chprov1 == 0.0f) {
                                sincosval.y = 1.f;
                                sincosval.x = 0.0f;
                            } else {
                                sincosval.y = aa / (Chprov1 * 327.68f);
                                sincosval.x = bb / (Chprov1 * 327.68f);
                            }


#ifdef _DEBUG
                            bool neg = false;
                            bool more_rgb = false;
                            //gamut control : Lab values are in gamut
                            Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                            //gamut control : Lab values are in gamut
                            Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f);
#endif

                            lab->L[i][j] = Lprov1 * 327.68f;
                            lab->a[i][j] = 327.68f * Chprov1 * sincosval.y;
                            lab->b[i][j] = 327.68f * Chprov1 * sincosval.x;
                        } else {
                            lab->L[i][j] = Ll;
                            lab->a[i][j] = aa;
                            lab->b[i][j] = bb;
                        }

#endif
                    }

#ifdef __SSE2__
                    // process line buffers
                    int k;
                    vfloat x, y, z;
                    vfloat c655d35 = F2V (655.35f);

                    for (k = 0; k < bufferLength; k += 4) {
                        Ciecam02::jch2xyz_ciecam02float ( x, y, z,
                                                          LVF (Jbuffer[k]), LVF (Cbuffer[k]), LVF (hbuffer[k]),
                                                          F2V (xw2), F2V (yw2), F2V (zw2),
                                                          F2V (f2), F2V (nc2), F2V (pow1n), F2V (nbbj), F2V (ncbj), F2V (flj), F2V (dj), F2V (awj), F2V (reccmcz));
                        x *= c655d35;
                        y *= c655d35;
                        z *= c655d35;
                        STVF (xbuffer[k], x);
                        STVF (ybuffer[k], y);
                        STVF (zbuffer[k], z);
                    }

                    // XYZ2Lab uses a lookup table. The function behind that lut is a cube root.
                    // SSE can't beat the speed of that lut, so it doesn't make sense to use SSE
                    for (int j = 0; j < width; j++) {
                        float Ll, aa, bb;
                        //convert xyz=>lab
                        Color::XYZ2Lab (xbuffer[j], ybuffer[j], zbuffer[j], Ll, aa, bb);

                        if (gamu == 1) {
                            float Lprov1, Chprov1;
                            Lprov1 = Ll / 327.68f;
                            Chprov1 = sqrtf (SQR (aa) + SQR (bb)) / 327.68f;
                            float2  sincosval;

                            if (Chprov1 == 0.0f) {
                                sincosval.y = 1.f;
                                sincosval.x = 0.0f;
                            } else {
                                sincosval.y = aa / (Chprov1 * 327.68f);
                                sincosval.x = bb / (Chprov1 * 327.68f);
                            }

#ifdef _DEBUG
                            bool neg = false;
                            bool more_rgb = false;
                            //gamut control : Lab values are in gamut
                            Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                            //gamut control : Lab values are in gamut
                            Color::gamutLchonly (sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f);
#endif
                            lab->L[i][j] = Lprov1 * 327.68f;
                            lab->a[i][j] = 327.68f * Chprov1 * sincosval.y;
                            lab->b[i][j] = 327.68f * Chprov1 * sincosval.x;
                        } else {
                            lab->L[i][j] = Ll;
                            lab->a[i][j] = aa;
                            lab->b[i][j] = bb;
                        }

                    }

#endif // __SSE2__
                }

            } //end parallelization

            //show CIECAM histograms
            if (ciedata) {
                //update histogram J and Q
                //update histogram J
                hist16JCAM.compressTo (histLCAM);

                //update color histogram M,s,C
                hist16_CCAM.compressTo (histCCAM);
            }
        }
    }
}
//end CIECAM

void ImProcFunctions::moyeqt (Imagefloat* working, float &moyS, float &eqty)
{
    BENCHFUN

    int tHh = working->getHeight();
    int tWw = working->getWidth();
    double moy = 0.0;
    double sqrs = 0.0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:moy,sqrs) schedule(dynamic,16)
#endif

    for (int i = 0; i < tHh; i++) {
        for (int j = 0; j < tWw; j++) {
            float s = Color::rgb2s (CLIP (working->r (i, j)), CLIP (working->g (i, j)), CLIP (working->b (i, j)));
            moy += s;
            sqrs += SQR (s);
        }
    }

    moy /= (tHh * tWw);
    sqrs /= (tHh * tWw);
    eqty = sqrt (sqrs - SQR (moy));
    moyS = moy;
}

static inline void
filmlike_clip_rgb_tone (float *r, float *g, float *b, const float L)
{
    float r_ = *r > L ? L : *r;
    float b_ = *b > L ? L : *b;
    float g_ = b_ + ((r_ - b_) * (*g - *b) / (*r - *b));
    *r = r_;
    *g = g_;
    *b = b_;
}

static void
filmlike_clip (float *r, float *g, float *b)
{
    // This is Adobe's hue-stable film-like curve with a diagonal, ie only used for clipping. Can probably be further optimized.
    const float L = 65535.0;

    if (*r >= *g) {
        if (*g > *b) {         // Case 1: r >= g >  b
            filmlike_clip_rgb_tone (r, g, b, L);
        } else if (*b > *r) {  // Case 2: b >  r >= g
            filmlike_clip_rgb_tone (b, r, g, L);
        } else if (*b > *g) {  // Case 3: r >= b >  g
            filmlike_clip_rgb_tone (r, b, g, L);
        } else {               // Case 4: r >= g == b
            *r = *r > L ? L : *r;
            *g = *g > L ? L : *g;
            *b = *g;
        }
    } else {
        if (*r >= *b) {        // Case 5: g >  r >= b
            filmlike_clip_rgb_tone (g, r, b, L);
        } else if (*b > *g) {  // Case 6: b >  g >  r
            filmlike_clip_rgb_tone (b, g, r, L);
        } else {               // Case 7: g >= b >  r
            filmlike_clip_rgb_tone (g, b, r, L);
        }
    }
}

void ImProcFunctions::rgbProc (Imagefloat* working, LabImage* lab, PipetteBuffer *pipetteBuffer, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                               SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve, float satLimit, float satLimitOpacity, const ColorGradientCurve & ctColorCurve, const OpacityCurve & ctOpacityCurve, bool opautili,  LUTf & clToningcurve, LUTf & cl2Toningcurve,
                               const ToneCurve & customToneCurve1, const ToneCurve & customToneCurve2, const ToneCurve & customToneCurvebw1, const ToneCurve & customToneCurvebw2, double &rrm, double &ggm, double &bbm, float &autor, float &autog, float &autob, DCPProfile *dcpProf, const DCPProfile::ApplyState &asIn, LUTu &histToneCurve )
{
    rgbProc (working, lab, pipetteBuffer, hltonecurve, shtonecurve, tonecurve, shmap, sat, rCurve, gCurve, bCurve, satLimit, satLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2,  customToneCurvebw1, customToneCurvebw2, rrm, ggm, bbm, autor, autog, autob, params->toneCurve.expcomp, params->toneCurve.hlcompr, params->toneCurve.hlcomprthresh, dcpProf, asIn, histToneCurve);
}

// Process RGB image and convert to LAB space
void ImProcFunctions::rgbProc (Imagefloat* working, LabImage* lab, PipetteBuffer *pipetteBuffer, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                               SHMap* shmap, int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve, float satLimit, float satLimitOpacity, const ColorGradientCurve & ctColorCurve, const OpacityCurve & ctOpacityCurve, bool opautili, LUTf & clToningcurve, LUTf & cl2Toningcurve,
                               const ToneCurve & customToneCurve1, const ToneCurve & customToneCurve2,  const ToneCurve & customToneCurvebw1, const ToneCurve & customToneCurvebw2, double &rrm, double &ggm, double &bbm, float &autor, float &autog, float &autob, double expcomp, int hlcompr, int hlcomprthresh, DCPProfile *dcpProf, const DCPProfile::ApplyState &asIn, LUTu &histToneCurve )
{
    BENCHFUN
    Imagefloat *tmpImage = nullptr;

    Imagefloat* editImgFloat = nullptr;
    PlanarWhateverData<float>* editWhatever = nullptr;
    EditUniqueID editID = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if (editID != EUID_None) {
        switch  (pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType()) {
            case (BT_IMAGEFLOAT):
                editImgFloat = pipetteBuffer->getImgFloatBuffer();
                break;

            case (BT_LABIMAGE):
                break;

            case (BT_SINGLEPLANE_FLOAT):
                editWhatever = pipetteBuffer->getSinglePlaneBuffer();
                break;
        }
    }

    int h_th = 0, s_th = 0;

    if (shmap) {
        h_th = shmap->max_f - params->sh.htonalwidth * (shmap->max_f - shmap->avg) / 100;
        s_th = params->sh.stonalwidth * (shmap->avg - shmap->min_f) / 100;
    }

    bool processSH  = params->sh.enabled && shmap && (params->sh.highlights > 0 || params->sh.shadows > 0);
    bool processLCE = params->sh.enabled && shmap && params->sh.localcontrast > 0;
    double lceamount = params->sh.localcontrast / 200.0;

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params->icm.working);
    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.working);

    float toxyz[3][3] = {
        {
            static_cast<float> ( wprof[0][0] / Color::D50x),
            static_cast<float> ( wprof[0][1] / Color::D50x),
            static_cast<float> ( wprof[0][2] / Color::D50x)
        }, {
            static_cast<float> ( wprof[1][0]),
            static_cast<float> ( wprof[1][1]),
            static_cast<float> ( wprof[1][2])
        }, {
            static_cast<float> ( wprof[2][0] / Color::D50z),
            static_cast<float> ( wprof[2][1] / Color::D50z),
            static_cast<float> ( wprof[2][2] / Color::D50z)
        }
    };
    float maxFactorToxyz = max (toxyz[1][0], toxyz[1][1], toxyz[1][2]);
    float equalR = maxFactorToxyz / toxyz[1][0];
    float equalG = maxFactorToxyz / toxyz[1][1];
    float equalB = maxFactorToxyz / toxyz[1][2];

    //inverse matrix user select
    double wip[3][3] = {
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };

    double wp[3][3] = {
        {wprof[0][0], wprof[0][1], wprof[0][2]},
        {wprof[1][0], wprof[1][1], wprof[1][2]},
        {wprof[2][0], wprof[2][1], wprof[2][2]}
    };

    bool mixchannels = (params->chmixer.red[0] != 100 || params->chmixer.red[1] != 0     || params->chmixer.red[2] != 0   ||
                        params->chmixer.green[0] != 0 || params->chmixer.green[1] != 100 || params->chmixer.green[2] != 0 ||
                        params->chmixer.blue[0] != 0  || params->chmixer.blue[1] != 0    || params->chmixer.blue[2] != 100);

    FlatCurve* hCurve = nullptr;
    FlatCurve* sCurve = nullptr;
    FlatCurve* vCurve = nullptr;
    FlatCurve* bwlCurve = nullptr;

    FlatCurveType hCurveType = (FlatCurveType)params->hsvequalizer.hcurve.at (0);
    FlatCurveType sCurveType = (FlatCurveType)params->hsvequalizer.scurve.at (0);
    FlatCurveType vCurveType = (FlatCurveType)params->hsvequalizer.vcurve.at (0);
    FlatCurveType bwlCurveType = (FlatCurveType)params->blackwhite.luminanceCurve.at (0);
    bool hCurveEnabled = hCurveType > FCT_Linear;
    bool sCurveEnabled = sCurveType > FCT_Linear;
    bool vCurveEnabled = vCurveType > FCT_Linear;
    bool bwlCurveEnabled = bwlCurveType > FCT_Linear;

    // TODO: We should create a 'skip' value like for CurveFactory::complexsgnCurve (rtengine/curves.cc)
    if (hCurveEnabled) {
        hCurve = new FlatCurve (params->hsvequalizer.hcurve);

        if (hCurve->isIdentity()) {
            delete hCurve;
            hCurve = nullptr;
            hCurveEnabled = false;
        }
    }

    if (sCurveEnabled) {
        sCurve = new FlatCurve (params->hsvequalizer.scurve);

        if (sCurve->isIdentity()) {
            delete sCurve;
            sCurve = nullptr;
            sCurveEnabled = false;
        }
    }

    if (vCurveEnabled) {
        vCurve = new FlatCurve (params->hsvequalizer.vcurve);

        if (vCurve->isIdentity()) {
            delete vCurve;
            vCurve = nullptr;
            vCurveEnabled = false;
        }
    }

    if (bwlCurveEnabled) {
        bwlCurve = new FlatCurve (params->blackwhite.luminanceCurve);

        if (bwlCurve->isIdentity()) {
            delete bwlCurve;
            bwlCurve = nullptr;
            bwlCurveEnabled = false;
        }
    }

    std::shared_ptr<HaldCLUT> hald_clut;
    bool clutAndWorkingProfilesAreSame = false;
    TMatrix xyz2clut = {}, clut2xyz = {};
#ifdef __SSE2__
    vfloat v_work2xyz[3][3] ALIGNED16;
    vfloat v_xyz2clut[3][3] ALIGNED16;
    vfloat v_clut2xyz[3][3] ALIGNED16;
    vfloat v_xyz2work[3][3] ALIGNED16;
#endif

    if ( params->filmSimulation.enabled && !params->filmSimulation.clutFilename.empty() ) {
        hald_clut = CLUTStore::getInstance().getClut ( params->filmSimulation.clutFilename );

        if ( hald_clut ) {
            clutAndWorkingProfilesAreSame = hald_clut->getProfile() == params->icm.working;

            if ( !clutAndWorkingProfilesAreSame ) {
                xyz2clut = ICCStore::getInstance()->workingSpaceInverseMatrix ( hald_clut->getProfile() );
                clut2xyz = ICCStore::getInstance()->workingSpaceMatrix ( hald_clut->getProfile() );

#ifdef __SSE2__

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        v_work2xyz[i][j] = F2V (wprof[i][j]);
                        v_xyz2clut[i][j] = F2V (xyz2clut[i][j]);
                        v_xyz2work[i][j] = F2V (wiprof[i][j]);
                        v_clut2xyz[i][j] = F2V (clut2xyz[i][j]);
                    }
                }

#endif

            }
        }
    }

    const float film_simulation_strength = static_cast<float> (params->filmSimulation.strength) / 100.0f;

    const float exp_scale = pow (2.0, expcomp);
    const float comp = (max (0.0, expcomp) + 1.0) * hlcompr / 100.0;
    const float shoulder = ((65536.0 / max (1.0f, exp_scale)) * (hlcomprthresh / 200.0)) + 0.1;
    const float hlrange = 65536.0 - shoulder;
    const bool isProPhoto = (params->icm.working == "ProPhoto");
    // extracting datas from 'params' to avoid cache flush (to be confirmed)
    ToneCurveParams::eTCModeId curveMode = params->toneCurve.curveMode;
    ToneCurveParams::eTCModeId curveMode2 = params->toneCurve.curveMode2;
    bool highlight = params->toneCurve.hrenabled;//Get the value if "highlight reconstruction" is activated
    bool hasToneCurve1 = bool (customToneCurve1);
    bool hasToneCurve2 = bool (customToneCurve2);
    BlackWhiteParams::eTCModeId beforeCurveMode = params->blackwhite.beforeCurveMode;
    BlackWhiteParams::eTCModeId afterCurveMode = params->blackwhite.afterCurveMode;

    bool hasToneCurvebw1 = bool (customToneCurvebw1);
    bool hasToneCurvebw2 = bool (customToneCurvebw2);

    PerceptualToneCurveState ptc1ApplyState, ptc2ApplyState;

    if (hasToneCurve1 && curveMode == ToneCurveParams::TC_MODE_PERCEPTUAL) {
        const PerceptualToneCurve& userToneCurve = static_cast<const PerceptualToneCurve&> (customToneCurve1);
        userToneCurve.initApplyState (ptc1ApplyState, params->icm.working);
    }

    if (hasToneCurve2 && curveMode2 == ToneCurveParams::TC_MODE_PERCEPTUAL) {
        const PerceptualToneCurve& userToneCurve = static_cast<const PerceptualToneCurve&> (customToneCurve2);
        userToneCurve.initApplyState (ptc2ApplyState, params->icm.working);
    }

    bool hasColorToning = params->colorToning.enabled && bool (ctOpacityCurve) &&  bool (ctColorCurve);
    //  float satLimit = float(params->colorToning.satProtectionThreshold)/100.f*0.7f+0.3f;
    //  float satLimitOpacity = 1.f-(float(params->colorToning.saturatedOpacity)/100.f);
    float strProtect = (float (params->colorToning.strength) / 100.f);

    /*
    // Debug output - Color LUTf points
    if (ctColorCurve) {
        printf("\nColor curve:");
        for (size_t i=0; i<501; i++) {
            if (i==0 || i==250 || i==500)
                printf("\n(%.1f)[", float(i)/500.f);
            printf("%.3f ", ctColorCurve.lutHueCurve[float(i)]);
            if (i==0 || i==250 || i==500)
            printf("]\n");
        }
        printf("\n");
    }
    */

    /*
    // Debug output - Opacity LUTf points
    if (ctOpacityCurve) {
        printf("\nOpacity curve:");
        for (size_t i=0; i<501; i++) {
            if (i==0 || i==250 || i==500)
                printf("\n(%.1f)[", float(i)/500.f);
            printf("%.3f ", ctOpacityCurve.lutOpacityCurve[float(i)]);
            if (i==0 || i==250 || i==500)
            printf("]\n");
        }
        printf("\n");
    }
    */

    float RedLow = (100.f + float (params->colorToning.redlow)) / 100.f; //printf("Rel=%f\n",RedLow);
    float GreenLow = (100.f + float (params->colorToning.greenlow)) / 100.f; //printf("Gre=%f\n",GreenLow);
    float BlueLow = (100.f + float (params->colorToning.bluelow)) / 100.f; //printf("Blu=%f\n",BlueLow);
    float RedMed = (100.f + float (params->colorToning.redmed)) / 100.f;
    float GreenMed = (100.f + float (params->colorToning.greenmed)) / 100.f;
    float BlueMed = (100.f + float (params->colorToning.bluemed)) / 100.f;
    float RedHigh = (100.f + float (params->colorToning.redhigh)) / 100.f; //printf("RedH=%f\n",RedHigh);
    float GreenHigh = (100.f + float (params->colorToning.greenhigh)) / 100.f;
    float BlueHigh = (100.f + float (params->colorToning.bluehigh)) / 100.f;
    float SatLow = float (params->colorToning.shadowsColSat.value[0]) / 100.f;
    float SatHigh = float (params->colorToning.hlColSat.value[0]) / 100.f;

    float Balan = float (params->colorToning.balance);

    float chMixRR = float (params->chmixer.red[0]);
    float chMixRG = float (params->chmixer.red[1]);
    float chMixRB = float (params->chmixer.red[2]);
    float chMixGR = float (params->chmixer.green[0]);
    float chMixGG = float (params->chmixer.green[1]);
    float chMixGB = float (params->chmixer.green[2]);
    float chMixBR = float (params->chmixer.blue[0]);
    float chMixBG = float (params->chmixer.blue[1]);
    float chMixBB = float (params->chmixer.blue[2]);

    int shHighlights = params->sh.highlights;
    int shShadows = params->sh.shadows;
    bool blackwhite = params->blackwhite.enabled;
    bool complem = params->blackwhite.enabledcc;
    float bwr = float (params->blackwhite.mixerRed);
    float bwg = float (params->blackwhite.mixerGreen);
    float bwb = float (params->blackwhite.mixerBlue);
    float bwrgam = float (params->blackwhite.gammaRed);
    float bwggam = float (params->blackwhite.gammaGreen);
    float bwbgam = float (params->blackwhite.gammaBlue);
    float mixerOrange = float (params->blackwhite.mixerOrange);
    float mixerYellow = float (params->blackwhite.mixerYellow);
    float mixerCyan = float (params->blackwhite.mixerCyan);
    float mixerMagenta = float (params->blackwhite.mixerMagenta);
    float mixerPurple = float (params->blackwhite.mixerPurple);
    int algm = 0;

    if     (params->blackwhite.method == "Desaturation") {
        algm = 0;
    } else if (params->blackwhite.method == "LumEqualizer") {
        algm = 1;
    } else if (params->blackwhite.method == "ChannelMixer") {
        algm = 2;
    }

    float kcorec = 1.f;
    //gamma correction of each channel
    float gamvalr = 125.f;
    float gamvalg = 125.f;
    float gamvalb = 125.f;
    bool computeMixerAuto = params->blackwhite.autoc && (autor < -5000.f);

    if (bwrgam < 0) {
        gamvalr = 100.f;
    }

    if (bwggam < 0) {
        gamvalg = 100.f;
    }

    if (bwbgam < 0) {
        gamvalb = 100.f;
    }

    float gammabwr = 1.f;
    float gammabwg = 1.f;
    float gammabwb = 1.f;
    //if     (params->blackwhite.setting=="Ma" || params->blackwhite.setting=="Mr" || params->blackwhite.setting=="Fr" || params->blackwhite.setting=="Fa")  {
    {
        gammabwr = 1.f - bwrgam / gamvalr;
        gammabwg = 1.f - bwggam / gamvalg;
        gammabwb = 1.f - bwbgam / gamvalb;
    }
    bool hasgammabw = gammabwr != 1.f || gammabwg != 1.f || gammabwb != 1.f;

    if (hasColorToning || blackwhite || (params->dirpyrequalizer.cbdlMethod == "bef" && params->dirpyrequalizer.enabled)) {
        tmpImage = new Imagefloat (working->getWidth(), working->getHeight());
    }

    // For tonecurve histogram
    int toneCurveHistSize = histToneCurve ? histToneCurve.getSize() : 0;
    int histToneCurveCompression = 0;

    if (toneCurveHistSize > 0) {
        histToneCurve.clear();
        histToneCurveCompression = log2 (65536 / toneCurveHistSize);
    }

    // For tonecurve histogram
    const float lumimulf[3] = {static_cast<float> (lumimul[0]), static_cast<float> (lumimul[1]), static_cast<float> (lumimul[2])};


#define TS 112

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        char *buffer;
        char *editIFloatBuffer = nullptr;
        char *editWhateverBuffer = nullptr;

        buffer = (char *) malloc (3 * sizeof (float) * TS * TS + 20 * 64 + 63);
        char *data;
        data = (char*) ( ( uintptr_t (buffer) + uintptr_t (63)) / 64 * 64);

        float *rtemp = (float (*))data;
        float *gtemp = (float (*))         ((char*)rtemp + sizeof (float) * TS * TS + 4 * 64);
        float *btemp = (float (*))         ((char*)gtemp + sizeof (float) * TS * TS + 8 * 64);
        int istart;
        int jstart;
        int tW;
        int tH;

        // Allocating buffer for the PipetteBuffer
        float *editIFloatTmpR = nullptr, *editIFloatTmpG = nullptr, *editIFloatTmpB = nullptr, *editWhateverTmp = nullptr;

        if (editImgFloat) {
            editIFloatBuffer = (char *) malloc (3 * sizeof (float) * TS * TS + 20 * 64 + 63);
            data = (char*) ( ( uintptr_t (editIFloatBuffer) + uintptr_t (63)) / 64 * 64);

            editIFloatTmpR = (float (*))data;
            editIFloatTmpG = (float (*))         ((char*)editIFloatTmpR + sizeof (float) * TS * TS + 4 * 64);
            editIFloatTmpB = (float (*))         ((char*)editIFloatTmpG + sizeof (float) * TS * TS + 8 * 64);
        }

        if (editWhatever) {
            editWhateverBuffer = (char *) malloc (sizeof (float) * TS * TS + 20 * 64 + 63);
            data = (char*) ( ( uintptr_t (editWhateverBuffer) + uintptr_t (63)) / 64 * 64);

            editWhateverTmp = (float (*))data;
        }

        float out_rgbx[4 * TS] ALIGNED16; // Line buffer for CLUT

        LUTu histToneCurveThr;

        if (toneCurveHistSize > 0) {
            histToneCurveThr (toneCurveHistSize);
            histToneCurveThr.clear();
        }

#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif

        for (int ii = 0; ii < working->getHeight(); ii += TS)
            for (int jj = 0; jj < working->getWidth(); jj += TS) {
                istart = ii;
                jstart = jj;
                tH = min (ii + TS, working->getHeight());
                tW = min (jj + TS, working->getWidth());


                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                        rtemp[ti * TS + tj] = working->r (i, j);
                        gtemp[ti * TS + tj] = working->g (i, j);
                        btemp[ti * TS + tj] = working->b (i, j);
                    }
                }

                if (mixchannels) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float r = rtemp[ti * TS + tj];
                            float g = gtemp[ti * TS + tj];
                            float b = btemp[ti * TS + tj];

                            //if (i==100 & j==100) printf("rgbProc input R= %f  G= %f  B= %f  \n",r,g,b);
                            float rmix = (r * chMixRR + g * chMixRG + b * chMixRB) / 100.f;
                            float gmix = (r * chMixGR + g * chMixGG + b * chMixGB) / 100.f;
                            float bmix = (r * chMixBR + g * chMixBG + b * chMixBB) / 100.f;

                            rtemp[ti * TS + tj] = rmix;
                            gtemp[ti * TS + tj] = gmix;
                            btemp[ti * TS + tj] = bmix;
                        }
                    }
                }

                if (processSH || processLCE) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                            float r = rtemp[ti * TS + tj];
                            float g = gtemp[ti * TS + tj];
                            float b = btemp[ti * TS + tj];

                            double mapval = 1.0 + shmap->map[i][j];
                            double factor = 1.0;

                            if (processSH) {
                                if (mapval > h_th) {
                                    factor = (h_th + (100.0 - shHighlights) * (mapval - h_th) / 100.0) / mapval;
                                } else if (mapval < s_th) {
                                    factor = (s_th - (100.0 - shShadows) * (s_th - mapval) / 100.0) / mapval;
                                }
                            }

                            if (processLCE) {
                                double sub = lceamount * (mapval - factor * (r * lumimul[0] + g * lumimul[1] + b * lumimul[2]));
                                rtemp[ti * TS + tj] = factor * r - sub;
                                gtemp[ti * TS + tj] = factor * g - sub;
                                btemp[ti * TS + tj] = factor * b - sub;
                            } else {
                                rtemp[ti * TS + tj] = factor * r;
                                gtemp[ti * TS + tj] = factor * g;
                                btemp[ti * TS + tj] = factor * b;
                            }
                        }
                    }
                }

                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        float r = rtemp[ti * TS + tj];
                        float g = gtemp[ti * TS + tj];
                        float b = btemp[ti * TS + tj];

                        //TODO: proper treatment of out-of-gamut colors
                        //float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
                        float tonefactor = ((r < MAXVALF ? hltonecurve[r] : CurveFactory::hlcurve (exp_scale, comp, hlrange, r) ) +
                                            (g < MAXVALF ? hltonecurve[g] : CurveFactory::hlcurve (exp_scale, comp, hlrange, g) ) +
                                            (b < MAXVALF ? hltonecurve[b] : CurveFactory::hlcurve (exp_scale, comp, hlrange, b) ) ) / 3.0;

                        // note: tonefactor includes exposure scaling, that is here exposure slider and highlight compression takes place
                        rtemp[ti * TS + tj] = r * tonefactor;
                        gtemp[ti * TS + tj] = g * tonefactor;
                        btemp[ti * TS + tj] = b * tonefactor;
                    }
                }

                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        float r = rtemp[ti * TS + tj];
                        float g = gtemp[ti * TS + tj];
                        float b = btemp[ti * TS + tj];

                        //shadow tone curve
                        float Y = (0.299f * r + 0.587f * g + 0.114f * b);
                        float tonefactor = shtonecurve[Y];
                        rtemp[ti * TS + tj] = rtemp[ti * TS + tj] * tonefactor;
                        gtemp[ti * TS + tj] = gtemp[ti * TS + tj] * tonefactor;
                        btemp[ti * TS + tj] = btemp[ti * TS + tj] * tonefactor;
                    }
                }

                if (dcpProf) {
                    dcpProf->step2ApplyTile (rtemp, gtemp, btemp, tW - jstart, tH - istart, TS, asIn);
                }

                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                        float r = rtemp[ti * TS + tj];
                        float g = gtemp[ti * TS + tj];
                        float b = btemp[ti * TS + tj];

                        // clip out of gamut colors, without distorting color too bad
                        if (r < 0) {
                            r = 0;
                        }

                        if (g < 0) {
                            g = 0;
                        }

                        if (b < 0) {
                            b = 0;
                        }

                        if (r > 65535 || g > 65535 || b > 65535) {
                            filmlike_clip (&r, &g, &b);
                        }

                        rtemp[ti * TS + tj] = r;
                        gtemp[ti * TS + tj] = g;
                        btemp[ti * TS + tj] = b;
                    }
                }

                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        //brightness/contrast
                        rtemp[ti * TS + tj] = tonecurve[ rtemp[ti * TS + tj] ];
                        gtemp[ti * TS + tj] = tonecurve[ gtemp[ti * TS + tj] ];
                        btemp[ti * TS + tj] = tonecurve[ btemp[ti * TS + tj] ];

                        if (histToneCurveThr) {
                            int y = CLIP<int> (lumimulf[0] * Color::gamma2curve[rtemp[ti * TS + tj]] + lumimulf[1] * Color::gamma2curve[gtemp[ti * TS + tj]] + lumimulf[2] * Color::gamma2curve[btemp[ti * TS + tj]]);
                            histToneCurveThr[y >> histToneCurveCompression]++;
                        }
                    }
                }

                if (editID == EUID_ToneCurve1) {  // filling the pipette buffer
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editIFloatTmpR[ti * TS + tj] = Color::gamma2curve[rtemp[ti * TS + tj]] / 65535.f;
                            editIFloatTmpG[ti * TS + tj] = Color::gamma2curve[gtemp[ti * TS + tj]] / 65535.f;
                            editIFloatTmpB[ti * TS + tj] = Color::gamma2curve[btemp[ti * TS + tj]] / 65535.f;
                        }
                    }
                }

                if (hasToneCurve1) {
                    if (curveMode == ToneCurveParams::TC_MODE_STD) { // Standard
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                const StandardToneCurve& userToneCurve = static_cast<const StandardToneCurve&> (customToneCurve1);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode == ToneCurveParams::TC_MODE_FILMLIKE) { // Adobe like
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                const AdobeToneCurve& userToneCurve = static_cast<const AdobeToneCurve&> (customToneCurve1);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode == ToneCurveParams::TC_MODE_SATANDVALBLENDING) { // apply the curve on the saturation and value channels
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                const SatAndValueBlendingToneCurve& userToneCurve = static_cast<const SatAndValueBlendingToneCurve&> (customToneCurve1);
                                rtemp[ti * TS + tj] = CLIP<float> (rtemp[ti * TS + tj]);
                                gtemp[ti * TS + tj] = CLIP<float> (gtemp[ti * TS + tj]);
                                btemp[ti * TS + tj] = CLIP<float> (btemp[ti * TS + tj]);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode == ToneCurveParams::TC_MODE_WEIGHTEDSTD) { // apply the curve to the rgb channels, weighted
                        const WeightedStdToneCurve& userToneCurve = static_cast<const WeightedStdToneCurve&> (customToneCurve1);

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                rtemp[ti * TS + tj] = CLIP<float> (rtemp[ti * TS + tj]);
                                gtemp[ti * TS + tj] = CLIP<float> (gtemp[ti * TS + tj]);
                                btemp[ti * TS + tj] = CLIP<float> (btemp[ti * TS + tj]);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode == ToneCurveParams::TC_MODE_LUMINANCE) { // apply the curve to the luminance channel
                        const LuminanceToneCurve& userToneCurve = static_cast<const LuminanceToneCurve&> (customToneCurve1);

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                rtemp[ti * TS + tj] = CLIP<float> (rtemp[ti * TS + tj]);
                                gtemp[ti * TS + tj] = CLIP<float> (gtemp[ti * TS + tj]);
                                btemp[ti * TS + tj] = CLIP<float> (btemp[ti * TS + tj]);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode == ToneCurveParams::TC_MODE_PERCEPTUAL) { // apply curve while keeping color appearance constant
                        const PerceptualToneCurve& userToneCurve = static_cast<const PerceptualToneCurve&> (customToneCurve1);

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                rtemp[ti * TS + tj] = CLIP<float> (rtemp[ti * TS + tj]);
                                gtemp[ti * TS + tj] = CLIP<float> (gtemp[ti * TS + tj]);
                                btemp[ti * TS + tj] = CLIP<float> (btemp[ti * TS + tj]);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], ptc1ApplyState);
                            }
                        }
                    }
                }

                if (editID == EUID_ToneCurve2) {  // filling the pipette buffer
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editIFloatTmpR[ti * TS + tj] = Color::gamma2curve[rtemp[ti * TS + tj]] / 65535.f;
                            editIFloatTmpG[ti * TS + tj] = Color::gamma2curve[gtemp[ti * TS + tj]] / 65535.f;
                            editIFloatTmpB[ti * TS + tj] = Color::gamma2curve[btemp[ti * TS + tj]] / 65535.f;
                        }
                    }
                }

                if (hasToneCurve2) {
                    if (curveMode2 == ToneCurveParams::TC_MODE_STD) { // Standard
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                const StandardToneCurve& userToneCurve = static_cast<const StandardToneCurve&> (customToneCurve2);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode2 == ToneCurveParams::TC_MODE_FILMLIKE) { // Adobe like
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                const AdobeToneCurve& userToneCurve = static_cast<const AdobeToneCurve&> (customToneCurve2);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode2 == ToneCurveParams::TC_MODE_SATANDVALBLENDING) { // apply the curve on the saturation and value channels
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                const SatAndValueBlendingToneCurve& userToneCurve = static_cast<const SatAndValueBlendingToneCurve&> (customToneCurve2);
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode2 == ToneCurveParams::TC_MODE_WEIGHTEDSTD) { // apply the curve to the rgb channels, weighted
                        const WeightedStdToneCurve& userToneCurve = static_cast<const WeightedStdToneCurve&> (customToneCurve2);

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode2 == ToneCurveParams::TC_MODE_LUMINANCE) { // apply the curve to the luminance channel
                        const LuminanceToneCurve& userToneCurve = static_cast<const LuminanceToneCurve&> (customToneCurve2);

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    } else if (curveMode2 == ToneCurveParams::TC_MODE_PERCEPTUAL) { // apply curve while keeping color appearance constant
                        const PerceptualToneCurve& userToneCurve = static_cast<const PerceptualToneCurve&> (customToneCurve2);

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                userToneCurve.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], ptc2ApplyState);
                            }
                        }
                    }
                }

                if (editID == EUID_RGB_R) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[rtemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                } else if (editID == EUID_RGB_G) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[gtemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                } else if (editID == EUID_RGB_B) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[btemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                }

                if (rCurve || gCurve || bCurve) { // if any of the RGB curves is engaged
                    if (!params->rgbCurves.lumamode) { // normal RGB mode

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                // individual R tone curve
                                if (rCurve) {
                                    rtemp[ti * TS + tj] = rCurve[ rtemp[ti * TS + tj] ];
                                }

                                // individual G tone curve
                                if (gCurve) {
                                    gtemp[ti * TS + tj] = gCurve[ gtemp[ti * TS + tj] ];
                                }

                                // individual B tone curve
                                if (bCurve) {
                                    btemp[ti * TS + tj] = bCurve[ btemp[ti * TS + tj] ];
                                }
                            }
                        }
                    } else { //params->rgbCurves.lumamode==true (Luminosity mode)
                        // rCurve.dump("r_curve");//debug

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                // rgb values before RGB curves
                                float r = rtemp[ti * TS + tj] ;
                                float g = gtemp[ti * TS + tj] ;
                                float b = btemp[ti * TS + tj] ;
                                //convert to Lab to get a&b before RGB curves
                                float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
                                float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                                float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;

                                float fx = x < MAXVALF ? Color::cachef[x] : 327.68f * std::cbrt (x / MAXVALF);
                                float fy = y < MAXVALF ? Color::cachef[y] : 327.68f * std::cbrt (y / MAXVALF);
                                float fz = z < MAXVALF ? Color::cachef[z] : 327.68f * std::cbrt (z / MAXVALF);

                                float a_1 = 500.0f * (fx - fy);
                                float b_1 = 200.0f * (fy - fz);

                                // rgb values after RGB curves
                                if (rCurve) {
                                    float rNew = rCurve[r];
                                    r += (rNew - r) * equalR;
                                }

                                if (gCurve) {
                                    float gNew = gCurve[g];
                                    g += (gNew - g) * equalG;
                                }

                                if (bCurve) {
                                    float bNew = bCurve[b];
                                    b += (bNew - b) * equalB;
                                }

                                // Luminosity after
                                // only Luminance in Lab
                                float newy = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                                float newfy = newy < MAXVALF ? Color::cachef[newy] : 327.68f * std::cbrt (newy / MAXVALF);
                                float L_2 = 116.0f * newfy - 5242.88f;

                                //gamut control
                                if (settings->rgbcurveslumamode_gamut) {
                                    float Lpro = L_2 / 327.68f;
                                    float Chpro = sqrtf (SQR (a_1) + SQR (b_1)) / 327.68f;
                                    float HH = NAN; // we set HH to NAN, because then it will be calculated in Color::gamutLchonly only if needed
//                                    float HH = xatan2f(b_1, a_1);
                                    // According to mathematical laws we can get the sin and cos of HH by simple operations even if we don't calculate HH
                                    float2 sincosval;

                                    if (Chpro == 0.0f) {
                                        sincosval.y = 1.0f;
                                        sincosval.x = 0.0f;
                                    } else {
                                        sincosval.y = a_1 / (Chpro * 327.68f);
                                        sincosval.x = b_1 / (Chpro * 327.68f);
                                    }

#ifdef _DEBUG
                                    bool neg = false;
                                    bool more_rgb = false;
                                    //gamut control : Lab values are in gamut
                                    Color::gamutLchonly (HH, sincosval, Lpro, Chpro, rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                                    //gamut control : Lab values are in gamut
                                    Color::gamutLchonly (HH, sincosval, Lpro, Chpro, rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], wip, highlight, 0.15f, 0.96f);
#endif
                                    //end of gamut control
                                } else {
                                    float x_, y_, z_;
                                    //calculate RGB with L_2 and old value of a and b
                                    Color::Lab2XYZ (L_2, a_1, b_1, x_, y_, z_) ;
                                    Color::xyz2rgb (x_, y_, z_, rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], wip);
                                }
                            }
                        }
                    }
                }

                if (editID == EUID_HSV_H || editID == EUID_HSV_S || editID == EUID_HSV_V) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float h, s, v;
                            Color::rgb2hsv (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], h, s, v);
                            editWhateverTmp[ti * TS + tj] = h;
                        }
                    }
                }

                if (sat != 0 || hCurveEnabled || sCurveEnabled || vCurveEnabled) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                            const float satby100 = sat / 100.f;
                            float r = rtemp[ti * TS + tj];
                            float g = gtemp[ti * TS + tj];
                            float b = btemp[ti * TS + tj];
                            float h, s, v;
                            Color::rgb2hsv (r, g, b, h, s, v);

                            if (sat > 0) {
                                s = (1.f - satby100) * s + satby100 * (1.f - SQR (SQR (1.f - min (s, 1.0f))));

                                if (s < 0.f) {
                                    s = 0.f;
                                }
                            } else { /*if (sat < 0)*/
                                s *= 1.f + satby100;
                            }

                            //HSV equalizer
                            if (hCurveEnabled) {
                                h = (hCurve->getVal (double (h)) - 0.5) * 2.f + h;

                                if (h > 1.0f) {
                                    h -= 1.0f;
                                } else if (h < 0.0f) {
                                    h += 1.0f;
                                }
                            }

                            if (sCurveEnabled) {
                                //shift saturation
                                float satparam = (sCurve->getVal (double (h)) - 0.5) * 2;

                                if (satparam > 0.00001f) {
                                    s = (1.f - satparam) * s + satparam * (1.f - SQR (1.f - min (s, 1.0f)));

                                    if (s < 0.f) {
                                        s = 0.f;
                                    }
                                } else if (satparam < -0.00001f) {
                                    s *= 1.f + satparam;
                                }

                            }

                            if (vCurveEnabled) {
                                if (v < 0) {
                                    v = 0;    // important
                                }

                                //shift value
                                float valparam = vCurve->getVal ((double)h) - 0.5f;
                                valparam *= (1.f - SQR (SQR (1.f - min (s, 1.0f))));

                                if (valparam > 0.00001f) {
                                    v = (1.f - valparam) * v + valparam * (1.f - SQR (1.f - min (v, 1.0f))); // SQR (SQR  to increase action and avoid artefacts

                                    if (v < 0) {
                                        v = 0;
                                    }
                                } else {
                                    if (valparam < -0.00001f) {
                                        v *= (1.f + valparam);    //1.99 to increase action
                                    }
                                }

                            }

                            Color::hsv2rgb (h, s, v, rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                        }
                    }
                }

                if (isProPhoto) { // this is a hack to avoid the blue=>black bug (Issue 2141)
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float r = rtemp[ti * TS + tj];
                            float g = gtemp[ti * TS + tj];

                            if (r == 0.0f || g == 0.0f) {
                                float b = btemp[ti * TS + tj];
                                float h, s, v;
                                Color::rgb2hsv (r, g, b, h, s, v);
                                s *= 0.99f;
                                Color::hsv2rgb (h, s, v, rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                            }
                        }
                    }
                }

                if (hasColorToning && !blackwhite) {
                    if (params->colorToning.method == "Splitlr") {
                        float balanS, balanH;
                        float reducac = 0.4f;
                        int preser = 0;

                        if (params->colorToning.lumamode) {
                            preser = 1;
                        }

                        balanS = 1.f + Balan / 100.f; //balan between 0 and 2
                        balanH = 1.f - Balan / 100.f;
                        float rh, gh, bh;
                        float rl, gl, bl;
                        float xh, yh, zh;
                        float xl, yl, zl;
                        float iplow, iphigh;
                        iplow = (float)ctColorCurve.low;
                        iphigh = (float)ctColorCurve.high;
                        //2 colours
                        ctColorCurve.getVal (iphigh, xh, yh, zh);
                        ctColorCurve.getVal (iplow, xl, yl, zl);

                        Color::xyz2rgb (xh, yh, zh, rh, gh, bh, wip);
                        Color::xyz2rgb (xl, yl, zl, rl, gl, bl, wip);
                        //reteave rgb value with s and l =1
                        retreavergb (rl, gl, bl);
                        retreavergb (rh, gh, bh);
                        //printf("rl=%f gl=%f bl=%f\n",rl,gl,bl);

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                float r = rtemp[ti * TS + tj];
                                float g = gtemp[ti * TS + tj];
                                float b = btemp[ti * TS + tj];
                                float ro, go, bo;
                                int mode = 0;
                                toning2col (r, g, b, ro, go, bo, iplow, iphigh, rl, gl, bl, rh, gh, bh, SatLow, SatHigh, balanS, balanH, reducac, mode, preser, strProtect);
                                rtemp[ti * TS + tj] = ro;
                                gtemp[ti * TS + tj] = go;
                                btemp[ti * TS + tj] = bo;
                            }
                        }
                    }

                    // color toning with colour
                    else if (params->colorToning.method == "Splitco") {
                        /*
                        #if 1
                                            for (int i=istart,ti=0; i<tH; i++,ti++) {
                                                for (int j=jstart,tj=0; j<tW; j++,tj++) {
                                                    float r = rtemp[ti*TS+tj];
                                                    float g = gtemp[ti*TS+tj];
                                                    float b = btemp[ti*TS+tj];
                                                    float ro,go,bo;
                                                    labtoning (r, g, b, ro, go, bo, 1, 0, 1, satLimit, satLimitOpacity, ctColorCurve, ctOpacityCurve, clToningcurve, cl2Toningcurve, 0.f, 1.f, wp, wip);
                                                    rtemp[ti*TS+tj] = CLIP(ro);//I used CLIP because there is a little bug in gamutLchonly that return 65536.ii intead of 65535 ==> crash
                                                    gtemp[ti*TS+tj] = CLIP(go);
                                                    btemp[ti*TS+tj] = CLIP(bo);
                                                }
                                            }
                        #else
                        */
                        float reducac = 0.3f;
                        int preser = 0;

                        //bool execbal = params->colorToning.method=="Splitbal";
                        if (params->colorToning.lumamode) {
                            preser = 1;
                        }

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                float r = rtemp[ti * TS + tj];
                                float g = gtemp[ti * TS + tj];
                                float b = btemp[ti * TS + tj];
                                float lumbefore = 0.299f * r + 0.587f * g + 0.114f * b;
                                float ro, go, bo;
                                int mode = 0;
                                toningsmh (r, g, b, ro, go, bo, RedLow, GreenLow, BlueLow, RedMed, GreenMed, BlueMed, RedHigh, GreenHigh, BlueHigh, reducac, mode, preser, strProtect);
                                float lumafter = 0.299f * ro + 0.587f * go + 0.114f * bo;
                                float preserv = 1.f;

                                if (preser == 1) {
                                    preserv = lumbefore / lumafter;
                                }

                                ro *= preserv;
                                go *= preserv;
                                bo *= preserv;
                                ro = CLIP (ro);
                                go = CLIP (go);
                                bo = CLIP (bo);
                                rtemp[ti * TS + tj] = ro;
                                gtemp[ti * TS + tj] = go;
                                btemp[ti * TS + tj] = bo;
                            }
                        }

//#endif
                    }

                    //colortoning with shift color XYZ or Lch
                    else if (params->colorToning.method == "Lab" && opautili) {
                        int algm = 0;
                        bool twocol = true;//true=500 color   false=2 color
                        int metchrom = 0;

                        if      (params->colorToning.twocolor == "Std"  ) {
                            metchrom = 0;
                        } else if (params->colorToning.twocolor == "All"  ) {
                            metchrom = 1;
                        } else if (params->colorToning.twocolor == "Separ") {
                            metchrom = 2;
                        } else if (params->colorToning.twocolor == "Two"  ) {
                            metchrom = 3;
                        }

                        if (metchrom == 3) {
                            twocol = false;
                        }

                        float iplow = 0.f, iphigh = 0.f;

                        if (!twocol) {
                            iplow = (float)ctColorCurve.low;
                            iphigh = (float)ctColorCurve.high;
                        }

                        int twoc = 0; //integer instead of bool to let more possible choice...other than 2 and 500.

                        if (!twocol) {
                            twoc = 0;    // 2 colours
                        } else {
                            twoc = 1;    // 500 colours
                        }

                        if      (params->colorToning.method == "Lab") {
                            algm = 1;
                        } else if (params->colorToning.method == "Lch") {
                            algm = 2;    //in case of
                        }

                        if (algm <= 2) {
                            for (int i = istart, ti = 0; i < tH; i++, ti++) {
                                for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                    float r = rtemp[ti * TS + tj];
                                    float g = gtemp[ti * TS + tj];
                                    float b = btemp[ti * TS + tj];
                                    float ro, go, bo;
                                    labtoning (r, g, b, ro, go, bo, algm, metchrom, twoc, satLimit, satLimitOpacity, ctColorCurve, ctOpacityCurve, clToningcurve, cl2Toningcurve, iplow, iphigh, wp, wip);
                                    rtemp[ti * TS + tj] = CLIP (ro); //I used CLIP because there is a little bug in gamutLchonly that return 65536.ii intead of 65535 ==> crash
                                    gtemp[ti * TS + tj] = CLIP (go);
                                    btemp[ti * TS + tj] = CLIP (bo);
                                }
                            }
                        }
                    } else if (params->colorToning.method.substr (0, 3) == "RGB" && opautili) {
                        // color toning
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                float r = rtemp[ti * TS + tj];
                                float g = gtemp[ti * TS + tj];
                                float b = btemp[ti * TS + tj];

                                // Luminance = (0.299f*r + 0.587f*g + 0.114f*b)

                                float h, s, l;
                                Color::rgb2hsl (r, g, b, h, s, l);

                                float l_ = Color::gamma_srgb (l * 65535.f) / 65535.f;

                                // get the opacity and tweak it to preserve saturated colors
                                float opacity;

                                if (ctOpacityCurve) {
                                    opacity = (1.f - min<float> (s / satLimit, 1.f) * (1.f - satLimitOpacity)) * ctOpacityCurve.lutOpacityCurve[l_ * 500.f];
                                }

                                if (!ctOpacityCurve) {
                                    opacity = 0.f;
                                }

                                float r2, g2, b2;
                                ctColorCurve.getVal (l_, r2, g2, b2); // get the color from the color curve

                                float h2, s2, l2;
                                Color::rgb2hsl (r2, g2, b2, h2, s2, l2); // transform this new color to hsl

                                Color::hsl2rgb (h2, s + ((1.f - s) * (1.f - l) * 0.7f), l, r2, g2, b2);

                                rtemp[ti * TS + tj] = r + (r2 - r) * opacity; // merge the color to the old color, depending on the opacity
                                gtemp[ti * TS + tj] = g + (g2 - g) * opacity;
                                btemp[ti * TS + tj] = b + (b2 - b) * opacity;
                            }
                        }
                    }
                }

                // filling the pipette buffer
                if (editID == EUID_BlackWhiteBeforeCurve) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editIFloatTmpR[ti * TS + tj] = Color::gamma2curve[rtemp[ti * TS + tj]] / 65535.f;
                            editIFloatTmpG[ti * TS + tj] = Color::gamma2curve[gtemp[ti * TS + tj]] / 65535.f;
                            editIFloatTmpB[ti * TS + tj] = Color::gamma2curve[btemp[ti * TS + tj]] / 65535.f;
                        }
                    }
                } else if (editID == EUID_BlackWhiteLuminance) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float X, Y, Z, L, aa, bb;
                            //rgb=>lab
                            Color::rgbxyz (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], X, Y, Z, wp);
                            //convert Lab
                            Color::XYZ2Lab (X, Y, Z, L, aa, bb);
                            //end rgb=>lab
                            float HH = xatan2f (bb, aa); // HH hue in -3.141  +3.141

                            editWhateverTmp[ti * TS + tj] = float (Color::huelab_to_huehsv2 (HH));
                        }
                    }
                }

                //black and white
                if (blackwhite) {
                    if (hasToneCurvebw1) {
                        if (beforeCurveMode == BlackWhiteParams::TC_MODE_STD_BW) { // Standard
                            for (int i = istart, ti = 0; i < tH; i++, ti++) {
                                for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                    const StandardToneCurvebw& userToneCurvebw = static_cast<const StandardToneCurvebw&> (customToneCurvebw1);
                                    userToneCurvebw.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                                }
                            }
                        } else if (beforeCurveMode == BlackWhiteParams::TC_MODE_FILMLIKE_BW) { // Adobe like
                            for (int i = istart, ti = 0; i < tH; i++, ti++) {
                                for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                    const AdobeToneCurvebw& userToneCurvebw = static_cast<const AdobeToneCurvebw&> (customToneCurvebw1);
                                    userToneCurvebw.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                                }
                            }
                        } else if (beforeCurveMode == BlackWhiteParams::TC_MODE_SATANDVALBLENDING_BW) { // apply the curve on the saturation and value channels
                            for (int i = istart, ti = 0; i < tH; i++, ti++) {
                                for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                    const SatAndValueBlendingToneCurvebw& userToneCurvebw = static_cast<const SatAndValueBlendingToneCurvebw&> (customToneCurvebw1);
                                    rtemp[ti * TS + tj] = CLIP<float> (rtemp[ti * TS + tj]);
                                    gtemp[ti * TS + tj] = CLIP<float> (gtemp[ti * TS + tj]);
                                    btemp[ti * TS + tj] = CLIP<float> (btemp[ti * TS + tj]);
                                    userToneCurvebw.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                                }
                            }
                        } else if (beforeCurveMode == BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW) { // apply the curve to the rgb channels, weighted
                            for (int i = istart, ti = 0; i < tH; i++, ti++) {
                                for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                    const WeightedStdToneCurvebw& userToneCurvebw = static_cast<const WeightedStdToneCurvebw&> (customToneCurvebw1);
                                    rtemp[ti * TS + tj] = CLIP<float> (rtemp[ti * TS + tj]);
                                    gtemp[ti * TS + tj] = CLIP<float> (gtemp[ti * TS + tj]);
                                    btemp[ti * TS + tj] = CLIP<float> (btemp[ti * TS + tj]);

                                    userToneCurvebw.Apply (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                                }
                            }
                        }
                    }

                    if (algm == 0) { //lightness
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                                float r = rtemp[ti * TS + tj];
                                float g = gtemp[ti * TS + tj];
                                float b = btemp[ti * TS + tj];

                                // --------------------------------------------------

                                // Method 1: Luminosity (code taken from Gimp)
                                /*
                                float maxi = max(r, g, b);
                                float mini = min(r, g, b);
                                r = g = b = (maxi+mini)/2;
                                */

                                // Method 2: Luminance (former RT code)
                                r = g = b = (0.299f * r + 0.587f * g + 0.114f * b);

                                // --------------------------------------------------

#ifndef __SSE2__

                                //gamma correction: pseudo TRC curve
                                if (hasgammabw) {
                                    Color::trcGammaBW (r, g, b, gammabwr, gammabwg, gammabwb);
                                }

#endif
                                rtemp[ti * TS + tj] = r;
                                gtemp[ti * TS + tj] = g;
                                btemp[ti * TS + tj] = b;
                            }

#ifdef __SSE2__

                            if (hasgammabw) {
                                //gamma correction: pseudo TRC curve
                                Color::trcGammaBWRow (&rtemp[ti * TS], &gtemp[ti * TS], &btemp[ti * TS], tW - jstart, gammabwr, gammabwg, gammabwb);
                            }

#endif

                        }
                    } else if (algm == 1) { //Luminance mixer in Lab mode to avoid artefacts
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                //rgb => xyz
                                float X, Y, Z;
                                Color::rgbxyz (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], X, Y, Z, wp);
                                //xyz => Lab
                                float L, aa, bb;
                                Color::XYZ2Lab (X, Y, Z, L, aa, bb);
                                float CC = sqrtf (SQR (aa) + SQR (bb)) / 327.68f; //CC chromaticity in 0..180 or more
                                float HH = xatan2f (bb, aa); // HH hue in -3.141  +3.141
                                float2 sincosval;

                                if (CC == 0.0f) {
                                    sincosval.y = 1.f;
                                    sincosval.x = 0.0f;
                                } else {
                                    sincosval.y = aa / (CC * 327.68f);
                                    sincosval.x = bb / (CC * 327.68f);
                                }

                                if (bwlCurveEnabled) {
                                    L /= 32768.f;
                                    double hr = Color::huelab_to_huehsv2 (HH);
                                    float valparam = float ((bwlCurve->getVal (hr) - 0.5f) * 2.0f); //get l_r=f(H)
                                    float kcc = (CC / 70.f); //take Chroma into account...70 "middle" of chromaticity (arbitrary and simple), one can imagine other algorithme
                                    //reduct action for low chroma and increase action for high chroma
                                    valparam *= kcc;

                                    if (valparam > 0.f) {
                                        L = (1.f - valparam) * L + valparam * (1.f - SQR (SQR (SQR (SQR (1.f - min (L, 1.0f)))))); // SQR (SQR((SQR)  to increase action in low light
                                    } else {
                                        L *= (1.f + valparam);    //for negative
                                    }

                                    L *= 32768.f;
                                }

                                float RR, GG, BB;
                                L /= 327.68f;
#ifdef _DEBUG
                                bool neg = false;
                                bool more_rgb = false;
                                //gamut control : Lab values are in gamut
                                Color::gamutLchonly (HH, sincosval, L, CC, RR, GG, BB, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                                //gamut control : Lab values are in gamut
                                Color::gamutLchonly (HH, sincosval, L, CC, RR, GG, BB, wip, highlight, 0.15f, 0.96f);
#endif
                                L *= 327.68f;
                                //convert l => rgb
                                Color::L2XYZ (L, X, Y, Z);
                                float newRed; // We use the red channel for bw
                                Color::xyz2r (X, Y, Z, newRed, wip);
                                rtemp[ti * TS + tj] = gtemp[ti * TS + tj] = btemp[ti * TS + tj] = newRed;
#ifndef __SSE2__

                                if (hasgammabw) {
                                    //gamma correction: pseudo TRC curve
                                    Color::trcGammaBW (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], gammabwr, gammabwg, gammabwb);
                                }

#endif
                            }

#ifdef __SSE2__

                            if (hasgammabw) {
                                //gamma correction: pseudo TRC curve
                                Color::trcGammaBWRow (&rtemp[ti * TS], &gtemp[ti * TS], &btemp[ti * TS], tW - jstart, gammabwr, gammabwg, gammabwb);
                            }

#endif
                        }
                    }
                }


                // Film Simulations
                if (hald_clut) {

                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        if (!clutAndWorkingProfilesAreSame) {
                            // Convert from working to clut profile
                            int j = jstart;
                            int tj = 0;

#ifdef __SSE2__

                            for (; j < tW - 3; j += 4, tj += 4) {
                                vfloat sourceR = LVF (rtemp[ti * TS + tj]);
                                vfloat sourceG = LVF (gtemp[ti * TS + tj]);
                                vfloat sourceB = LVF (btemp[ti * TS + tj]);

                                vfloat x;
                                vfloat y;
                                vfloat z;
                                Color::rgbxyz (sourceR, sourceG, sourceB, x, y, z, v_work2xyz);
                                Color::xyz2rgb (x, y, z, sourceR, sourceG, sourceB, v_xyz2clut);

                                STVF (rtemp[ti * TS + tj], sourceR);
                                STVF (gtemp[ti * TS + tj], sourceG);
                                STVF (btemp[ti * TS + tj], sourceB);
                            }

#endif

                            for (; j < tW; j++, tj++) {
                                float &sourceR = rtemp[ti * TS + tj];
                                float &sourceG = gtemp[ti * TS + tj];
                                float &sourceB = btemp[ti * TS + tj];

                                float x, y, z;
                                Color::rgbxyz ( sourceR, sourceG, sourceB, x, y, z, wprof );
                                Color::xyz2rgb (x, y, z, sourceR, sourceG, sourceB, xyz2clut);
                            }
                        }

                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float &sourceR = rtemp[ti * TS + tj];
                            float &sourceG = gtemp[ti * TS + tj];
                            float &sourceB = btemp[ti * TS + tj];

                            // Apply gamma sRGB (default RT)
                            sourceR = Color::gamma_srgbclipped (sourceR);
                            sourceG = Color::gamma_srgbclipped (sourceG);
                            sourceB = Color::gamma_srgbclipped (sourceB);
                        }

                        const std::size_t line_offset = ti * TS;
                        hald_clut->getRGB (
                            film_simulation_strength,
                            std::min (TS, tW - jstart),
                            rtemp + line_offset,
                            gtemp + line_offset,
                            btemp + line_offset,
                            out_rgbx
                        );

                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float &sourceR = rtemp[ti * TS + tj];
                            float &sourceG = gtemp[ti * TS + tj];
                            float &sourceB = btemp[ti * TS + tj];

                            // Apply inverse gamma sRGB
                            sourceR = Color::igamma_srgb (out_rgbx[tj * 4 + 0]);
                            sourceG = Color::igamma_srgb (out_rgbx[tj * 4 + 1]);
                            sourceB = Color::igamma_srgb (out_rgbx[tj * 4 + 2]);
                        }

                        if (!clutAndWorkingProfilesAreSame) {
                            // Convert from clut to working profile
                            int j = jstart;
                            int tj = 0;

#ifdef __SSE2__

                            for (; j < tW - 3; j += 4, tj += 4) {
                                vfloat sourceR = LVF (rtemp[ti * TS + tj]);
                                vfloat sourceG = LVF (gtemp[ti * TS + tj]);
                                vfloat sourceB = LVF (btemp[ti * TS + tj]);

                                vfloat x;
                                vfloat y;
                                vfloat z;
                                Color::rgbxyz (sourceR, sourceG, sourceB, x, y, z, v_clut2xyz);
                                Color::xyz2rgb (x, y, z, sourceR, sourceG, sourceB, v_xyz2work);

                                STVF (rtemp[ti * TS + tj], sourceR);
                                STVF (gtemp[ti * TS + tj], sourceG);
                                STVF (btemp[ti * TS + tj], sourceB);
                            }

#endif

                            for (; j < tW; j++, tj++) {
                                float &sourceR = rtemp[ti * TS + tj];
                                float &sourceG = gtemp[ti * TS + tj];
                                float &sourceB = btemp[ti * TS + tj];

                                float x, y, z;
                                Color::rgbxyz (sourceR, sourceG, sourceB, x, y, z, clut2xyz);
                                Color::xyz2rgb ( x, y, z, sourceR, sourceG, sourceB, wiprof );
                            }
                        }
                    }
                }


                if (!blackwhite) {
                    // ready, fill lab
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                            // filling the pipette buffer by the content of the temp pipette buffers
                            if (editImgFloat) {
                                editImgFloat->r (i, j) = editIFloatTmpR[ti * TS + tj];
                                editImgFloat->g (i, j) = editIFloatTmpG[ti * TS + tj];
                                editImgFloat->b (i, j) = editIFloatTmpB[ti * TS + tj];
                            } else if (editWhatever) {
                                editWhatever->v (i, j) = editWhateverTmp[ti * TS + tj];
                            }

                            float r = rtemp[ti * TS + tj];
                            float g = gtemp[ti * TS + tj];
                            float b = btemp[ti * TS + tj];

                            float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
                            float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                            float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;

                            float fx, fy, fz;

                            fx = (x < 65535.0f ? Color::cachef[x] : 327.68f * std::cbrt (x / MAXVALF));
                            fy = (y < 65535.0f ? Color::cachef[y] : 327.68f * std::cbrt (y / MAXVALF));
                            fz = (z < 65535.0f ? Color::cachef[z] : 327.68f * std::cbrt (z / MAXVALF));

                            lab->L[i][j] = (116.0f *  fy - 5242.88f); //5242.88=16.0*327.68;
                            lab->a[i][j] = (500.0f * (fx - fy) );
                            lab->b[i][j] = (200.0f * (fy - fz) );

                            //test for color accuracy
                            /*
                            float fy = (0.00862069 * lab->L[i][j])/327.68 + 0.137932; // (L+16)/116
                            float fx = (0.002 * lab->a[i][j])/327.68 + fy;
                            float fz = fy - (0.005 * lab->b[i][j])/327.68;

                            float x_ = 65535*Lab2xyz(fx)*Color::D50x;
                            float y_ = 65535*Lab2xyz(fy);
                            float z_ = 65535*Lab2xyz(fz)*Color::D50z;

                            int R,G,B;
                            xyz2srgb(x_,y_,z_,R,G,B);
                            r=(float)R; g=(float)G; b=(float)B;
                            float xxx=1;
                            */
                        }
                    }
                } else { // black & white
                    // Auto channel mixer needs whole image, so we now copy to tmpImage and close the tiled processing
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            // filling the pipette buffer by the content of the temp pipette buffers
                            if (editImgFloat) {
                                editImgFloat->r (i, j) = editIFloatTmpR[ti * TS + tj];
                                editImgFloat->g (i, j) = editIFloatTmpG[ti * TS + tj];
                                editImgFloat->b (i, j) = editIFloatTmpB[ti * TS + tj];
                            } else if (editWhatever) {
                                editWhatever->v (i, j) = editWhateverTmp[ti * TS + tj];
                            }

                            tmpImage->r (i, j) = rtemp[ti * TS + tj];
                            tmpImage->g (i, j) = gtemp[ti * TS + tj];
                            tmpImage->b (i, j) = btemp[ti * TS + tj];
                        }
                    }
                }
            }

        free (buffer);

        if (editIFloatBuffer) {
            free (editIFloatBuffer);
        }

        if (editWhateverBuffer) {
            free (editWhateverBuffer);
        }

#ifdef _OPENMP
        #pragma omp critical
        {
            if (toneCurveHistSize > 0) {
                histToneCurve += histToneCurveThr;
            }
        }
#endif // _OPENMP
    }

    // starting a new tile processing with a 'reduction' clause for the auto mixer computing
    if (blackwhite) {//channel-mixer
        int tW = working->getWidth();
        int tH = working->getHeight();

        if (algm == 2) { //channel-mixer
            //end auto chmix
            if (computeMixerAuto) {
                // auto channel-mixer
                double nr = 0;
                double ng = 0;
                double nb = 0;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic, 16) reduction(+:nr,ng,nb)
#endif

                for (int i = 0; i < tH; i++) {
                    for (int j = 0; j < tW; j++) {
                        nr += tmpImage->r (i, j);
                        ng += tmpImage->g (i, j);
                        nb += tmpImage->b (i, j);
                    }
                }

                double srgb = nr + ng + nb;
                double knr = srgb / nr;
                double kng = srgb / ng;
                double knb = srgb / nb;
                double sk = knr + kng + knb;
                autor = (float) (100.0 * knr / sk);
                autog = (float) (100.0 * kng / sk);
                autob = (float) (100.0 * knb / sk);

            }

            if (params->blackwhite.autoc) {
                // auto channel-mixer
                bwr = autor;
                bwg = autog;
                bwb = autob;
                mixerOrange  = 33.f;
                mixerYellow  = 33.f;
                mixerMagenta = 33.f;
                mixerPurple  = 33.f;
                mixerCyan    = 33.f;
            }

            float filcor;
            Color::computeBWMixerConstants (params->blackwhite.setting, params->blackwhite.filter, params->blackwhite.algo, filcor,
                                            bwr, bwg, bwb, mixerOrange, mixerYellow, mixerCyan, mixerPurple, mixerMagenta,
                                            params->blackwhite.autoc, complem, kcorec, rrm, ggm, bbm);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 16)
#endif

            for (int i = 0; i < tH; i++) {
                for (int j = 0; j < tW; j++) {

                    //mix channel
                    tmpImage->r (i, j) = tmpImage->g (i, j) = tmpImage->b (i, j) = CLIP ((bwr * tmpImage->r (i, j) + bwg * tmpImage->g (i, j) + bwb * tmpImage->b (i, j)) * kcorec);

#ifndef __SSE2__

                    //gamma correction: pseudo TRC curve
                    if (hasgammabw) {
                        Color::trcGammaBW (tmpImage->r (i, j), tmpImage->g (i, j), tmpImage->b (i, j), gammabwr, gammabwg, gammabwb);
                    }

#endif
                }

#ifdef __SSE2__

                if (hasgammabw) {
                    //gamma correction: pseudo TRC curve
                    Color::trcGammaBWRow (tmpImage->r (i), tmpImage->g (i), tmpImage->b (i), tW, gammabwr, gammabwg, gammabwb);
                }

#endif
            }
        }

        if (editID == EUID_BlackWhiteAfterCurve) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 5)
#endif

            for (int i = 0; i < tH; i++) {
                for (int j = 0; j < tW; j++) {
                    editWhatever->v (i, j) = Color::gamma2curve[tmpImage->r (i, j)] / 65535.f; // assuming that r=g=b
                }
            }
        }

        if (hasToneCurvebw2) {

            if (afterCurveMode == BlackWhiteParams::TC_MODE_STD_BW) { // Standard
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic, 5)
#endif

                for (int i = 0; i < tH; i++) {
                    for (int j = 0; j < tW; j++) {
                        const StandardToneCurvebw& userToneCurve = static_cast<const StandardToneCurvebw&> (customToneCurvebw2);
                        userToneCurve.Apply (tmpImage->r (i, j), tmpImage->g (i, j), tmpImage->b (i, j));
                    }
                }
            } else if (afterCurveMode == BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW) { // apply the curve to the rgb channels, weighted
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic, 5)
#endif

                for (int i = 0; i < tH; i++) { //for ulterior usage if bw data modified
                    for (int j = 0; j < tW; j++) {
                        const WeightedStdToneCurvebw& userToneCurve = static_cast<const WeightedStdToneCurvebw&> (customToneCurvebw2);

                        tmpImage->r (i, j) = CLIP<float> (tmpImage->r (i, j));
                        tmpImage->g (i, j) = CLIP<float> (tmpImage->g (i, j));
                        tmpImage->b (i, j) = CLIP<float> (tmpImage->b (i, j));

                        userToneCurve.Apply (tmpImage->r (i, j), tmpImage->g (i, j), tmpImage->b (i, j));
                    }
                }
            }
        }

        //colortoning with black and white
        if (hasColorToning) {
            if (params->colorToning.method == "Splitco") {
                /*
                                    #if 1
                                    for (int i=istart,ti=0; i<tH; i++,ti++) {
                                        for (int j=jstart,tj=0; j<tW; j++,tj++) {
                                            float r = rtemp[ti*TS+tj];
                                            float g = gtemp[ti*TS+tj];
                                            float b = btemp[ti*TS+tj];
                                            float ro,go,bo;
                                            labtoning (r, g, b, ro, go, bo, 1, 0, 1, satLimit, satLimitOpacity, ctColorCurve, ctOpacityCurve, clToningcurve, cl2Toningcurve, 0.f, 1.f, wp, wip);
                                            rtemp[ti*TS+tj] = CLIP(ro);//I used CLIP because there is a little bug in gamutLchonly that return 65536.ii intead of 65535 ==> crash
                                            gtemp[ti*TS+tj] = CLIP(go);
                                            btemp[ti*TS+tj] = CLIP(bo);
                                        }
                                    }
                #else
                */
                int preser = 0;

                if (params->colorToning.lumamode) {
                    preser = 1;
                }

                float reducac = 0.3f;
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic, 5)
#endif

                for (int i = 0; i < tH; i++) {
                    for (int j = 0; j < tW; j++) {
                        float r = tmpImage->r (i, j);
                        float g = tmpImage->g (i, j);
                        float b = tmpImage->b (i, j);

                        float lumbefore = 0.299f * r + 0.587f * g + 0.114f * b;

                        if (lumbefore < 65000.f  && lumbefore > 500.f) { //reduct artifacts for highlights an extrem shadows
                            float ro, go, bo;
                            int mode = 1;
                            toningsmh (r, g, b, ro, go, bo, RedLow, GreenLow, BlueLow, RedMed, GreenMed, BlueMed, RedHigh, GreenHigh, BlueHigh, reducac, mode, preser, strProtect);
                            float lumafter = 0.299f * ro + 0.587f * go + 0.114f * bo;
                            float preserv = 1.f;

                            if (preser == 1) {
                                preserv = lumbefore / lumafter;
                            }

                            ro *= preserv;
                            go *= preserv;
                            bo *= preserv;
                            ro = CLIP (ro);
                            go = CLIP (go);
                            bo = CLIP (bo);
                            tmpImage->r (i, j) = ro;
                            tmpImage->g (i, j) = go;
                            tmpImage->b (i, j) = bo;
                        }
                    }
                }

//#endif
            }

            else if (params->colorToning.method == "Splitlr") {
                float balanS, balanH;
                float reducac =  0.4f;
                int preser = 0;

                if (params->colorToning.lumamode) {
                    preser = 1;
                }

                balanS = 1.f + Balan / 100.f; //balan between 0 and 2
                balanH = 1.f - Balan / 100.f;
                float rh, gh, bh;
                float rl, gl, bl;
                float xh, yh, zh;
                float xl, yl, zl;
                float iplow, iphigh;
                iplow = (float)ctColorCurve.low;
                iphigh = (float)ctColorCurve.high;

                //2 colours
                ctColorCurve.getVal (iphigh, xh, yh, zh);
                ctColorCurve.getVal (iplow, xl, yl, zl);

                Color::xyz2rgb (xh, yh, zh, rh, gh, bh, wip);
                Color::xyz2rgb (xl, yl, zl, rl, gl, bl, wip);

                //retrieve rgb value with s and l =1
                retreavergb (rl, gl, bl);
                retreavergb (rh, gh, bh);
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic, 5)
#endif

                for (int i = 0; i < tH; i++) {
                    for (int j = 0; j < tW; j++) {
                        float r = tmpImage->r (i, j);
                        float g = tmpImage->g (i, j);
                        float b = tmpImage->b (i, j);

                        float ro, go, bo;
                        int mode = 1;
                        toning2col (r, g, b, ro, go, bo, iplow, iphigh, rl, gl, bl, rh, gh, bh, SatLow, SatHigh, balanS, balanH, reducac, mode, preser, strProtect);
                        tmpImage->r (i, j) = ro;
                        tmpImage->g (i, j) = go;
                        tmpImage->b (i, j) = bo;
                    }
                }
            }

            //colortoning with shift color Lab
            else if (params->colorToning.method == "Lab"  && opautili) {
                int algm = 0;
                bool twocol = true;
                int metchrom = 0;

                if      (params->colorToning.twocolor == "Std"  ) {
                    metchrom = 0;
                } else if (params->colorToning.twocolor == "All"  ) {
                    metchrom = 1;
                } else if (params->colorToning.twocolor == "Separ") {
                    metchrom = 2;
                } else if (params->colorToning.twocolor == "Two"  ) {
                    metchrom = 3;
                }

                if (metchrom == 3) {
                    twocol = false;
                }

                float iplow = 0.f, iphigh = 0.f;

                if (!twocol) {
                    iplow = (float)ctColorCurve.low;
                    iphigh = (float)ctColorCurve.high;

                }

                int twoc = 0; //integer instead of bool to let more possible choice...other than 2 and 500.

                if (!twocol) {
                    twoc = 0;    // 2 colours
                } else {
                    twoc = 1;    // 500 colours
                }

                if     (params->colorToning.method == "Lab") {
                    algm = 1;
                } else if (params->colorToning.method == "Lch") {
                    algm = 2;    //in case of
                }

                if (algm <= 2) {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic, 5)
#endif

                    for (int i = 0; i < tH; i++) {
                        for (int j = 0; j < tW; j++) {
                            float r = tmpImage->r (i, j);
                            float g = tmpImage->g (i, j);
                            float b = tmpImage->b (i, j);
                            float ro, bo, go;
                            labtoning (r, g, b, ro, go, bo, algm, metchrom,  twoc, satLimit, satLimitOpacity, ctColorCurve,  ctOpacityCurve, clToningcurve, cl2Toningcurve,  iplow, iphigh,  wp,  wip);
                            tmpImage->r (i, j) = CLIP (ro);
                            tmpImage->g (i, j) = CLIP (go);
                            tmpImage->b (i, j) = CLIP (bo);

                        }
                    }
                }
            }

            else if (params->colorToning.method.substr (0, 3) == "RGB"  && opautili) {
                // color toning
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic, 5)
#endif

                for (int i = 0; i < tH; i++) {
                    for (int j = 0; j < tW; j++) {
                        float r = tmpImage->r (i, j);
                        float g = tmpImage->g (i, j);
                        float b = tmpImage->b (i, j);

                        // Luminance = (0.299f*r + 0.587f*g + 0.114f*b)

                        float h, s, l;
                        Color::rgb2hsl (r, g, b, h, s, l);

                        float l_ = Color::gamma_srgb (l * 65535.f) / 65535.f;

                        // get the opacity and tweak it to preserve saturated colors
                        float opacity = ctOpacityCurve.lutOpacityCurve[l_ * 500.f] / 4.f;

                        float r2, g2, b2;
                        ctColorCurve.getVal (l_, r2, g2, b2); // get the color from the color curve

                        float h2, s2, l2;
                        Color::rgb2hsl (r2, g2, b2, h2, s2, l2); // transform this new color to hsl

                        Color::hsl2rgb (h2, s2, l, r2, g2, b2);

                        tmpImage->r (i, j) = r + (r2 - r) * opacity;
                        tmpImage->g (i, j) = g + (g2 - g) * opacity;
                        tmpImage->b (i, j) = b + (b2 - b) * opacity;
                    }
                }
            }
        }

        // filling the pipette buffer by the content of the temp pipette buffers
        // due to optimization, we have to test now if the pipette has been filled in the second tile loop, by
        // testing editID
        /*if (editImgFloat) {
            for (int i=istart,ti=0; i<tH; i++,ti++)
                for (int j=jstart,tj=0; j<tW; j++,tj++) {
                    editImgFloat->r(i,j) = editIFloatTmpR[ti*TS+tj];
                    editImgFloat->g(i,j) = editIFloatTmpG[ti*TS+tj];
                    editImgFloat->b(i,j) = editIFloatTmpB[ti*TS+tj];
                }
        }
        else*/
        /*
        if (editWhatever && (editID==EUID_BlackWhiteAfterCurve)) {
            for (int i=istart,ti=0; i<tH; i++,ti++)
                for (int j=jstart,tj=0; j<tW; j++,tj++) {
                    editWhatever->v(i,j) = editWhateverTmp[ti*TS+tj];
                }
        }
        */

        // ready, fill lab (has to be the same code than the "fill lab" above!)

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 5)
#endif

        for (int i = 0; i < tH; i++) {
            for (int j = 0; j < tW; j++) {
                float r = tmpImage->r (i, j);
                float g = tmpImage->g (i, j);
                float b = tmpImage->b (i, j);

                float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
                float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;

                float fx, fy, fz;

                fx = (x < MAXVALF ? Color::cachef[x] : 327.68f * std::cbrt (x / MAXVALF));
                fy = (y < MAXVALF ? Color::cachef[y] : 327.68f * std::cbrt (y / MAXVALF));
                fz = (z < MAXVALF ? Color::cachef[z] : 327.68f * std::cbrt (z / MAXVALF));

                lab->L[i][j] = 116.0f *  fy - 5242.88f; //5242.88=16.0*327.68;
                lab->a[i][j] = 500.0f * (fx - fy);
                lab->b[i][j] = 200.0f * (fy - fz);
            }
        }


    }

    if (tmpImage) {
        delete tmpImage;
    }

    if (hCurveEnabled) {
        delete hCurve;
    }

    if (sCurveEnabled) {
        delete sCurve;
    }

    if (vCurveEnabled) {
        delete vCurve;
    }
}

/**
* @brief retreave RGB value with maximum saturation
* @param r red input and in exit new r
* @param g green input and in exit new g
* @param b blue input and in exit new b
**/
void ImProcFunctions::retreavergb (float &r, float &g, float &b)
{
    float mini = min (r, g, b);
    float maxi = max (r, g, b);
    float kkm = 65535.f / maxi;

    if (b == mini && r == maxi) {
        r = 65535.f;
        g = kkm * (g - b);
        b = 0.f;
    } else if (b == mini && g == maxi) {
        g = 65535.f;
        r = kkm * (r - b);
        b = 0.f;
    } else if (g == mini && r == maxi) {
        r = 65535.f;
        b = kkm * (b - g);
        g = 0.f;
    } else if (g == mini && b == maxi) {
        b = 65535.f;
        r = kkm * (r - g);
        g = 0.f;
    } else if (r == mini && b == maxi) {
        b = 65535.f;
        g = kkm * (g - r);
        r = 0.f;
    } else if (r == mini && g == maxi) {
        g = 65535.f;
        b = kkm * (b - r);
        r = 0.f;
    }
}

/**
* @brief Interpolate by decreasing with a parabol k = aa*v*v + bb*v +c  v[0..1]
* @param reducac val ue of the reduction in the middle of the range
* @param vinf value [0..1] for beginning decrease
* @param aa second degree parameter
* @param bb first degree parameter
* @param cc third parameter
**/
void ImProcFunctions::secondeg_end (float reducac, float vinf, float &aa, float &bb, float &cc)
{
    float zrd = reducac; //value at me  linear =0.5
    float v0 = vinf; //max shadows
    float me = (1.f + v0) / 2.f; //"median" value = (v0 + 1.=/2)
    //float a1=1.f-v0;
    float a2 = me - v0;
    float a3 = 1.f - v0 * v0;
    float a4 = me * me - v0 * v0;
    aa = (1.f + (zrd - 1.f) * (1 - v0) / a2) / (a4 * (1.f - v0) / a2 - a3);
    bb = - (1.f + a3 * aa) / (1.f - v0);
    cc = - (aa + bb);
}

/**
* @brief Interpolate by increasing with a parabol k = aa*v*v + bb*v  v[0..1]
* @param reducac val ue of the reduction in the middle of the range
* @param vend value [0..1] for beginning increase
* @param aa second degree parameter
* @param bb first degree parameter
**/
void ImProcFunctions::secondeg_begin (float reducac, float vend, float &aam, float &bbm)
{
    float zrmd = reducac; //linear = 0.5
    float v0m = vend;
    float mem = vend / 2.f; //(0. + 0.8)/2.f
    aam = (1.f - zrmd * v0m / mem) / (v0m * v0m - mem * v0m); //
    bbm = (1.f - aam * v0m * v0m) / v0m;
}


/**
* @brief color toning with 9 sliders shadows middletones highlight
* @param r red input values [0..65535]
* @param g green input values [0..65535]
* @param b blue input values [0..65535]
* @param ro red output values [0..65535]
* @param go green output values [0..65535]
* @param bo blue output values [0..65535]
* @param RedLow    [0..1] value after transformations of sliders [-100..100] for shadows
* @param GreenLow  [0..1] value after transformations of sliders [-100..100] for shadows
* @param BlueLow   [0..1] value after transformations of sliders [-100..100] for shadows
* @param RedMed    [0..1] value after transformations of sliders [-100..100] for midtones
* @param GreenMed  [0..1] value after transformations of sliders [-100..100] for midtones
* @param BlueMed   [0..1] value after transformations of sliders [-100..100] for midtones
* @param RedHigh   [0..1] value after transformations of sliders [-100..100] for highlights
* @param GreenHigh [0..1] value after transformations of sliders [-100..100] for highlights
* @param BlueHigh  [0..1] value after transformations of sliders [-100..100] for highlights
* @param reducac value of the reduction in the middle of the range for second degree increse or decrease action
* @param mode ?
* @param preser whether to preserve luminance (if 1) or not
**/
void ImProcFunctions::toningsmh (float r, float g, float b, float &ro, float &go, float &bo, float RedLow, float GreenLow, float BlueLow, float RedMed, float GreenMed, float BlueMed, float RedHigh, float GreenHigh, float BlueHigh, float reducac, int mode, int preser, float strProtect)
{
    float bmu = mode == 1 ? 0.5f : 0.4f;
    float RedL   = 1.f + (RedLow  - 1.f) * 0.4f;
    float GreenL = 1.f + (GreenLow - 1.f) * 0.4f;
    float BlueL  = 1.f + (BlueLow - 1.f) * bmu;
    float h, s, v;
    Color::rgb2hsv (r, g, b, h, s, v);
    float ksat = 1.f;
    float ksatlow = 1.f;
//  float s_0=0.55f;
//  float s_1=0.85f;
    /*
    if(mode==0) {//color
        if(s < s_0) ksat=SQR((1.f/s_0)*s);
        if(s > s_1) ksat=SQR((1.f/(s_1-1.f))*s - (1.f/(s_1-1.f)));
    }
    */
    float kl = 1.f;
    float rlo = 1.f; //0.4  0.5
    float rlm = 1.5f; //1.1
    float rlh = 2.2f; //1.1
    float rlob = bmu; //for BW old mode

    if (mode == 0) { //color
        rlo *= pow_F (strProtect, 0.4f); //0.5 ==>  0.75
        rlh *= pow_F (strProtect, 0.4f);
        rlm *= pow_F (strProtect, 0.4f);
    } else { //bw coefficient to preserve same results as before for satlimtopacity = 0.5 (default)
        rlo = strProtect * 0.8f; //0.4
        rlob = strProtect; //0.5
        rlm = strProtect * 2.2f; //1.1
        rlh = strProtect * 2.4f; //1.2
    }

    if (mode == 0) {
        rlob = rlo;
    }


    //fixed value of reducac=0.3
    //secondeg_end (reducac, v0, aa, bb, cc);
    if (mode == 1) {
        reducac = 0.5f; //black and white mode

        if (v > 0.15f) {
            kl = (-1.f / 0.85f) * v + (1.f) / 0.85f;    //Low light ==> decrease action after v=0.15
        }
    } else { //color
        float v0 = 0.15f;
        //second degree
        float aa, bb, cc;
        secondeg_end (reducac, v0, aa, bb, cc);
        float aab, bbb;
        secondeg_begin (0.7f, v0, aab, bbb);

        if (v > v0) {
            kl = aa * v * v + bb * v + cc;    //verified ==> exact
        } else if (mode == 0) {
            kl = aab * v * v + bbb * v;    //ksatlow=ksat;
        }
    }

    if (RedLow != 1.f) {
        RedL = 1.f + (RedLow - 1.f) * kl * ksat * rlo; //0.4

        if (RedLow >= 1.f) {
            g -= 20000.f * (RedL - 1.f) * ksatlow;
            b -= 20000.f * (RedL - 1.f) * ksatlow;
        } else {
            r += 20000.f * (RedL - 1.f) * ksatlow;
        }

        r = CLIP (r);
        g = CLIP (g);
        b = CLIP (b);
    }

    if (GreenLow != 1.f) {
        GreenL = 1.f + (GreenLow - 1.f) * kl * ksat * rlo; //0.4

        if (GreenLow >= 1.f) {
            r -= 20000.f * (GreenL - 1.f) * ksatlow;
            b -= 20000.f * (GreenL - 1.f) * ksatlow;
        } else {
            g += 20000.f * (GreenL - 1.f) * ksatlow;
        }

        r = CLIP (r);
        b = CLIP (b);
        g = CLIP (g);
    }

    if (BlueLow != 1.f) {
        BlueL = 1.f + (BlueLow - 1.f) * kl * ksat * rlob;

        if (BlueLow >= 1.f) {
            r -= 20000.f * (BlueL - 1.f) * ksatlow;
            g -= 20000.f * (BlueL - 1.f) * ksatlow;
        } else {
            b += 20000.f * (BlueL - 1.f) * ksatlow;
        }

        r = CLIP (r);
        g = CLIP (g);
        b = CLIP (b);
    }

    // mid tones
    float km;
    float v0m = 0.5f; //max action

    if (v < v0m) {
        float aam, bbm;
        float vend = v0m;
        secondeg_begin (reducac, vend, aam, bbm);
        km = aam * v * v + bbm * v; //verification = good
    } else {
        float v0mm = 0.5f; //max
        float aamm, bbmm, ccmm;
        secondeg_end (reducac, v0mm, aamm, bbmm, ccmm);
        km = aamm * v * v + bbmm * v + ccmm; //verification good
    }

    float RedM = 1.f + (RedMed - 1.f) * rlm;

    if (RedMed != 1.f) {
        RedM = 1.f + (RedMed - 1.f) * km * rlm;

        if (RedMed >= 1.f) {
            r += 20000.f * (RedM - 1.f);
            g -= 10000.f * (RedM - 1.f);
            b -= 10000.f * (RedM - 1.f);
            r = CLIP (r);
            g = CLIP (g);
            b = CLIP (b);
        } else {
            r += 10000.f * (RedM - 1.f);
            g -= 20000.f * (RedM - 1.f);
            b -= 20000.f * (RedM - 1.f);
            r = CLIP (r);
            g = CLIP (g);
            b = CLIP (b);
        }
    }

    float GreenM = 1.f + (GreenMed - 1.f) * rlm;

    if (GreenMed != 1.f) {
        GreenM = 1.f + (GreenMed - 1.f) * km * rlm;

        if (GreenMed >= 1.f) {
            r -= 10000.f * (GreenM - 1.f);
            g += 20000.f * (GreenM - 1.f);
            b -= 10000.f * (GreenM - 1.f);
            r = CLIP (r);
            g = CLIP (g);
            b = CLIP (b);
        } else {
            r -= 20000.f * (GreenM - 1.f);
            g += 10000.f * (GreenM - 1.f);
            b -= 20000.f * (GreenM - 1.f);
            r = CLIP (r);
            g = CLIP (g);
            b = CLIP (b);
        }
    }

    float BlueM = 1.f + (BlueMed - 1.f) * rlm;

    if (BlueMed != 1.f) {
        BlueM = 1.f + (BlueMed - 1.f) * km * rlm;

        if (BlueMed >= 1.f) {
            r -= 10000.f * (BlueM - 1.f);
            g -= 10000.f * (BlueM - 1.f);
            b += 20000.f * (BlueM - 1.f);
            r = CLIP (r);
            g = CLIP (g);
            b = CLIP (b);
        } else {
            r -= 20000.f * (BlueM - 1.f);
            g -= 20000.f * (BlueM - 1.f);
            b += 10000.f * (BlueM - 1.f);
            r = CLIP (r);
            g = CLIP (g);
            b = CLIP (b);
        }
    }

    //high tones
    float kh;
    kh = 1.f;
    float v00 = 0.8f; //max action
    float aa0, bb0;
    secondeg_begin (reducac, v00, aa0, bb0);
//  float hmu=1.5f;
//  if(mode==1) hmu=1.2f;//for BW old mode

    if (v > v00) {
        kh = (-1.f / (1.f - v00)) * v + (1.f) / (1.f - v00);    //High tones
    } else {
        kh = aa0 * v * v + bb0 * v;    //verification = good
    }

    float RedH = 1.f + (RedHigh - 1.f) * rlh;
    float GreenH = 1.f + (GreenHigh - 1.f) * rlh;
    float BlueH = 1.f + (BlueHigh - 1.f) * rlh; //1.2

    if (RedHigh != 1.f) {
        RedH = 1.f + (RedHigh - 1.f) * kh * rlh; //1.2

        if (RedHigh >= 1.f) {
            r += 20000.f * (RedH - 1.f);
            r = CLIP (r);
        } else {
            g -= 20000.f * (RedH - 1.f);
            b -= 20000.f * (RedH - 1.f);
        }

        g = CLIP (g);
        b = CLIP (b);
    }

    if (GreenHigh != 1.f) {
        GreenH = 1.f + (GreenHigh - 1.f) * kh * rlh; //1.2

        if (GreenHigh >= 1.f) {
            g += 20000.f * (GreenH - 1.f);
            g = CLIP (g);
        } else {
            r -= 20000.f * (GreenH - 1.f);
            b -= 20000.f * (GreenH - 1.f);
        }

        r = CLIP (r);
        b = CLIP (b);
    }

    if (BlueHigh != 1.f) {
        BlueH = 1.f + (BlueHigh - 1.f) * kh * rlh; //1.2

        if (BlueHigh >= 1.f) {
            b += 20000.f * (BlueH - 1.f);
            b = CLIP (b);
        } else {
            r -= 20000.f * (BlueH - 1.f);
            g -= 20000.f * (BlueH - 1.f);
        }

        r = CLIP (r);
        g = CLIP (g);
    }

    ro = r;
    go = g;
    bo = b;
}

/**
* @brief color toning with 2 colors - 2 sliders saturation shadows and highlight and one balance
* @param r g b input values [0..65535]
* @param ro go bo output values [0..65535]
* @param iplow iphigh [0..1] from curve color - value of luminance shadows and highlights
* @param rl gl bl [0..65535] - color of reference shadow
* @param rh gh bh [0..65535] - color of reference highlight
* @param SatLow SatHigh [0..1] from sliders saturation shadows and highlight
* @param balanS [0..1] balance for shadows (one slider)
* @param balanH [0..1] balance for highlights (same slider than for balanS)
* @param reducac value of the reduction in the middle of the range for second degree, increase or decrease action
**/
void ImProcFunctions::toning2col (float r, float g, float b, float &ro, float &go, float &bo, float iplow, float iphigh, float rl, float gl, float bl, float rh, float gh, float bh, float SatLow, float SatHigh, float balanS, float balanH, float reducac, int mode, int preser, float strProtect)
{
    float lumbefore = 0.299f * r + 0.587f * g + 0.114f * b;
    float h, s, l;
    Color::rgb2hsl (r, g, b, h, s, l);
    float v;
    Color::rgb2hsv (r, g, b, h, s, v);
    float ksat = 1.f;
    float ksatlow = 1.f;
    /*
        if(mode==0) {//color
            if(s < s_0) ksat=SQR((1.f/s_0)*s);
            if(s > s_1) ksat=SQR((1.f/(s_1-1.f))*s - (1.f/(s_1-1.f)));
        }
        */
    float kl = 1.f;
    float rlo = 1.f;
    float rlh = 2.2f;
    rlo *= pow_F (strProtect, 0.4f); //0.5 ==> 0.75  transfered value for more action
    rlh *= pow_F (strProtect, 0.4f);
    //low tones
    //second degree
    float aa, bb, cc;
    //fixed value of reducac =0.4;
    secondeg_end (reducac, iplow, aa, bb, cc);
    float aab, bbb;

    secondeg_begin (0.7f, iplow, aab, bbb);

    if (v > iplow) {
        kl = aa * v * v + bb * v + cc;
    } else if (mode == 0) {
        kl = aab * v * v + bbb * v;
    }


    if (SatLow > 0.f) {
        //rl gl bl
        float krl = rl / (rl + gl + bl);
        float kgl = gl / (rl + gl + bl);
        float kbl = bl / (rl + gl + bl);
        float RedL, GreenL, BlueL;

        if (g < 20000.f || b < 20000.f || r < 20000.f) {
            float kmgb = min (r, g, b);   //I have tested ...0.85 compromise...
            kl *= pow ((kmgb / 20000.f), 0.85f);
        }

        RedL = 1.f + (SatLow * krl) * kl * ksat * rlo * balanS; //0.4

        if (krl > 0.f) {
            g -= 20000.f * (RedL - 1.f) * ksatlow;
            b -= 20000.f * (RedL - 1.f) * ksatlow;
        }

        g = CLIP (g);
        b = CLIP (b);

        GreenL = 1.f + (SatLow * kgl) * kl * ksat * rlo * balanS; //0.4

        if (kgl > 0.f) {
            r -= 20000.f * (GreenL - 1.f) * ksatlow;
            b -= 20000.f * (GreenL - 1.f) * ksatlow;
        }

        r = CLIP (r);
        b = CLIP (b);

        BlueL = 1.f + (SatLow * kbl) * kl * ksat * rlo * balanS; //0.4

        if (kbl > 0.f) {
            r -= 20000.f * (BlueL - 1.f) * ksatlow;
            g -= 20000.f * (BlueL - 1.f) * ksatlow;
        }

        r = CLIP (r);
        g = CLIP (g);
    }

    //high tones
    float kh = 1.f;
    float aa0, bb0;
    //fixed value of reducac ==0.4;
    secondeg_begin (reducac, iphigh, aa0, bb0);

    if (v > iphigh) {
        kh = (-1.f / (1.f - iphigh)) * v + (1.f) / (1.f - iphigh);    //Low light ==> decrease action after iplow
    } else {
        kh = aa0 * v * v + bb0 * v;
    }


    if (g > 45535.f || b > 45535.f || r > 45535.f) {
        float kmgb = max (r, g, b);
        float cora = 1.f / (45535.f - 65535.f);
        float corb = 1.f - cora * 45535.f;
        float cor = kmgb * cora + corb;
        kh *= cor;
        /* best algo if necessary with non linear response...little differences and more time!
        float aa=1.f /(pow(45535.f,0.65f) - pow(65535.f,0.65f));
        float bb=1.f-aa*pow(45535.f,0.65f);
        float cor=aa*pow(kmbg,0.65f)+bb;
        kh*=cor;*/
    }


    if (SatHigh > 0.f) {
        float RedH, GreenH, BlueH;
        float krh = rh / (rh + gh + bh);
        float kgh = gh / (rh + gh + bh);
        float kbh = bh / (rh + gh + bh);
        RedH = 1.f + (SatHigh * krh) * kh * rlh * balanH; //1.2

        if (krh > 0.f) {
            r += 20000.f * (RedH - 1.f);
            r = CLIP (r);
        }

        g = CLIP (g);
        b = CLIP (b);

        GreenH = 1.f + (SatHigh * kgh) * kh * rlh * balanH; //1.2

        if (kgh > 0.f) {
            g += 20000.f * (GreenH - 1.f);
            g = CLIP (g);
        }

        r = CLIP (r);
        b = CLIP (b);
        BlueH = 1.f + (SatHigh * kbh) * kh * rlh * balanH; //1.2

        if (kbh > 0.f) {
            b += 20000.f * (BlueH - 1.f);
            b = CLIP (b);
        }

        r = CLIP (r);
        g = CLIP (g);
    }

    float lumafter = 0.299f * r + 0.587f * g + 0.114f * b;
    float preserv = 1.f;

    if (preser == 1) {
        preserv = lumbefore / lumafter;
    }

    //float preserv=lumbefore/lumafter;
    ro = r;
    go = g;
    bo = b;
    ro *= preserv;
    go *= preserv;
    bo *= preserv;
    ro = CLIP (ro);
    go = CLIP (go);
    bo = CLIP (bo);
}

/**
* @brief color toning with interpolation in mode Lab
* @param r g b input values [0..65535]
* @param ro go bo output values [0..65535]
* @param algm  metchrom twoc - methods
* @param ctColorCurve curve 500 colors
* @param ctOpacityCurve curve standard 'ab'
* @param clToningcurve  curve special 'ab' and 'a'
* @param cl2Toningcurve curve special 'b'
* @param iplow iphigh [0..1] luminance
* @param wp wip 3x3 matrix and inverse conversion rgb XYZ
**/
void ImProcFunctions::labtoning (float r, float g, float b, float &ro, float &go, float &bo, int algm, int metchrom, int twoc, float satLimit, float satLimitOpacity, const ColorGradientCurve & ctColorCurve, const OpacityCurve & ctOpacityCurve, LUTf & clToningcurve, LUTf & cl2Toningcurve, float iplow, float iphigh, double wp[3][3], double wip[3][3]  )
{
    float realL;
    float h, s, l;
    Color::rgb2hsl (r, g, b, h, s, l);
    float x2, y2, z2;
    float xl, yl, zl;

    if (twoc != 1) {
        l = (Color::gammatab_13_2[     l * 65535.f]) / 65535.f; //to compensate L from Lab
        iphigh = (Color::gammatab_13_2[iphigh * 65535.f]) / 65535.f;
        iplow  = (Color::gammatab_13_2[ iplow * 65535.f]) / 65535.f;
    }

    if (twoc == 1) {
        ctColorCurve.getVal (l, x2, y2, z2);
    } else {
        ctColorCurve.getVal (iphigh, x2, y2, z2);
        ctColorCurve.getVal (iplow, xl, yl, zl);
    }

    realL = l;


    //float opacity = ctOpacityCurve.lutOpacityCurve[l*500.f];
    //if(params->blackwhite.enabled){satLimit=80.f;satLimitOpacity=30.f;}//force BW

    // get the opacity and tweak it to preserve saturated colors
    //float l_ = Color::gamma_srgb(l*65535.f)/65535.f;
    float opacity;
    opacity = (1.f - min<float> (s / satLimit, 1.f) * (1.f - satLimitOpacity)) * ctOpacityCurve.lutOpacityCurve[l * 500.f];
    float opacity2 = (1.f - min<float> (s / satLimit, 1.f) * (1.f - satLimitOpacity));

    //float ro, go, bo;
    bool chr = true;
    bool lum = false;
    float lm = l;
    float chromat, luma;

    if (clToningcurve[lm * 65535.f] / (lm * 65535.f) < 1.f) {
        chromat = (clToningcurve[ (lm) * 65535.f] / (lm * 65535.f)) - 1.f;  //special effect
    } else {
        chromat = 1.f - SQR (SQR ((lm * 65535.f) / clToningcurve[ (lm) * 65535.f])); //apply C=f(L) acts  on 'a' and 'b'
    }

    if (cl2Toningcurve[lm * 65535.f] / (lm * 65535.f) < 1.f) {
        luma = (cl2Toningcurve[ (lm) * 65535.f] / (lm * 65535.f)) - 1.f;  //special effect
    } else {
        luma = 1.f - SQR (SQR ((lm * 65535.f) / (cl2Toningcurve[ (lm) * 65535.f]))); //apply C2=f(L) acts only on 'b'
    }

    int todo = 1;

    if (algm == 1) {
        Color::interpolateRGBColor (realL, iplow, iphigh, algm, opacity, twoc, metchrom, chr, lum, chromat, luma, r, g, b, xl, yl, zl, x2, y2, z2, todo, wp, wip, ro, go, bo);
    } else {
        Color::interpolateRGBColor (realL, iplow, iphigh, algm, opacity2, twoc, metchrom, chr, lum, chromat, luma, r, g, b, xl, yl, zl, x2, y2, z2, todo, wp, wip, ro, go, bo);
    }
}


void ImProcFunctions::luminanceCurve (LabImage* lold, LabImage* lnew, LUTf & curve)
{

    int W = lold->W;
    int H = lold->H;

    #pragma omp parallel for if (multiThread)

    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            float Lin = lold->L[i][j];
            //if (Lin>0 && Lin<65535)
            lnew->L[i][j] = curve[Lin];
        }
}



SSEFUNCTION void ImProcFunctions::chromiLuminanceCurve (PipetteBuffer *pipetteBuffer, int pW, LabImage* lold, LabImage* lnew, LUTf & acurve, LUTf & bcurve, LUTf & satcurve, LUTf & lhskcurve, LUTf & clcurve, LUTf & curve, bool utili, bool autili, bool butili, bool ccutili, bool cclutili, bool clcutili, LUTu &histCCurve, LUTu &histLCurve)
{
    int W = lold->W;
    int H = lold->H;
    // lhskcurve.dump("lh_curve");
    //init Flatcurve for C=f(H)

    PlanarWhateverData<float>* editWhatever = nullptr;
    EditUniqueID editID = EUID_None;
    bool editPipette = false;

    if (pipetteBuffer) {
        editID = pipetteBuffer->getEditID();

        if (editID != EUID_None) {

            switch  (pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType()) {
                case (BT_IMAGEFLOAT):
                    break;

                case (BT_LABIMAGE):
                    break;

                case (BT_SINGLEPLANE_FLOAT):
                    editPipette = true;
                    editWhatever = pipetteBuffer->getSinglePlaneBuffer();
                    break;
            }
        }
    }

    FlatCurve* chCurve = nullptr;// curve C=f(H)
    bool chutili = false;

    if (params->labCurve.chromaticity > -100) {
        chCurve = new FlatCurve (params->labCurve.chcurve);

        if (!chCurve || chCurve->isIdentity()) {
            if (chCurve) {
                delete chCurve;
                chCurve = nullptr;
            }
        }//do not use "Munsell" if Chcurve not used
        else {
            chutili = true;
        }
    }

    FlatCurve* lhCurve = nullptr;//curve L=f(H)
    bool lhutili = false;

    if (params->labCurve.chromaticity > -100) {
        lhCurve = new FlatCurve (params->labCurve.lhcurve);

        if (!lhCurve || lhCurve->isIdentity()) {
            if (lhCurve) {
                delete lhCurve;
                lhCurve = nullptr;
            }
        }//do not use "Munsell" if Chcurve not used
        else {
            lhutili = true;
        }
    }

    FlatCurve* hhCurve = nullptr;//curve H=f(H)
    bool hhutili = false;

    if (params->labCurve.chromaticity > -100) {
        hhCurve = new FlatCurve (params->labCurve.hhcurve);

        if (!hhCurve || hhCurve->isIdentity()) {
            if (hhCurve) {
                delete hhCurve;
                hhCurve = nullptr;
            }
        }//do not use "Munsell" if Chcurve not used
        else {
            hhutili = true;
        }
    }


    LUTu hist16Clad;
    LUTu hist16Llad;

    //preparate for histograms CIECAM
    if (pW != 1) { //only with improccoordinator
        hist16Clad (65536);
        hist16Clad.clear();
        hist16Llad (65536);
        hist16Llad.clear();

    }

#ifdef _DEBUG
    MyTime t1e, t2e;
    t1e.set();
    // init variables to display Munsell corrections
    MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif

    float adjustr = 1.0f;

//  if(params->labCurve.avoidclip ){
    // parameter to adapt curve C=f(C) to gamut

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

    // reference to the params structure has to be done outside of the parallelization to avoid CPU cache problem
    const bool highlight = params->toneCurve.hrenabled; //Get the value if "highlight reconstruction" is activated
    const int chromaticity = params->labCurve.chromaticity;
    const float chromapro = (chromaticity + 100.0f) / 100.0f;
    const bool bwonly = params->blackwhite.enabled && !params->colorToning.enabled;
    bool bwq = false;
//  if(params->ppVersion > 300  && params->labCurve.chromaticity == - 100) bwq = true;
    // const bool bwToning = params->labCurve.chromaticity == - 100  /*|| params->blackwhite.method=="Ch" || params->blackwhite.enabled */ || bwonly;
    const bool bwToning = bwq  /*|| params->blackwhite.method=="Ch" || params->blackwhite.enabled */ || bwonly;
    //if(chromaticity==-100) chromaticity==-99;
    const bool LCredsk = params->labCurve.lcredsk;
    const bool ccut = ccutili;
    const bool clut = clcutili;
    const double rstprotection = 100. - params->labCurve.rstprotection; // Red and Skin Tones Protection
    // avoid color shift is disabled when bwToning is activated and enabled if gamut is true in colorappearanace
    const bool avoidColorShift = (params->labCurve.avoidcolorshift || (params->colorappearance.gamut && params->colorappearance.enabled)) && !bwToning ;
    const float protectRed = (float)settings->protectred;
    const double protectRedH = settings->protectredh;
    float protect_red, protect_redh;
    protect_red = protectRed;//default=60  chroma: one can put more or less if necessary...in 'option'  40...160

    if (protect_red < 20.0f) {
        protect_red = 20.0;    // avoid too low value
    }

    if (protect_red > 180.0f) {
        protect_red = 180.0;    // avoid too high value
    }

    protect_redh = float (protectRedH); //default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0

    if (protect_redh < 0.1f) {
        protect_redh = 0.1f;    //avoid divide by 0 and negatives values
    }

    if (protect_redh > 1.0f) {
        protect_redh = 1.0f;    //avoid too big values
    }

    float protect_redhcur = protectRedH;//default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1

    if (protect_redhcur < 0.1f) {
        protect_redhcur = 0.1f;    //avoid divide by 0 and negatives values:minimal protection for transition
    }

    if (protect_redhcur > 3.5f) {
        protect_redhcur = 3.5f;    //avoid too big values
    }

    //increase saturation after denoise : ...approximation
    float factnoise = 1.f;

    if (params->dirpyrDenoise.enabled) {
        factnoise = (1.f + params->dirpyrDenoise.chroma / 500.f); //levels=5
//      if(yyyy) factnoise=(1.f+params->dirpyrDenoise.chroma/100.f);//levels=7
    }

    const float scaleConst = 100.0f / 100.1f;


    const bool gamutLch = settings->gamutLch;
    const float amountchroma = (float) settings->amchroma;

    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.working);
    const double wip[3][3] = {
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params->icm.working);
    const double wp[3][3] = {
        {wprof[0][0], wprof[0][1], wprof[0][2]},
        {wprof[1][0], wprof[1][1], wprof[1][2]},
        {wprof[2][0], wprof[2][1], wprof[2][2]}
    };

#ifdef _DEBUG
    #pragma omp parallel default(shared) firstprivate(lold, lnew, MunsDebugInfo, pW) if (multiThread)
#else
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float HHBuffer[W] ALIGNED16;
        float CCBuffer[W] ALIGNED16;
#endif
        #pragma omp for schedule(dynamic, 16)

        for (int i = 0; i < H; i++) {
            if (avoidColorShift)

                // only if user activate Lab adjustments
                if (autili || butili || ccutili ||  cclutili || chutili || lhutili || hhutili || clcutili || utili || chromaticity) {
                    Color::LabGamutMunsell (lold->L[i], lold->a[i], lold->b[i], W, /*corMunsell*/true, /*lumaMuns*/false, params->toneCurve.hrenabled, /*gamut*/true, wip, multiThread);
                }

#ifdef __SSE2__

            // precalculate some values using SSE
            if (bwToning || (!autili && !butili)) {
                __m128 c327d68v = _mm_set1_ps (327.68f);
                __m128 av, bv;
                int k;

                for (k = 0; k < W - 3; k += 4) {
                    av = LVFU (lold->a[i][k]);
                    bv = LVFU (lold->b[i][k]);
                    STVF (HHBuffer[k], xatan2f (bv, av));
                    STVF (CCBuffer[k], _mm_sqrt_ps (SQRV (av) + SQRV (bv)) / c327d68v);
                }

                for (; k < W; k++) {
                    HHBuffer[k] = xatan2f (lold->b[i][k], lold->a[i][k]);
                    CCBuffer[k] = sqrt (SQR (lold->a[i][k]) + SQR (lold->b[i][k])) / 327.68f;
                }
            }

#endif // __SSE2__

            for (int j = 0; j < W; j++) {
                const float Lin = lold->L[i][j];
                float LL = Lin / 327.68f;
                float CC;
                float HH;
                float Chprov;
                float Chprov1;
                float memChprov;
                float2 sincosval;

                if (bwToning) { // this values will be also set when bwToning is false some lines down
#ifdef __SSE2__
                    // use precalculated values from above
                    HH = HHBuffer[j];
                    CC = CCBuffer[j];
#else
                    HH = xatan2f (lold->b[i][j], lold->a[i][j]);
                    CC = sqrt (SQR (lold->a[i][j]) + SQR (lold->b[i][j])) / 327.68f;
#endif

                    // According to mathematical laws we can get the sin and cos of HH by simple operations
                    if (CC == 0.0f) {
                        sincosval.y = 1.0f;
                        sincosval.x = 0.0f;
                    } else {
                        sincosval.y = lold->a[i][j] / (CC * 327.68f);
                        sincosval.x = lold->b[i][j] / (CC * 327.68f);
                    }

                    Chprov = CC;
                    Chprov1 = CC;
                    memChprov = Chprov;
                }

                if (editPipette && editID == EUID_Lab_LCurve) {
                    editWhatever->v (i, j) = LIM01<float> (Lin / 32768.0f);  // Lab L pipette
                }

                lnew->L[i][j] = curve[Lin];

                float Lprov1 = (lnew->L[i][j]) / 327.68f;

                if (editPipette) {
                    if (editID == EUID_Lab_aCurve) { // Lab a pipette
                        float chromapipa = lold->a[i][j] + (32768.f * 1.28f);
                        editWhatever->v (i, j) = LIM01<float> ((chromapipa) / (65536.f * 1.28f));
                    } else if (editID == EUID_Lab_bCurve) { //Lab b pipette
                        float chromapipb = lold->b[i][j] + (32768.f * 1.28f);
                        editWhatever->v (i, j) = LIM01<float> ((chromapipb) / (65536.f * 1.28f));
                    }
                }

                float atmp, btmp;

                atmp = lold->a[i][j];

                if (autili) {
                    atmp = acurve[atmp + 32768.0f] - 32768.0f;    // curves Lab a
                }

                btmp = lold->b[i][j];

                if (butili) {
                    btmp = bcurve[btmp + 32768.0f] - 32768.0f;    // curves Lab b
                }

                if (!bwToning) { //take into account modification of 'a' and 'b'
#ifdef __SSE2__
                    if (!autili && !butili) {
                        // use precalculated values from above
                        HH = HHBuffer[j];
                        CC = CCBuffer[j];
                    } else {
                        CC = sqrt (SQR (atmp) + SQR (btmp)) / 327.68f;
                        HH = xatan2f (btmp, atmp);
                    }

#else
                    CC = sqrt (SQR (atmp) + SQR (btmp)) / 327.68f;
                    HH = xatan2f (btmp, atmp);
#endif

                    // According to mathematical laws we can get the sin and cos of HH by simple operations
                    //float2  sincosval;
                    if (CC == 0.f) {
                        sincosval.y = 1.f;
                        sincosval.x = 0.f;
                    } else {
                        sincosval.y = atmp / (CC * 327.68f);
                        sincosval.x = btmp / (CC * 327.68f);
                    }

                    Chprov = CC;
                    Chprov1 = CC;
                    memChprov = Chprov;
                } // now new values of lold with 'a' and 'b'

                if (editPipette)
                    if (editID == EUID_Lab_LHCurve || editID == EUID_Lab_CHCurve || editID == EUID_Lab_HHCurve) {//H pipette
                        float valpar = Color::huelab_to_huehsv2 (HH);
                        editWhatever->v (i, j) = valpar;
                    }

                if (lhutili) {  // L=f(H)
                    const float ClipLevel = 65535.f;
                    float l_r;//Luminance Lab in 0..1
                    l_r = Lprov1 / 100.f;
                    {
                        float valparam = float ((lhCurve->getVal (Color::huelab_to_huehsv2 (HH)) - 0.5f)); //get l_r=f(H)
                        float valparamneg;
                        valparamneg = valparam;
                        float kcc = (CC / amountchroma); //take Chroma into account...40 "middle low" of chromaticity (arbitrary and simple), one can imagine other algorithme
                        //reduce action for low chroma and increase action for high chroma
                        valparam *= 2.f * kcc;
                        valparamneg *= kcc; //slightly different for negative

                        if (valparam > 0.f) {
                            l_r = (1.f - valparam) * l_r + valparam * (1.f - SQR (((SQR (1.f - min (l_r, 1.0f))))));
                        } else
                            //for negative
                        {
                            float khue = 1.9f; //in reserve in case of!
                            l_r *= (1.f + khue * valparamneg);
                        }
                    }

                    Lprov1 = l_r * 100.f;

                    float Chprov2 = sqrt (SQR (atmp) + SQR (btmp)) / 327.68f;
                    //Gamut control especialy fot negative values slightly different of gamutlchonly
                    bool inRGB;

                    do {
                        inRGB = true;
                        float aprov1 = Chprov2 * sincosval.y;
                        float bprov1 = Chprov2 * sincosval.x;

                        float fy = (0.00862069f * Lprov1 ) + 0.137932f;
                        float fx = (0.002f * aprov1) + fy;
                        float fz = fy - (0.005f * bprov1);

                        float x_ = 65535.0f * Color::f2xyz (fx) * Color::D50x;
                        float z_ = 65535.0f * Color::f2xyz (fz) * Color::D50z;
                        float y_ = (Lprov1 > Color::epskap) ? 65535.0 * fy * fy * fy : 65535.0 * Lprov1 / Color::kappa;
                        float R, G, B;
                        Color::xyz2rgb (x_, y_, z_, R, G, B, wip);

                        if (R < 0.0f || G < 0.0f || B < 0.0f) {
                            if (Lprov1 < 0.1f) {
                                Lprov1 = 0.1f;
                            }

                            Chprov2 *= 0.95f;
                            inRGB = false;
                        } else if (!highlight && (R > ClipLevel || G > ClipLevel || B > ClipLevel)) {
                            if (Lprov1 > 99.98f) {
                                Lprov1 = 99.98f;
                            }

                            Chprov2 *= 0.95f;
                            inRGB = false;
                        }
                    } while (!inRGB);

                    atmp = 327.68f * Chprov2 * sincosval.y;
                    btmp = 327.68f * Chprov2 * sincosval.x;
                }

//          calculate C=f(H)
                if (chutili) {
                    double hr = Color::huelab_to_huehsv2 (HH);
                    float chparam = float ((chCurve->getVal (hr) - 0.5f) * 2.0f); //get C=f(H)
                    float chromaChfactor = 1.0f + chparam;
                    atmp *= chromaChfactor;//apply C=f(H)
                    btmp *= chromaChfactor;
                }

                if (hhutili) {  // H=f(H)
                    //hue Lab in -PI +PI
                    float valparam = float ((hhCurve->getVal (Color::huelab_to_huehsv2 (HH)) - 0.5f) * 1.7f) + HH; //get H=f(H)  1.7 optimisation !
                    HH = valparam;
                    sincosval = xsincosf (HH);
                }

                if (!bwToning) {
                    float factorskin, factorsat, factorskinext;

                    if (chromapro > 1.f) {
                        float scale = scaleConst;//reduction in normal zone
                        float scaleext = 1.f;//reduction in transition zone
                        Color::scalered ( rstprotection, chromapro, 0.0, HH, protect_redh, scale, scaleext);//1.0
                        float interm = (chromapro - 1.f);
                        factorskin = 1.f + (interm * scale);
                        factorskinext = 1.f + (interm * scaleext);
                    } else {
                        factorskin = chromapro ; // +(chromapro)*scale;
                        factorskinext = chromapro ;// +(chromapro)*scaleext;
                    }

                    factorsat = chromapro * factnoise;

                    //simulate very approximative gamut f(L) : with pyramid transition
                    float dred /*=55.f*/;//C red value limit

                    if     (Lprov1 < 25.f) {
                        dred = 40.f;
                    } else if (Lprov1 < 30.f) {
                        dred = 3.f * Lprov1 - 35.f;
                    } else if (Lprov1 < 70.f) {
                        dred = 55.f;
                    } else if (Lprov1 < 75.f) {
                        dred = -3.f * Lprov1 + 265.f;
                    } else {
                        dred = 40.f;
                    }

                    // end pyramid

                    // Test if chroma is in the normal range first
                    Color::transitred ( HH, Chprov1, dred, factorskin, protect_red, factorskinext, protect_redh, factorsat, factorsat);
                    atmp *= factorsat;
                    btmp *= factorsat;

                    if (editPipette && editID == EUID_Lab_CLCurve) {
                        editWhatever->v (i, j) = LIM01<float> (LL / 100.f);  // Lab C=f(L) pipette
                    }

                    if (clut) { // begin C=f(L)
                        float factorskin, factorsat, factor, factorskinext;
                        float chromaCfactor = (clcurve[LL * 655.35f]) / (LL * 655.35f); //apply C=f(L)
                        float curf = 0.7f; //empirical coeff because curve is more progressive
                        float scale = 100.0f / 100.1f; //reduction in normal zone for curve C
                        float scaleext = 1.0f; //reduction in transition zone for curve C
                        float protect_redcur, protect_redhcur; //perhaps the same value than protect_red and protect_redh
                        float deltaHH;//HH value transition for C curve
                        protect_redcur = curf * protectRed; //default=60  chroma: one can put more or less if necessary...in 'option'  40...160==> curf =because curve is more progressive

                        if (protect_redcur < 20.0f) {
                            protect_redcur = 20.0;    // avoid too low value
                        }

                        if (protect_redcur > 180.0f) {
                            protect_redcur = 180.0;    // avoid too high value
                        }

                        protect_redhcur = curf * float (protectRedH); //default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0 ==> curf =because curve is more progressive

                        if (protect_redhcur < 0.1f) {
                            protect_redhcur = 0.1f;    //avoid divide by 0 and negatives values
                        }

                        if (protect_redhcur > 1.0f) {
                            protect_redhcur = 1.0f;    //avoid too big values
                        }

                        deltaHH = protect_redhcur; //transition hue

                        if (chromaCfactor > 0.0) {
                            Color::scalered ( rstprotection, chromaCfactor, 0.0, HH, deltaHH, scale, scaleext);    //1.0
                        }

                        if (chromaCfactor > 1.0) {
                            float interm = (chromaCfactor - 1.0f) * 100.0f;
                            factorskin = 1.0f + (interm * scale) / 100.0f;
                            factorskinext = 1.0f + (interm * scaleext) / 100.0f;
                        } else {
                            factorskin = chromaCfactor; // +(1.0f-chromaCfactor)*scale;
                            factorskinext = chromaCfactor ; //+(1.0f-chromaCfactor)*scaleext;
                        }

                        factorsat = chromaCfactor;
                        factor = factorsat;
                        Color::transitred ( HH, Chprov1, dred, factorskin, protect_redcur, factorskinext, deltaHH, factorsat, factor);
                        atmp *= factor;
                        btmp *= factor;
                    }

                    // end C=f(L)
                    //  if (editID == EUID_Lab_CLCurve)
                    //      editWhatever->v(i,j) = LIM01<float>(Lprov2/100.f);// Lab C=f(L) pipette

                    // I have placed C=f(C) after all C treatments to assure maximum amplitude of "C"
                    if (editPipette && editID == EUID_Lab_CCurve) {
                        float chromapip = sqrt (SQR (atmp) + SQR (btmp) + 0.001f);
                        editWhatever->v (i, j) = LIM01<float> ((chromapip) / (48000.f));
                    }//Lab C=f(C) pipette

                    if (ccut) {
                        float factorskin, factorsat, factor, factorskinext;
                        float chroma = sqrt (SQR (atmp) + SQR (btmp) + 0.001f);
                        float chromaCfactor = (satcurve[chroma * adjustr]) / (chroma * adjustr); //apply C=f(C)
                        float curf = 0.7f; //empirical coeff because curve is more progressive
                        float scale = 100.0f / 100.1f; //reduction in normal zone for curve CC
                        float scaleext = 1.0f; //reduction in transition zone for curve CC
                        float protect_redcur, protect_redhcur; //perhaps the same value than protect_red and protect_redh
                        float deltaHH;//HH value transition for CC curve
                        protect_redcur = curf * protectRed; //default=60  chroma: one can put more or less if necessary...in 'option'  40...160==> curf =because curve is more progressive

                        if (protect_redcur < 20.0f) {
                            protect_redcur = 20.0;    // avoid too low value
                        }

                        if (protect_redcur > 180.0f) {
                            protect_redcur = 180.0;    // avoid too high value
                        }

                        protect_redhcur = curf * float (protectRedH); //default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0 ==> curf =because curve is more progressive

                        if (protect_redhcur < 0.1f) {
                            protect_redhcur = 0.1f;    //avoid divide by 0 and negatives values
                        }

                        if (protect_redhcur > 1.0f) {
                            protect_redhcur = 1.0f;    //avoid too big values
                        }

                        deltaHH = protect_redhcur; //transition hue

                        if (chromaCfactor > 0.0) {
                            Color::scalered ( rstprotection, chromaCfactor, 0.0, HH, deltaHH, scale, scaleext);    //1.0
                        }

                        if (chromaCfactor > 1.0) {
                            float interm = (chromaCfactor - 1.0f) * 100.0f;
                            factorskin = 1.0f + (interm * scale) / 100.0f;
                            factorskinext = 1.0f + (interm * scaleext) / 100.0f;
                        } else {
                            //factorskin= chromaCfactor*scale;
                            //factorskinext=chromaCfactor*scaleext;
                            factorskin = chromaCfactor; // +(1.0f-chromaCfactor)*scale;
                            factorskinext = chromaCfactor ; //+(1.0f-chromaCfactor)*scaleext;

                        }

                        factorsat = chromaCfactor;
                        factor = factorsat;
                        Color::transitred ( HH, Chprov1, dred, factorskin, protect_redcur, factorskinext, deltaHH, factorsat, factor);
                        atmp *= factor;
                        btmp *= factor;
                    }
                }

                // end chroma C=f(C)

                //update histogram C
                if (pW != 1) { //only with improccoordinator
                    int posp = (int)sqrt (atmp * atmp + btmp * btmp);
                    hist16Clad[posp]++;
                }

                if (editPipette && editID == EUID_Lab_LCCurve) {
                    float chromapiplc = sqrt (SQR (atmp) + SQR (btmp) + 0.001f);
                    editWhatever->v (i, j) = LIM01<float> ((chromapiplc) / (48000.f));
                }//Lab L=f(C) pipette


                if (cclutili && !bwToning) {    //apply curve L=f(C) for skin and rd...but also for extended color ==> near green and blue (see 'curf')

                    const float xx = 0.25f; //soft : between 0.2 and 0.4
                    float skdeltaHH;

                    skdeltaHH = protect_redhcur; //transition hue

                    float skbeg = -0.05f; //begin hue skin
                    float skend = 1.60f; //end hue skin
                    const float chrmin = 50.0f; //to avoid artifact, because L curve is not a real curve for luminance
                    float aa, bb;
                    float zz = 0.0f;
                    float yy = 0.0f;

                    if (Chprov1 < chrmin) {
                        yy = SQR (Chprov1 / chrmin) * xx;
                    } else {
                        yy = xx;    //avoid artifact for low C
                    }

                    if (!LCredsk) {
                        skbeg = -3.1415;
                        skend = 3.14159;
                        skdeltaHH = 0.001f;
                    }

                    if (HH > skbeg && HH < skend ) {
                        zz = yy;
                    } else if (HH > skbeg - skdeltaHH && HH <= skbeg) { //transition
                        aa = yy / skdeltaHH;
                        bb = -aa * (skbeg - skdeltaHH);
                        zz = aa * HH + bb;
                    } else if (HH >= skend && HH < skend + skdeltaHH) { //transition
                        aa = -yy / skdeltaHH;
                        bb = -aa * (skend + skdeltaHH);
                        zz = aa * HH + bb;
                    }

                    float chroma = sqrt (SQR (atmp) + SQR (btmp) + 0.001f);
                    float Lc = (lhskcurve[chroma * adjustr]) / (chroma * adjustr); //apply L=f(C)
                    Lc = (Lc - 1.0f) * zz + 1.0f; //reduct action
                    Lprov1 *= Lc; //adjust luminance
                }

                //update histo LC
                if (pW != 1) { //only with improccoordinator
                    int posl = Lprov1 * 327.68f;
                    hist16Llad[posl]++;
                }

                Chprov1 = sqrt (SQR (atmp) + SQR (btmp)) / 327.68f;

                // labCurve.bwtoning option allows to decouple modulation of a & b curves by saturation
                // with bwtoning enabled the net effect of a & b curves is visible
                if (bwToning) {
                    atmp -= lold->a[i][j];
                    btmp -= lold->b[i][j];
                }

                if (avoidColorShift) {
                    //gamutmap Lch ==> preserve Hue,but a little slower than gamutbdy for high values...and little faster for low values
                    if (gamutLch) {
                        float R, G, B;

#ifdef _DEBUG
                        bool neg = false;
                        bool more_rgb = false;
                        //gamut control : Lab values are in gamut
                        Color::gamutLchonly (HH, sincosval, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                        //gamut control : Lab values are in gamut
                        Color::gamutLchonly (HH, sincosval, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif
                        lnew->L[i][j] = Lprov1 * 327.68f;
//                  float2 sincosval = xsincosf(HH);
                        lnew->a[i][j] = 327.68f * Chprov1 * sincosval.y;
                        lnew->b[i][j] = 327.68f * Chprov1 * sincosval.x;
                    } else {
                        //use gamutbdy
                        //Luv limiter
                        float Y, u, v;
                        Color::Lab2Yuv (lnew->L[i][j], atmp, btmp, Y, u, v);
                        //Yuv2Lab includes gamut restriction map
                        Color::Yuv2Lab (Y, u, v, lnew->L[i][j], lnew->a[i][j], lnew->b[i][j], wp);
                    }

                    if (utili || autili || butili || ccut || clut || cclutili || chutili || lhutili || hhutili || clcutili || chromaticity) {
                        float correctionHue = 0.f; // Munsell's correction
                        float correctlum = 0.f;

                        Lprov1 = lnew->L[i][j] / 327.68f;
                        Chprov = sqrt (SQR (lnew->a[i][j]) + SQR (lnew->b[i][j])) / 327.68f;

#ifdef _DEBUG
                        Color::AllMunsellLch (/*lumaMuns*/true, Lprov1, LL, HH, Chprov, memChprov, correctionHue, correctlum, MunsDebugInfo);
#else
                        Color::AllMunsellLch (/*lumaMuns*/true, Lprov1, LL, HH, Chprov, memChprov, correctionHue, correctlum);
#endif

                        if (correctionHue != 0.f || correctlum != 0.f) {
                            if (fabs (correctionHue) < 0.015f) {
                                HH += correctlum;    // correct only if correct Munsell chroma very little.
                            }

                            /*      if((HH>0.0f && HH < 1.6f)   && memChprov < 70.0f) HH+=correctlum;//skin correct
                                    else if(fabs(correctionHue) < 0.3f) HH+=0.08f*correctlum;
                                    else if(fabs(correctionHue) < 0.2f) HH+=0.25f*correctlum;
                                    else if(fabs(correctionHue) < 0.1f) HH+=0.35f*correctlum;
                                    else if(fabs(correctionHue) < 0.015f) HH+=correctlum;   // correct only if correct Munsell chroma very little.
                            */
                            sincosval = xsincosf (HH + correctionHue);
                        }

                        lnew->a[i][j] = 327.68f * Chprov * sincosval.y; // apply Munsell
                        lnew->b[i][j] = 327.68f * Chprov * sincosval.x;
                    }
                } else {
//              if(Lprov1 > maxlp) maxlp=Lprov1;
//              if(Lprov1 < minlp) minlp=Lprov1;
                    if (!bwToning) {
                        lnew->L[i][j] = Lprov1 * 327.68f;
//                  float2 sincosval = xsincosf(HH);
                        lnew->a[i][j] = 327.68f * Chprov1 * sincosval.y;
                        lnew->b[i][j] = 327.68f * Chprov1 * sincosval.x;
                    } else {
                        //Luv limiter only
                        lnew->a[i][j] = atmp;
                        lnew->b[i][j] = btmp;
                    }
                }
            }
        }
    } // end of parallelization

    if (pW != 1) { //only with improccoordinator
        //update histogram C  with data chromaticity and not with CC curve
        hist16Clad.compressTo (histCCurve);
        //update histogram L with data luminance
        hist16Llad.compressTo (histLCurve);
    }

#ifdef _DEBUG

    if (settings->verbose) {
        t2e.set();
        printf ("Color::AllMunsellLch (correction performed in %d usec):\n", t2e.etime (t1e));
        printf ("   Munsell chrominance: MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%u\n", MunsDebugInfo->maxdhue[0],    MunsDebugInfo->maxdhue[1],    MunsDebugInfo->maxdhue[2],    MunsDebugInfo->maxdhue[3],    MunsDebugInfo->depass);
        printf ("   Munsell luminance  : MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%u\n", MunsDebugInfo->maxdhuelum[0], MunsDebugInfo->maxdhuelum[1], MunsDebugInfo->maxdhuelum[2], MunsDebugInfo->maxdhuelum[3], MunsDebugInfo->depassLum);
    }

    delete MunsDebugInfo;
#endif

    if (chCurve) {
        delete chCurve;
    }

    if (lhCurve) {
        delete lhCurve;
    }

    if (hhCurve) {
        delete hhCurve;
    }

    //  t2e.set();
    //  printf("Chromil took %d nsec\n",t2e.etime(t1e));
}


//#include "cubic.cc"

void ImProcFunctions::colorCurve (LabImage* lold, LabImage* lnew)
{

    /*    LUT<double> cmultiplier(181021);

        double boost_a = ((float)params->colorBoost.amount + 100.0) / 100.0;
        double boost_b = ((float)params->colorBoost.amount + 100.0) / 100.0;

        double c, amul = 1.0, bmul = 1.0;
        if (boost_a > boost_b) {
            c = boost_a;
            if (boost_a > 0)
                bmul = boost_b / boost_a;
        }
        else {
            c = boost_b;
            if (boost_b > 0)
                amul = boost_a / boost_b;
        }

        if (params->colorBoost.enable_saturationlimiter && c>1.0) {
            // re-generate color multiplier lookup table
            double d = params->colorBoost.saturationlimit / 3.0;
            double alpha = 0.5;
            double threshold1 = alpha * d;
            double threshold2 = c*d*(alpha+1.0) - d;
            for (int i=0; i<=181020; i++) { // lookup table stores multipliers with a 0.25 chrominance resolution
                double chrominance = (double)i/4.0;
                if (chrominance < threshold1)
                    cmultiplier[i] = c;
                else if (chrominance < d)
                    cmultiplier[i] = (c / (2.0*d*(alpha-1.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
                else if (chrominance < threshold2)
                    cmultiplier[i] = (1.0 / (2.0*d*(c*(alpha+1.0)-2.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
                else
                    cmultiplier[i] = 1.0;
            }
        }

        float eps = 0.001;
        double shift_a = params->colorShift.a + eps, shift_b = params->colorShift.b + eps;

        float** oa = lold->a;
        float** ob = lold->b;

        #pragma omp parallel for if (multiThread)
        for (int i=0; i<lold->H; i++)
            for (int j=0; j<lold->W; j++) {

                double wanted_c = c;
                if (params->colorBoost.enable_saturationlimiter && c>1) {
                    float chroma = (float)(4.0 * sqrt((oa[i][j]+shift_a)*(oa[i][j]+shift_a) + (ob[i][j]+shift_b)*(ob[i][j]+shift_b)));
                    wanted_c = cmultiplier [chroma];
                }

                double real_c = wanted_c;
                if (wanted_c >= 1.0 && params->colorBoost.avoidclip) {
                    double cclip = 100000.0;
                    double cr = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, 3.079935, -1.5371515, -0.54278342);
                    double cg = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, -0.92123418, 1.87599, 0.04524418);
                    double cb = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, 0.052889682, -0.20404134, 1.15115166);
                    if (cr>1.0 && cr<cclip) cclip = cr;
                    if (cg>1.0 && cg<cclip) cclip = cg;
                    if (cb>1.0 && cb<cclip) cclip = cb;
                    if (cclip<100000.0) {
                        real_c = -cclip + 2.0*cclip / (1.0+exp(-2.0*wanted_c/cclip));
                        if (real_c<1.0)
                            real_c = 1.0;
                    }
                }

                float nna = ((oa[i][j]+shift_a) * real_c * amul);
                float nnb = ((ob[i][j]+shift_b) * real_c * bmul);
                lnew->a[i][j] = LIM(nna,-32000.0f,32000.0f);
                lnew->b[i][j] = LIM(nnb,-32000.0f,32000.0f);
            }
    */
    //delete [] cmultiplier;
}

void ImProcFunctions::impulsedenoise (LabImage* lab)
{

    if (params->impulseDenoise.enabled && lab->W >= 8 && lab->H >= 8)

    {
        impulse_nr (lab, (float)params->impulseDenoise.thresh / 20.0 );
    }
}

void ImProcFunctions::impulsedenoisecam (CieImage* ncie, float **buffers[3])
{

    if (params->impulseDenoise.enabled && ncie->W >= 8 && ncie->H >= 8)

    {
        impulse_nrcam (ncie, (float)params->impulseDenoise.thresh / 20.0, buffers );
    }
}

void ImProcFunctions::defringe (LabImage* lab)
{

    if (params->defringe.enabled && lab->W >= 8 && lab->H >= 8)

    {
        PF_correct_RT (lab, lab, params->defringe.radius, params->defringe.threshold);
    }
}

void ImProcFunctions::defringecam (CieImage* ncie)
{
    if (params->defringe.enabled && ncie->W >= 8 && ncie->H >= 8) {
        PF_correct_RTcam (ncie, ncie, params->defringe.radius, params->defringe.threshold);
    }
}

void ImProcFunctions::badpixcam (CieImage* ncie, double rad, int thr, int mode, float b_l, float t_l, float t_r, float b_r, float skinprot, float chrom, int hotbad)
{
    if (ncie->W >= 8 && ncie->H >= 8) {
        Badpixelscam (ncie, ncie, rad, thr, mode, b_l, t_l, t_r, b_r, skinprot, chrom, hotbad);
    }
}

void ImProcFunctions::badpixlab (LabImage* lab, double rad, int thr, int mode, float b_l, float t_l, float t_r, float b_r, float skinprot, float chrom)
{
    if (lab->W >= 8 && lab->H >= 8) {
        BadpixelsLab (lab, lab, rad, thr, mode, b_l, t_l, t_r, b_r, skinprot, chrom);
    }
}

void ImProcFunctions::dirpyrequalizer (LabImage* lab, int scale)
{
    if (params->dirpyrequalizer.enabled && lab->W >= 8 && lab->H >= 8) {
        float b_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[0]) / 100.0f;
        float t_l = static_cast<float> (params->dirpyrequalizer.hueskin.value[1]) / 100.0f;
        float b_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[2]) / 100.0f;
        float t_r = static_cast<float> (params->dirpyrequalizer.hueskin.value[3]) / 100.0f;
        int choice = 0; //I have  not disabled this statement in case of ! always 0
        //      if     (params->dirpyrequalizer.algo=="FI") choice=0;
        //      else if(params->dirpyrequalizer.algo=="LA") choice=1;
        float artifact = (float) settings->artifact_cbdl;

        if (artifact > 6.f) {
            artifact = 6.f;
        }

        if (artifact < 0.f) {
            artifact = 1.f;
        }

        float chrom = 50.f;

        if (params->dirpyrequalizer.gamutlab) {
            ImProcFunctions::badpixlab (lab, artifact, 5, 3, b_l, t_l, t_r, b_r, params->dirpyrequalizer.skinprotect, chrom);    //for artifacts
        }

        //dirpyrLab_equalizer(lab, lab, params->dirpyrequalizer.mult);
        dirpyr_equalizer (lab->L, lab->L, lab->W, lab->H, lab->a, lab->b, lab->a, lab->b, params->dirpyrequalizer.mult, params->dirpyrequalizer.threshold, params->dirpyrequalizer.skinprotect, params->dirpyrequalizer.gamutlab,  b_l, t_l, t_r, b_r, choice, scale);
    }
}
void ImProcFunctions::EPDToneMapCIE (CieImage *ncie, float a_w, float c_, float w_h, int Wid, int Hei, int begh, int endh, float minQ, float maxQ, unsigned int Iterates, int skip)
{

    if (!params->epd.enabled) {
        return;
    }

    if (params->wavelet.enabled  && params->wavelet.tmrs != 0) {
        return;
    }

    float stren = params->epd.strength;
    float edgest = params->epd.edgeStopping;
    float sca = params->epd.scale;
    float gamm = params->epd.gamma;
    float rew = params->epd.reweightingIterates;
    float Qpro = ( 4.0 / c_)  * ( a_w + 4.0 ) ; //estimate Q max if J=100.0
    float *Qpr = ncie->Q_p[0];

    if (settings->verbose) {
        printf ("minQ=%f maxQ=%f  Qpro=%f\n", minQ, maxQ, Qpro);
    }

    if (maxQ > Qpro) {
        Qpro = maxQ;
    }

    EdgePreservingDecomposition epd (Wid, Hei);

    #pragma omp parallel for

    for (int i = 0; i < Hei; i++)
        for (int j = 0; j < Wid; j++) {
            ncie->Q_p[i][j] = gamm * ncie->Q_p[i][j] / (Qpro);
        }

    float Compression = expf (-stren);      //This modification turns numbers symmetric around 0 into exponents.
    float DetailBoost = stren;

    if (stren < 0.0f) {
        DetailBoost = 0.0f;    //Go with effect of exponent only if uncompressing.
    }

    //Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
    if (Iterates == 0) {
        Iterates = (unsigned int) (edgest * 15.0);
    }

    //Jacques Desmis : always Iterates=5 for compatibility images between preview and output

    epd.CompressDynamicRange (Qpr, sca / (float)skip, edgest, Compression, DetailBoost, Iterates, rew);

    //Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
    float s = (1.0f + 38.7889f) * powf (Compression, 1.5856f) / (1.0f + 38.7889f * powf (Compression, 1.5856f));
#ifndef _DEBUG
    #pragma omp parallel for schedule(dynamic,10)
#endif

    for (int i = 0; i < Hei; i++)
        for (int j = 0; j < Wid; j++) {
            ncie->Q_p[i][j] = (ncie->Q_p[i][j] * Qpro) / gamm;
            ncie->M_p[i][j] *= s;
        }

    /*
        float *Qpr2 = new float[Wid*((heir)+1)];

            for (int i=heir; i<Hei; i++)
                for (int j=0; j<Wid; j++) { Qpr2[(i-heir)*Wid+j]=ncie->Q_p[i][j];}
        if(minQ>0.0) minQ=0.0;//normaly minQ always > 0...
    //  EdgePreservingDecomposition epd = EdgePreservingDecomposition(Wid, Hei);
    //EdgePreservingDecomposition epd = EdgePreservingDecomposition(Wid, Hei/2);
        for(i = N2; i != N; i++)
    //  for(i = begh*Wid; i != N; i++)
            //Qpr[i] = (Qpr[i]-minQ)/(maxQ+1.0);
            Qpr2[i-N2] = (Qpr2[i-N2]-minQ)/(Qpro+1.0);

        float Compression2 = expf(-stren);      //This modification turns numbers symmetric around 0 into exponents.
        float DetailBoost2 = stren;
        if(stren < 0.0f) DetailBoost2 = 0.0f;   //Go with effect of exponent only if uncompressing.

        //Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
        if(Iterates == 0) Iterates = (unsigned int)(edgest*15.0);


        epd.CompressDynamicRange(Qpr2, sca/(float)skip, edgest, Compression2, DetailBoost2, Iterates, rew, Qpr2);

        //Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
         float s2 = (1.0f + 38.7889f)*powf(Compression, 1.5856f)/(1.0f + 38.7889f*powf(Compression, 1.5856f));
            for (int i=heir; i<Hei; i++)
        //  for (int i=begh; i<endh; i++)
                for (int j=0; j<Wid; j++) {
                ncie->Q_p[i][j]=Qpr2[(i-heir)*Wid+j]*Qpro + minQ;
            //  Qpr[i*Wid+j]=Qpr[i*Wid+j]*maxQ + minQ;
            //  ncie->J_p[i][j]=(100.0* Qpr[i*Wid+j]*Qpr[i*Wid+j]) /(w_h*w_h);

                ncie->M_p[i][j]*=s2;
            }
                    delete [] Qpr2;

    */
}


//Map tones by way of edge preserving decomposition. Is this the right way to include source?
//#include "EdgePreservingDecomposition.cc"
void ImProcFunctions::EPDToneMap (LabImage *lab, unsigned int Iterates, int skip)
{
    //Hasten access to the parameters.
//  EPDParams *p = (EPDParams *)(&params->epd);

    //Enabled? Leave now if not.
//  if(!p->enabled) return;
    if (!params->epd.enabled) {
        return;
    }

    if (params->wavelet.enabled  && params->wavelet.tmrs != 0) {
        return;
    }

    float stren = params->epd.strength;
    float edgest = params->epd.edgeStopping;
    float sca = params->epd.scale;
    float gamm = params->epd.gamma;
    float rew = params->epd.reweightingIterates;
    //Pointers to whole data and size of it.
    float *L = lab->L[0];
    float *a = lab->a[0];
    float *b = lab->b[0];
    size_t N = lab->W * lab->H;
    EdgePreservingDecomposition epd (lab->W, lab->H);

    //Due to the taking of logarithms, L must be nonnegative. Further, scale to 0 to 1 using nominal range of L, 0 to 15 bit.
    float minL = FLT_MAX;
    float maxL = 0.f;
    #pragma omp parallel
    {
        float lminL = FLT_MAX;
        float lmaxL = 0.f;
        #pragma omp for

        for (size_t i = 0; i < N; i++) {
            if (L[i] < lminL) {
                lminL = L[i];
            }

            if (L[i] > lmaxL) {
                lmaxL = L[i];
            }
        }

        #pragma omp critical
        {
            if (lminL < minL) {
                minL = lminL;
            }

            if (lmaxL > maxL) {
                maxL = lmaxL;
            }
        }
    }

    if (minL > 0.0f) {
        minL = 0.0f;    //Disable the shift if there are no negative numbers. I wish there were just no negative numbers to begin with.
    }

    #pragma omp parallel for

    for (size_t i = 0; i < N; ++i)
        //{L[i] = (L[i] - minL)/32767.0f;
    {
        L[i] = (L[i] - minL) / maxL;
        L[i] *= gamm;
    }

    //Some interpretations.
    float Compression = expf (-stren);      //This modification turns numbers symmetric around 0 into exponents.
    float DetailBoost = stren;

    if (stren < 0.0f) {
        DetailBoost = 0.0f;    //Go with effect of exponent only if uncompressing.
    }

    //Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
    if (Iterates == 0) {
        Iterates = (unsigned int) (edgest * 15.0f);
    }

    /* Debuggery. Saves L for toying with outside of RT.
    char nm[64];
    sprintf(nm, "%ux%ufloat.bin", lab->W, lab->H);
    FILE *f = fopen(nm, "wb");
    fwrite(L, N, sizeof(float), f);
    fclose(f);*/

    epd.CompressDynamicRange (L, sca / float (skip), edgest, Compression, DetailBoost, Iterates, rew);

    //Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
    float s = (1.0f + 38.7889f) * powf (Compression, 1.5856f) / (1.0f + 38.7889f * powf (Compression, 1.5856f));
#ifdef _OPENMP
    #pragma omp parallel for            // removed schedule(dynamic,10)
#endif

    for (size_t ii = 0; ii < N; ++ii) {
        a[ii] *= s;
        b[ii] *= s;
        L[ii] = L[ii] * maxL * (1.f / gamm) + minL;
    }
}


void ImProcFunctions::getAutoExp  (const LUTu &histogram, int histcompr, double defgain, double clip,
                                   double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh)
{

    float scale = 65536.0f;
    float midgray = 0.1842f; //middle gray in linear gamma =1 50% luminance

    int imax = 65536 >> histcompr;
    int overex = 0;
    float sum = 0.f, hisum = 0.f, losum = 0.f;
    float ave = 0.f, hidev = 0.f, lodev = 0.f;

    //find average luminance
    histogram.getSumAndAverage (sum, ave);

    //find median of luminance
    int median = 0, count = histogram[0];

    while (count < sum / 2) {
        median++;
        count += histogram[median];
    }

    if (median == 0 || ave < 1.f) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }

    // compute std dev on the high and low side of median
    // and octiles of histogram
    float octile[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}, ospread = 0.f;
    count = 0;

    int i = 0;

    for (; i < min ((int)ave, imax); i++) {
        if (count < 8) {
            octile[count] += histogram[i];

            if (octile[count] > sum / 8.f || (count == 7 && octile[count] > sum / 16.f)) {
                octile[count] = xlog (1. + (float)i) / log (2.f);
                count++;// = min(count+1,7);
            }
        }

        //lodev += SQR(ave-i)*histogram[i];
        lodev += (xlog (ave + 1.f) - xlog ((float)i + 1.)) * histogram[i];
        losum += histogram[i];
    }

    for (; i < imax; i++) {
        if (count < 8) {
            octile[count] += histogram[i];

            if (octile[count] > sum / 8.f || (count == 7 && octile[count] > sum / 16.f)) {
                octile[count] = xlog (1. + (float)i) / log (2.f);
                count++;// = min(count+1,7);
            }
        }

        //hidev += SQR(i-ave)*histogram[i];
        hidev += (xlog ((float)i + 1.) - xlog (ave + 1.f)) * histogram[i];
        hisum += histogram[i];

    }

    if (losum == 0 || hisum == 0) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }

//    lodev = (lodev / (log(2.f) * losum));
//    hidev = (hidev / (log(2.f) * hisum));

    if (octile[6] > log1p ((float)imax) / log2 (2.f)) { //if very overxposed image
        octile[6] = 1.5f * octile[5] - 0.5f * octile[4];
        overex = 2;
    }

    if (octile[7] > log1p ((float)imax) / log2 (2.f)) { //if overexposed
        octile[7] = 1.5f * octile[6] - 0.5f * octile[5];
        overex = 1;
    }

    // store values of octile[6] and octile[7] for calculation of exposure compensation
    // if we don't do this and the pixture is underexposed, calculation of exposure compensation assumes
    // that it's overexposed and calculates the wrong direction
    float oct6, oct7;
    oct6 = octile[6];
    oct7 = octile[7];


    for (int i = 1; i < 8; i++) {
        if (octile[i] == 0.0f) {
            octile[i] = octile[i - 1];
        }
    }

    // compute weighted average separation of octiles
    // for future use in contrast setting
    for (int i = 1; i < 6; i++) {
        ospread += (octile[i + 1] - octile[i]) / max (0.5f, (i > 2 ? (octile[i + 1] - octile[3]) : (octile[3] - octile[i])));
    }

    ospread /= 5.f;

    if (ospread <= 0.f) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }


    // compute clipping points based on the original histograms (linear, without exp comp.)
    unsigned int clipped = 0;
    int rawmax = (imax) - 1;

    while (histogram[rawmax] + clipped <= 0 && rawmax > 1) {
        clipped += histogram[rawmax];
        rawmax--;
    }

    //compute clipped white point
    unsigned int clippable = (int) (sum * clip / 100.f );
    clipped = 0;
    int whiteclip = (imax) - 1;

    while (whiteclip > 1 && (histogram[whiteclip] + clipped) <= clippable) {
        clipped += histogram[whiteclip];
        whiteclip--;
    }

    //compute clipped black point
    clipped = 0;
    int shc = 0;

    while (shc < whiteclip - 1 && histogram[shc] + clipped <= clippable) {
        clipped += histogram[shc];
        shc++;
    }

    //rescale to 65535 max
    rawmax <<= histcompr;
    whiteclip <<= histcompr;
    ave = ave * (1 << histcompr);
    median <<= histcompr;
    shc <<= histcompr;

//    //prevent division by 0
//    if (lodev == 0.f) {
//        lodev = 1.f;
//    }

    //compute exposure compensation as geometric mean of the amount that
    //sets the mean or median at middle gray, and the amount that sets the estimated top
    //of the histogram at or near clipping.
    //float expcomp1 = (log(/*(median/ave)*//*(hidev/lodev)*/midgray*scale/(ave-shc+midgray*shc))+log((hidev/lodev)))/log(2.f);
    float expcomp1 = (log (/*(median/ave)*//*(hidev/lodev)*/midgray * scale / (ave - shc + midgray * shc))) / log (2.f);
    float expcomp2;

    if (overex == 0) { // image is not overexposed
        expcomp2 = 0.5f * ( (15.5f - histcompr - (2.f * oct7 - oct6)) + log (scale / rawmax) / log (2.f) );
    } else {
        expcomp2 = 0.5f * ( (15.5f - histcompr - (2.f * octile[7] - octile[6])) + log (scale / rawmax) / log (2.f) );
    }

    if (fabs (expcomp1) - fabs (expcomp2) > 1.f) { //for great expcomp
        expcomp = (expcomp1 * fabs (expcomp2) + expcomp2 * fabs (expcomp1)) / (fabs (expcomp1) + fabs (expcomp2));
    } else {
        expcomp = 0.5 * (double)expcomp1 + 0.5 * (double) expcomp2; //for small expcomp
    }

    float gain = exp ((float)expcomp * log (2.f));

    float corr = sqrt (gain * scale / rawmax);
    black = (int) shc * corr;


    //now tune hlcompr to bring back rawmax to 65535
    hlcomprthresh = 33;
    //this is a series approximation of the actual formula for comp,
    //which is a transcendental equation
    float comp = (gain * ((float)whiteclip) / scale - 1.f) * 2.3f; // 2.3 instead of 2 to increase slightly comp
    hlcompr = (int) (100.*comp / (max (0.0, expcomp) + 1.0));
    hlcompr = max (0, min (100, hlcompr));

    //now find brightness if gain didn't bring ave to midgray using
    //the envelope of the actual 'control cage' brightness curve for simplicity
    float midtmp = gain * sqrt (median * ave) / scale;

    if (midtmp < 0.1f) {
        bright = (midgray - midtmp) * 15.0 / (midtmp);
    } else {
        bright = (midgray - midtmp) * 15.0 / (0.10833 - 0.0833 * midtmp);
    }

    bright = 0.25 */*(median/ave)*(hidev/lodev)*/max (0, bright);

    //compute contrast that spreads the average spacing of octiles
    contr = (int) 50.0f * (1.1f - ospread);
    contr = max (0, min (100, contr));
    //take gamma into account
    double whiteclipg = (int) (CurveFactory::gamma2 (whiteclip * corr / 65536.0) * 65536.0);

    float gavg = 0.;

    float val = 0.f;
    float increment = corr * (1 << histcompr);

    for (int i = 0; i < 65536 >> histcompr; i++) {
        gavg += histogram[i] * Color::gamma2curve[val];
        val += increment;
    }

    gavg /= sum;

    if (black < gavg) {
        int maxwhiteclip = (gavg - black) * 4 / 3 + black; // dont let whiteclip be such large that the histogram average goes above 3/4

        if (whiteclipg < maxwhiteclip) {
            whiteclipg = maxwhiteclip;
        }
    }

    whiteclipg = CurveFactory::igamma2 ((float) (whiteclipg / 65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter

    //corection with gamma
    black = (int) ((65535 * black) / whiteclipg);
    //expcomp = log(65535.0 / (whiteclipg)) / log(2.0);

    //diagnostics
    //printf ("**************** AUTO LEVELS ****************\n");
    /*
    if (settings->verbose) {
        printf("sum=%i clip=%f clippable=%i  clipWh=%i  clipBl=%i\n",somm, clip, clippable,clipwh, clipbl);
        printf ("expcomp1= %f   expcomp2= %f gain= %f  expcomp=%f\n",expcomp1,expcomp2,gain,expcomp);
        printf ("expo=%f\n",expo);
        printf ("median: %i  average: %f    median/average: %f\n",median,ave, median/ave);
        printf ("average: %f\n",ave);
        printf("comp=%f hlcomp: %i\n",comp, hlcompr);
        printf ("median/average: %f\n",median/ave);
        printf ("lodev: %f   hidev: %f      hidev/lodev: %f\n",lodev,hidev,hidev/lodev);
        printf ("lodev: %f\n",lodev);
        printf ("hidev: %f\n",hidev);
        printf ("imax=%d rawmax= %d  whiteclip= %d  gain= %f\n",imax,rawmax,whiteclip,gain);
        printf ("octiles: %f %f %f %f %f %f %f %f\n",octile[0],octile[1],octile[2],octile[3],octile[4],octile[5],octile[6],octile[7]);
        printf ("ospread= %f\n",ospread);
        printf ("overexp= %i\n",overex);
    }
    */
    /*
     // %%%%%%%%%% LEGACY AUTOEXPOSURE CODE %%%%%%%%%%%%%
     // black point selection is based on the linear result (yielding better visual results)
     black = (int)(shc * corr);
     // compute the white point of the exp. compensated gamma corrected image
     double whiteclipg = (int)(CurveFactory::gamma2 (whiteclip * corr / 65536.0) * 65536.0);

     // compute average intensity of the exp compensated, gamma corrected image
     double gavg = 0;
     for (int i=0; i<65536>>histcompr; i++)
     gavg += histogram[i] * CurveFactory::gamma2((int)(corr*(i<<histcompr)<65535 ? corr*(i<<histcompr) : 65535)) / sum;

     if (black < gavg) {
     int maxwhiteclip = (gavg - black) * 4 / 3 + black; // dont let whiteclip be such large that the histogram average goes above 3/4
     //double mavg = 65536.0 / (whiteclipg-black) * (gavg - black);
     if (whiteclipg < maxwhiteclip)
     whiteclipg = maxwhiteclip;
     }

     whiteclipg = CurveFactory::igamma2 ((float)(whiteclipg/65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter

     black = (int)((65535*black)/whiteclipg);
     expcomp = log(65535.0 / (whiteclipg)) / log(2.0);

     if (expcomp<0.0)   expcomp = 0.0;*/
    if (expcomp < -5.0) {
        expcomp = -5.0;
    }

    if (expcomp > 12.0) {
        expcomp = 12.0;
    }

    bright = max (-100, min (bright, 100));

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double ImProcFunctions::getAutoDistor  (const Glib::ustring &fname, int thumb_size)
{
    if (fname != "") {
        rtengine::RawMetaDataLocation ri;
        int w_raw = -1, h_raw = thumb_size;
        int w_thumb = -1, h_thumb = thumb_size;

        Thumbnail* thumb = rtengine::Thumbnail::loadQuickFromRaw (fname, ri, w_thumb, h_thumb, 1, FALSE);

        if (!thumb) {
            return 0.0;
        }

        Thumbnail* raw =   rtengine::Thumbnail::loadFromRaw      (fname, ri, w_raw, h_raw, 1, 1.0, FALSE, 0);

        if (!raw) {
            delete thumb;
            return 0.0;
        }

        if (h_thumb != h_raw) {
            delete thumb;
            delete raw;
            return 0.0;
        }

        int width;

        if (w_thumb > w_raw) {
            width = w_raw;
        } else {
            width = w_thumb;
        }

        unsigned char* thumbGray;
        unsigned char* rawGray;
        thumbGray = thumb->getGrayscaleHistEQ (width);
        rawGray = raw->getGrayscaleHistEQ (width);

        if (!thumbGray || !rawGray) {
            if (thumbGray) {
                delete thumbGray;
            }

            if (rawGray) {
                delete rawGray;
            }

            delete thumb;
            delete raw;
            return 0.0;
        }

        double dist_amount;
        int dist_result = calcDistortion (thumbGray, rawGray, width, h_thumb, 1, dist_amount);

        if (dist_result == -1) { // not enough features found, try increasing max. number of features by factor 4
            calcDistortion (thumbGray, rawGray, width, h_thumb, 4, dist_amount);
        }

        delete thumbGray;
        delete rawGray;
        delete thumb;
        delete raw;
        return dist_amount;
    } else {
        return 0.0;
    }
}

void ImProcFunctions::rgb2lab (const Imagefloat &src, LabImage &dst, const Glib::ustring &workingSpace)
{
    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix ( workingSpace );
    const float wp[3][3] = {
        {static_cast<float> (wprof[0][0]), static_cast<float> (wprof[0][1]), static_cast<float> (wprof[0][2])},
        {static_cast<float> (wprof[1][0]), static_cast<float> (wprof[1][1]), static_cast<float> (wprof[1][2])},
        {static_cast<float> (wprof[2][0]), static_cast<float> (wprof[2][1]), static_cast<float> (wprof[2][2])}
    };

    const int W = src.getWidth();
    const int H = src.getHeight();

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            float X, Y, Z;
            Color::rgbxyz (src.r (i, j), src.g (i, j), src.b (i, j), X, Y, Z, wp);
            //convert Lab
            Color::XYZ2Lab (X, Y, Z, dst.L[i][j], dst.a[i][j], dst.b[i][j]);
        }
    }
}

SSEFUNCTION void ImProcFunctions::lab2rgb (const LabImage &src, Imagefloat &dst, const Glib::ustring &workingSpace)
{
    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix ( workingSpace );
    const float wip[3][3] = {
        {static_cast<float> (wiprof[0][0]), static_cast<float> (wiprof[0][1]), static_cast<float> (wiprof[0][2])},
        {static_cast<float> (wiprof[1][0]), static_cast<float> (wiprof[1][1]), static_cast<float> (wiprof[1][2])},
        {static_cast<float> (wiprof[2][0]), static_cast<float> (wiprof[2][1]), static_cast<float> (wiprof[2][2])}
    };

    const int W = dst.getWidth();
    const int H = dst.getHeight();
#ifdef __SSE2__
    vfloat wipv[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            wipv[i][j] = F2V (wiprof[i][j]);
        }
    }

#endif

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int i = 0; i < H; i++) {
        int j = 0;
#ifdef __SSE2__

        for (; j < W - 3; j += 4) {
            vfloat X, Y, Z;
            vfloat R, G, B;
            Color::Lab2XYZ (LVFU (src.L[i][j]), LVFU (src.a[i][j]), LVFU (src.b[i][j]), X, Y, Z);
            Color::xyz2rgb (X, Y, Z, R, G, B, wipv);
            STVFU (dst.r (i, j), R);
            STVFU (dst.g (i, j), G);
            STVFU (dst.b (i, j), B);
        }

#endif

        for (; j < W; j++) {
            float X, Y, Z;
            Color::Lab2XYZ (src.L[i][j], src.a[i][j], src.b[i][j], X, Y, Z);
            Color::xyz2rgb (X, Y, Z, dst.r (i, j), dst.g (i, j), dst.b (i, j), wip);
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

}
