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
#include "rtengine.h"
#include "improcfun.h"
#include <glibmm.h>
#include "iccstore.h"
#include "iccmatrices.h"
#include "../rtgui/options.h"
#include "settings.h"
#include "curves.h"
#include "alignedbuffer.h"
#include "color.h"

namespace rtengine
{

extern const Settings* settings;

void ImProcFunctions::lab2monitorRgb (LabImage* lab, Image8* image)
{
    if (monitorTransform) {

        int W = lab->W;
        int H = lab->H;
        unsigned char * data = image->data;

        // cmsDoTransform is relatively expensive
#ifdef _OPENMP
        #pragma omp parallel firstprivate(lab, data, W, H)
#endif
        {
            AlignedBuffer<float> pBuf(3 * lab->W);
            float *buffer = pBuf.data;

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int i = 0; i < H; i++) {

                const int ix = i * 3 * W;
                int iy = 0;

                float* rL = lab->L[i];
                float* ra = lab->a[i];
                float* rb = lab->b[i];

                float fy, fx, fz, x_, y_, z_, LL;

                for (int j = 0; j < W; j++) {
                    buffer[iy++] = rL[j] / 327.68f;
                    buffer[iy++] = ra[j] / 327.68f;
                    buffer[iy++] = rb[j] / 327.68f;
                }

                if (!settings->HistogramWorking && output2monitorTransform && lab2outputTransform) {
                    AlignedBuffer<float> buf(3 * W);
                    cmsDoTransform (lab2outputTransform, buffer, buf.data, W);
                    cmsDoTransform (output2monitorTransform, buf.data, data + ix, W);
                } else {
                    cmsDoTransform (monitorTransform, buffer, data + ix, W);
                }
            }

        } // End of parallelization

    } else {

        int W = lab->W;
        int H = lab->H;
        unsigned char * data = image->data;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

        for (int i = 0; i < H; i++) {
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
            int ix = i * 3 * W;

            float R, G, B;
            float fy, fx, fz, x_, y_, z_, LL;

            for (int j = 0; j < W; j++) {

                //float L1=rL[j],a1=ra[j],b1=rb[j];//for testing

                fy = (0.00862069 * rL[j]) / 327.68 + 0.137932; // (L+16)/116
                fx = (0.002 * ra[j]) / 327.68 + fy;
                fz = fy - (0.005 * rb[j]) / 327.68;
                LL = rL[j] / 327.68;

                x_ = 65535.0 * Color::f2xyz(fx) * Color::D50x;
                //  y_ = 65535.0 * Color::f2xyz(fy);
                z_ = 65535.0 * Color::f2xyz(fz) * Color::D50z;
                y_ = (LL > Color::epskap) ? 65535.0 * fy * fy * fy : 65535.0 * LL / Color::kappa;

                Color::xyz2srgb(x_, y_, z_, R, G, B);

                /* copy RGB */
                //int R1=((int)gamma2curve[(R)])
                data[ix++] = ((int)Color::gamma2curve[CLIP(R)]) >> 8;
                data[ix++] = ((int)Color::gamma2curve[CLIP(G)]) >> 8;
                data[ix++] = ((int)Color::gamma2curve[CLIP(B)]) >> 8;
            }
        }
    }
}

Image8* ImProcFunctions::lab2rgb (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile, RenderingIntent intent, bool standard_gamma)
{
    //gamutmap(lab);

    if (cx < 0) {
        cx = 0;
    }

    if (cy < 0) {
        cy = 0;
    }

    if (cx + cw > lab->W) {
        cw = lab->W - cx;
    }

    if (cy + ch > lab->H) {
        ch = lab->H - cy;
    }

    Image8* image = new Image8 (cw, ch);

    cmsHPROFILE oprof = iccStore->getProfile (profile);

    if (oprof) {
        cmsHPROFILE oprofG = oprof;

        if (standard_gamma) {
            oprofG = ICCStore::makeStdGammaProfile(oprof);
        }

        lcmsMutex->lock ();
        cmsHPROFILE hLab  = cmsCreateLab4Profile(NULL);
        cmsHTRANSFORM hTransform = cmsCreateTransform (hLab, TYPE_Lab_DBL, oprofG, TYPE_RGB_8, intent,
                                   cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety
        cmsCloseProfile(hLab);
        lcmsMutex->unlock ();

        unsigned char *data = image->data;

        // cmsDoTransform is relatively expensive
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            AlignedBuffer<double> pBuf(3 * cw);
            double *buffer = pBuf.data;
            int condition = cy + ch;

#ifdef _OPENMP
            #pragma omp for firstprivate(lab) schedule(dynamic,16)
#endif

            for (int i = cy; i < condition; i++) {
                const int ix = i * 3 * cw;
                int iy = 0;
                float* rL = lab->L[i];
                float* ra = lab->a[i];
                float* rb = lab->b[i];

                for (int j = cx; j < cx + cw; j++) {
                    buffer[iy++] = rL[j] / 327.68f;
                    buffer[iy++] = ra[j] / 327.68f;
                    buffer[iy++] = rb[j] / 327.68f;
                }

                cmsDoTransform (hTransform, buffer, data + ix, cw);
            }
        } // End of parallelization

        cmsDeleteTransform(hTransform);

        if (oprofG != oprof) {
            cmsCloseProfile(oprofG);
        }
    } else {

        const auto xyz_rgb = iccStore->workingSpaceInverseMatrix (profile);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

        for (int i = cy; i < cy + ch; i++) {
            float R, G, B;
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
            int ix = 3 * i * cw;

            for (int j = cx; j < cx + cw; j++) {

                float fy = (0.00862069 * rL[j]) / 327.68 + 0.137932; // (L+16)/116
                float fx = (0.002 * ra[j]) / 327.68 + fy;
                float fz = fy - (0.005 * rb[j]) / 327.68;
                float LL = rL[j] / 327.68;

                float x_ = 65535.0 * Color::f2xyz(fx) * Color::D50x;
                //float y_ = 65535.0 * Color::f2xyz(fy);
                float z_ = 65535.0 * Color::f2xyz(fz) * Color::D50z;
                float y_ = (LL > Color::epskap) ? 65535.0 * fy * fy * fy : 65535.0 * LL / Color::kappa;

                Color::xyz2rgb(x_, y_, z_, R, G, B, xyz_rgb);

                image->data[ix++] = (int)Color::gamma2curve[CLIP(R)] >> 8;
                image->data[ix++] = (int)Color::gamma2curve[CLIP(G)] >> 8;
                image->data[ix++] = (int)Color::gamma2curve[CLIP(B)] >> 8;
            }
        }
    }

    return image;
}
// for default (not gamma)
Image16* ImProcFunctions::lab2rgb16 (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile, RenderingIntent intent, bool bw)
{

    //gamutmap(lab);

    if (cx < 0) {
        cx = 0;
    }

    if (cy < 0) {
        cy = 0;
    }

    if (cx + cw > lab->W) {
        cw = lab->W - cx;
    }

    if (cy + ch > lab->H) {
        ch = lab->H - cy;
    }

    Image16* image = new Image16 (cw, ch);
    cmsHPROFILE oprof = iccStore->getProfile (profile);



    if (oprof) {
        #pragma omp parallel for if (multiThread)

        for (int i = cy; i < cy + ch; i++) {
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
            short* xa = (short*)image->r(i - cy);
            short* ya = (short*)image->g(i - cy);
            short* za = (short*)image->b(i - cy);

            for (int j = cx; j < cx + cw; j++) {

                float fy = (0.0086206897f * rL[j]) / 327.68f + 0.1379310345f; // (L+16)/116
                float fx = (0.002 * ra[j]) / 327.68f + fy;
                float fz = fy - (0.005f * rb[j]) / 327.68f;
                float LL = rL[j] / 327.68f;

                float x_ = 65535.0f * (float) Color::f2xyz(fx) * Color::D50x;
                //float y_ = 65535.0 * Color::f2xyz(fy);
                float z_ = 65535.0f * (float) Color::f2xyz(fz) * Color::D50z;
                float y_ = (LL > Color::epskap) ? 65535.0f * fy * fy * fy : 65535.0f * LL / Color::kappa;

                xa[j - cx] =  CLIP((int)  round(x_));
                ya[j - cx] =  CLIP((int)  round(y_));
                za[j - cx] = CLIP((int)   round(z_));

                if(bw && y_ < 65535.f ) { //force Bw value and take highlight into account
                    xa[j - cx] = (int) round(y_ * Color::D50x );
                    za[j - cx] = (int) round(y_ * Color::D50z);
                }

            }
        }

        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16, oprof, TYPE_RGB_16, intent, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
        lcmsMutex->unlock ();

        image->ExecCMSTransform(hTransform);

        cmsDeleteTransform(hTransform);
    } else {
        #pragma omp parallel for if (multiThread)

        for (int i = cy; i < cy + ch; i++) {
            float R, G, B;
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];

            for (int j = cx; j < cx + cw; j++) {

                float fy = (0.0086206897f * rL[j]) / 327.68f + 0.1379310345f; // (L+16)/116
                float fx = (0.002f * ra[j]) / 327.68f + fy;
                float fz = fy - (0.005f * rb[j]) / 327.68f;
                float LL = rL[j] / 327.68f;

                float x_ = 65535.0f * (float) Color::f2xyz(fx) * Color::D50x;
                //float y_ = 65535.0 * Color::f2xyz(fy);
                float z_ = 65535.0f * (float) Color::f2xyz(fz) * Color::D50z;
                float y_ = (LL > Color::epskap) ? (float) 65535.0f * fy * fy * fy : 65535.0f * LL / Color::kappa;

                Color::xyz2srgb(x_, y_, z_, R, G, B);

                image->r(i - cy, j - cx) = (int)Color::gamma2curve[CLIP(R)];
                image->g(i - cy, j - cx) = (int)Color::gamma2curve[CLIP(G)];
                image->b(i - cy, j - cx) = (int)Color::gamma2curve[CLIP(B)];
            }
        }
    }

    return image;
}


// for gamma options (BT709...sRGB linear...)
Image16* ImProcFunctions::lab2rgb16b (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile, RenderingIntent intent, Glib::ustring profi, Glib::ustring gam,  bool freegamma, double gampos, double slpos, double &ga0, double &ga1, double &ga2, double &ga3, double &ga4, double &ga5, double &ga6, bool bw)
{

    //gamutmap(lab);

    if (cx < 0) {
        cx = 0;
    }

    if (cy < 0) {
        cy = 0;
    }

    if (cx + cw > lab->W) {
        cw = lab->W - cx;
    }

    if (cy + ch > lab->H) {
        ch = lab->H - cy;
    }

    Image16* image = new Image16 (cw, ch);
    float p1, p2, p3, p4, p5, p6; //primaries

    double g_a0, g_a1, g_a2, g_a3, g_a4, g_a5; //gamma parameters
    double pwr;
    double ts;
    ga6 = 0.0;
    pwr = 1.0 / gampos;
    ts = slpos;
    int mode = 0, imax = 0;

    int t50;
    int select_temp = 1; //5003K
    const double eps = 0.000000001; // not divide by zero

    //primaries for 7 working profiles ==> output profiles
    // eventually to adapt primaries  if RT used special profiles !
    if (profi == "WideGamut") {
        p1 = 0.7350;    //Widegamut primaries
        p2 = 0.2650;
        p3 = 0.1150;
        p4 = 0.8260;
        p5 = 0.1570;
        p6 = 0.0180;
        select_temp = 1;
    } else if (profi == "Adobe RGB") {
        p1 = 0.6400;    //Adobe primaries
        p2 = 0.3300;
        p3 = 0.2100;
        p4 = 0.7100;
        p5 = 0.1500;
        p6 = 0.0600;
        select_temp = 2;
    } else if (profi == "sRGB") {
        p1 = 0.6400;    // sRGB primaries
        p2 = 0.3300;
        p3 = 0.3000;
        p4 = 0.6000;
        p5 = 0.1500;
        p6 = 0.0600;
        select_temp = 2;
    } else if (profi == "BruceRGB") {
        p1 = 0.6400;    // Bruce primaries
        p2 = 0.3300;
        p3 = 0.2800;
        p4 = 0.6500;
        p5 = 0.1500;
        p6 = 0.0600;
        select_temp = 2;
    } else if (profi == "Beta RGB") {
        p1 = 0.6888;    // Beta primaries
        p2 = 0.3112;
        p3 = 0.1986;
        p4 = 0.7551;
        p5 = 0.1265;
        p6 = 0.0352;
        select_temp = 1;
    } else if (profi == "BestRGB") {
        p1 = 0.7347;    // Best primaries
        p2 = 0.2653;
        p3 = 0.2150;
        p4 = 0.7750;
        p5 = 0.1300;
        p6 = 0.0350;
        select_temp = 1;
    } else if (profi == "Rec2020") {
        p1 = 0.7080;    // Rec2020 primaries
        p2 = 0.2920;
        p3 = 0.1700;
        p4 = 0.7970;
        p5 = 0.1310;
        p6 = 0.0460;
        select_temp = 2;
    } else {
        p1 = 0.7347;    //ProPhoto and default primaries
        p2 = 0.2653;
        p3 = 0.1596;
        p4 = 0.8404;
        p5 = 0.0366;
        p6 = 0.0001;
        select_temp = 1;
    }

    if (!freegamma) {//if Free gamma not selected
        // gamma : ga0,ga1,ga2,ga3,ga4,ga5 by calcul
        if(gam == "BT709_g2.2_s4.5")      {
            ga0 = 2.22;    //BT709  2.2  4.5  - my prefered as D.Coffin
            ga1 = 0.909995;
            ga2 = 0.090005;
            ga3 = 0.222222;
            ga4 = 0.081071;
            ga5 = 0.0;
        } else if (gam == "sRGB_g2.4_s12.92")   {
            ga0 = 2.40;    //sRGB 2.4 12.92  - RT default as Lightroom
            ga1 = 0.947858;
            ga2 = 0.052142;
            ga3 = 0.077399;
            ga4 = 0.039293;
            ga5 = 0.0;
        } else if (gam == "High_g1.3_s3.35")    {
            ga0 = 1.3 ;    //for high dynamic images
            ga1 = 0.998279;
            ga2 = 0.001721;
            ga3 = 0.298507;
            ga4 = 0.005746;
            ga5 = 0.0;
        } else if (gam == "Low_g2.6_s6.9")   {
            ga0 = 2.6 ;    //gamma 2.6 variable : for low contrast images
            ga1 = 0.891161;
            ga2 = 0.108839;
            ga3 = 0.144928;
            ga4 = 0.076332;
            ga5 = 0.0;
        } else if (gam == "linear_g1.0")   {
            ga0 = 1.0;    //gamma=1 linear : for high dynamic images (cf : D.Coffin...)
            ga1 = 1.;
            ga2 = 0.;
            ga3 = 1. / eps;
            ga4 = 0.;
            ga5 = 0.0;
        } else if (gam == "standard_g2.2")   {
            ga0 = 2.2;    //gamma=2.2 (as gamma of Adobe, Widegamut...)
            ga1 = 1.;
            ga2 = 0.;
            ga3 = 1. / eps;
            ga4 = 0.;
            ga5 = 0.0;
        } else if (gam == "standard_g1.8")   {
            ga0 = 1.8;    //gamma=1.8  (as gamma of Prophoto)
            ga1 = 1.;
            ga2 = 0.;
            ga3 = 1. / eps;
            ga4 = 0.;
            ga5 = 0.0;
        }
    } else { //free gamma selected
        if(slpos == 0) {
            slpos = eps;
        }

        Color::calcGamma(pwr, ts, mode, imax, g_a0, g_a1, g_a2, g_a3, g_a4, g_a5); // call to calcGamma with selected gamma and slope : return parameters for LCMS2
        ga4 = g_a3 * ts;
        //printf("g_a0=%f g_a1=%f g_a2=%f g_a3=%f g_a4=%f\n", g_a0,g_a1,g_a2,g_a3,g_a4);
        ga0 = gampos;
        ga1 = 1. / (1.0 + g_a4);
        ga2 = g_a4 / (1.0 + g_a4);
        ga3 = 1. / slpos;
        ga5 = 0.0;
        //printf("ga0=%f ga1=%f ga2=%f ga3=%f ga4=%f\n", ga0,ga1,ga2,ga3,ga4);

    }

    if(select_temp == 1) {
        t50 = 5003;    // for Widegamut, Prophoto Best, Beta   D50
    } else if (select_temp == 2) {
        t50 = 6504;    // for sRGB, AdobeRGB, Bruce Rec2020 D65
    }

    cmsCIExyY       xyD;
    cmsCIExyYTRIPLE Primaries = {{p1, p2, 1.0},//red primaries
        {p3, p4, 1.0}, // green
        {p5, p6, 1.0} //blue
    };
    cmsToneCurve* GammaTRC[3];
    cmsFloat64Number Parameters[7];
    Parameters[0] = ga0;
    Parameters[1] = ga1;
    Parameters[2] = ga2;
    Parameters[3] = ga3;
    Parameters[4] = ga4;
    Parameters[5] = ga5;
    Parameters[6] = ga6;
// 7 parameters for smoother curves
    cmsWhitePointFromTemp(&xyD, t50);
    GammaTRC[0] = GammaTRC[1] = GammaTRC[2] =   cmsBuildParametricToneCurve(NULL, 5, Parameters);//5 = more smoother than 4
    cmsHPROFILE oprofdef = cmsCreateRGBProfileTHR(NULL, &xyD, &Primaries, GammaTRC); //oprofdef  become Outputprofile

    cmsFreeToneCurve(GammaTRC[0]);


    if (oprofdef) {
        #pragma omp parallel for if (multiThread)

        for (int i = cy; i < cy + ch; i++) {
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
            short* xa = (short*)image->r(i - cy);
            short* ya = (short*)image->g(i - cy);
            short* za = (short*)image->b(i - cy);

            for (int j = cx; j < cx + cw; j++) {

                float fy = (0.0086206897f * rL[j]) / 327.68f + 0.1379310345f; // (L+16)/116
                float fx = (0.002f * ra[j]) / 327.68f + fy;
                float fz = fy - (0.005f * rb[j]) / 327.68f;
                float LL = rL[j] / 327.68f;

                float x_ = 65535.0f * (float)Color::f2xyz(fx) * Color::D50x;
                //  float y_ = 65535.0 * Color::f2xyz(fy);
                float z_ = 65535.0f * (float)Color::f2xyz(fz) * Color::D50z;
                float y_ = (LL > Color::epskap) ? (float) 65535.0 * fy * fy * fy : 65535.0f * LL / Color::kappa;

                xa[j - cx] = CLIP((int) round(x_)) ;
                ya[j - cx] = CLIP((int) round(y_));
                za[j - cx] = CLIP((int) round(z_));

                if(bw && y_ < 65535.f) { //force Bw value and take highlight into account
                    xa[j - cx] = (int) round(y_ * Color::D50x);
                    za[j - cx] = (int) round(y_ * Color::D50z);
                }

            }
        }

        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16, oprofdef, TYPE_RGB_16, intent,  cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
        lcmsMutex->unlock ();

        image->ExecCMSTransform(hTransform);
        cmsDeleteTransform(hTransform);
    } else {
        //
        #pragma omp parallel for if (multiThread)
        for (int i = cy; i < cy + ch; i++) {
            float R, G, B;
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];

            for (int j = cx; j < cx + cw; j++) {

                float fy = (0.0086206897f * rL[j]) / 327.68f + 0.1379310345f; // (L+16)/116
                float fx = (0.002f * ra[j]) / 327.68f + fy;
                float fz = fy - (0.005f * rb[j]) / 327.68f;
                float LL = rL[j] / 327.68f;

                float x_ = 65535.0f * (float) Color::f2xyz(fx) * Color::D50x;
                //float y_ = 65535.0 * Color::f2xyz(fy);
                float z_ = 65535.0f * (float) Color::f2xyz(fz) * Color::D50z;
                float y_ = (LL > Color::epskap) ? (float) 65535.0 * fy * fy * fy : 65535.0f * LL / Color::kappa;

                Color::xyz2srgb(x_, y_, z_, R, G, B);

                image->r(i - cy, j - cx) = (int)Color::gamma2curve[CLIP(R)];
                image->g(i - cy, j - cx) = (int)Color::gamma2curve[CLIP(G)];
                image->b(i - cy, j - cx) = (int)Color::gamma2curve[CLIP(B)];
            }
        }
    }

    return image;
}

//#include "sRGBgamutbdy.cc"

}
