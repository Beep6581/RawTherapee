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
#include "rtengine.h"
#include "image8.h"
#include "imagefloat.h"
#include "labimage.h"
#include "improcfun.h"
#include <glibmm/ustring.h>
#include "iccstore.h"
#include "iccmatrices.h"
#include "settings.h"
#include "alignedbuffer.h"
#include "color.h"
#include "procparams.h"

namespace rtengine
{

extern void filmlike_clip(float *r, float *g, float *b);

namespace {

inline void copyAndClampLine(const float *src, unsigned char *dst, const int W)
{
    for (int j = 0; j < W * 3; ++j) {
        dst[j] = uint16ToUint8Rounded(CLIP(src[j] * MAXVALF));
    }
}


inline void copyAndClamp(const LabImage *src, unsigned char *dst, const double rgb_xyz[3][3], bool multiThread)
{
    int W = src->W;
    int H = src->H;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
    for (int i = 0; i < H; ++i) {
        float* rL = src->L[i];
        float* ra = src->a[i];
        float* rb = src->b[i];
        int ix = i * 3 * W;

        float R, G, B;
        float x_, y_, z_;

        for (int j = 0; j < W; ++j) {
            Color::Lab2XYZ(rL[j], ra[j], rb[j], x_, y_, z_ );
            Color::xyz2rgb(x_, y_, z_, R, G, B, rgb_xyz);

            dst[ix++] = uint16ToUint8Rounded(Color::gamma2curve[R]);
            dst[ix++] = uint16ToUint8Rounded(Color::gamma2curve[G]);
            dst[ix++] = uint16ToUint8Rounded(Color::gamma2curve[B]);
        }
    }
}

} // namespace

// Used in ImProcCoordinator::updatePreviewImage  (rtengine/improccoordinator.cc)
//         Crop::update                           (rtengine/dcrop.cc)
//         Thumbnail::processImage                (rtengine/rtthumbnail.cc)
//
// If monitorTransform, divide by 327.68 then apply monitorTransform (which can integrate soft-proofing)
// otherwise divide by 327.68, convert to xyz and apply the sRGB transform, before converting with gamma2curve
void ImProcFunctions::lab2monitorRgb(LabImage* lab, Image8* image)
{
    if (monitorTransform) {

        const int W = lab->W;
        const int H = lab->H;
        unsigned char * data = image->data;

        // cmsDoTransform is relatively expensive
#ifdef _OPENMP
        #pragma omp parallel firstprivate(lab, data, W, H)
#endif
        {
            AlignedBuffer<float> pBuf(3 * lab->W);

            AlignedBuffer<float> mBuf;
            AlignedBuffer<float> gwBuf1;
            AlignedBuffer<float> gwBuf2;

            if (gamutWarning) {
                gwBuf1.resize(3 * lab->W);
                gwBuf2.resize(3 * lab->W);
                mBuf.resize(3 * lab->W);
            }

            float *buffer = pBuf.data;
            float *outbuffer = gamutWarning ? mBuf.data : pBuf.data; // make in place transformations when gamutWarning is not needed

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int i = 0; i < H; i++) {

                const int ix = i * 3 * W;
                int iy = 0;

                float* rL = lab->L[i];
                float* ra = lab->a[i];
                float* rb = lab->b[i];

                for (int j = 0; j < W; j++) {
                    buffer[iy++] = rL[j] / 327.68f;
                    buffer[iy++] = ra[j] / 327.68f;
                    buffer[iy++] = rb[j] / 327.68f;
                }

                cmsDoTransform(monitorTransform, buffer, outbuffer, W);
                copyAndClampLine(outbuffer, data + ix, W);

                if (gamutWarning) {
                    gamutWarning->markLine(image, i, buffer, gwBuf1.data, gwBuf2.data);
                }
            }
        } // End of parallelization
    } else {
        copyAndClamp(lab, image->data, sRGB_xyz, multiThread);
    }
}



// Used in ImProcCoordinator::updatePreviewImage  (rtengine/improccoordinator.cc)
//         Crop::update                           (rtengine/dcrop.cc)
//
// Generate an Image8
//
// If output profile used, divide by 327.68 then apply the "profile" profile (eventually with a standard gamma)
// otherwise divide by 327.68, convert to xyz and apply the RGB transform, before converting with gamma2curve
Image8* ImProcFunctions::lab2rgb(LabImage* lab, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm, bool consider_histogram_settings)
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

    Image8* image = new Image8(cw, ch);
    Glib::ustring profile;

    bool standard_gamma;

    if (settings->HistogramWorking && consider_histogram_settings) {
        profile = icm.workingProfile;
        standard_gamma = true;
    } else {
        profile = icm.outputProfile;

        if (icm.outputProfile.empty() || icm.outputProfile == ColorManagementParams::NoICMString) {
            profile = "sRGB";
        }

        standard_gamma = false;
    }

    cmsHPROFILE oprof = ICCStore::getInstance()->getProfile(profile);

    if (oprof) {
        cmsHPROFILE oprofG = oprof;

        if (standard_gamma) {
            oprofG = ICCStore::makeStdGammaProfile(oprof);
        }

        cmsUInt32Number flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

        if (icm.outputBPC) {
            flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
        }

        lcmsMutex->lock();
        cmsHPROFILE LabIProf  = cmsCreateLab4Profile(nullptr);
        cmsHTRANSFORM hTransform = cmsCreateTransform (LabIProf, TYPE_Lab_DBL, oprofG, TYPE_RGB_FLT, icm.outputIntent, flags);  // NOCACHE is important for thread safety
        cmsCloseProfile(LabIProf);
        lcmsMutex->unlock();

        unsigned char *data = image->data;

        // cmsDoTransform is relatively expensive
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            AlignedBuffer<double> pBuf(3 * cw);
            AlignedBuffer<float> oBuf(3 * cw);
            double *buffer = pBuf.data;
            float *outbuffer = oBuf.data;
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

                cmsDoTransform (hTransform, buffer, outbuffer, cw);
                copyAndClampLine(outbuffer, data + ix, cw);
            }
        } // End of parallelization

        cmsDeleteTransform(hTransform);

        if (oprofG != oprof) {
            cmsCloseProfile(oprofG);
        }
    } else {
        const auto xyz_rgb = ICCStore::getInstance()->workingSpaceInverseMatrix(profile);
        copyAndClamp(lab, image->data, xyz_rgb, multiThread);
    }

    return image;
}


/** @brief Convert the final Lab image to the output RGB color space
 *
 * Used in processImage   (rtengine/simpleprocess.cc)
 *
 * Provide a pointer to a 7 floats array for "ga" (uninitialized ; this array will be filled with the gamma values) if you want
 * to use the custom gamma scenario. Those gamma values will correspond to the ones of the chosen standard output profile
 * (Prophoto if non standard output profile given)
 *
 * If "ga" is NULL, then we're considering standard gamma with the chosen output profile.
 *
 * Generate an Image16
 *
 * If a custom gamma profile can be created, divide by 327.68, convert to xyz and apply the custom gamma transform
 * otherwise divide by 327.68, convert to xyz and apply the sRGB transform, before converting with gamma2curve
 */
Imagefloat* ImProcFunctions::lab2rgbOut(LabImage* lab, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm)
{

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

    Imagefloat* image = new Imagefloat(cw, ch);
    cmsHPROFILE oprof = ICCStore::getInstance()->getProfile(icm.outputProfile);

    if (oprof) {
        cmsUInt32Number flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

        if (icm.outputBPC) {
            flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
        }

        lcmsMutex->lock();
        cmsHPROFILE iprof = cmsCreateLab4Profile(nullptr);
        cmsHTRANSFORM hTransform = cmsCreateTransform(iprof, TYPE_Lab_FLT, oprof, TYPE_RGB_FLT, icm.outputIntent, flags);
        lcmsMutex->unlock();

        image->ExecCMSTransform(hTransform, *lab, cx, cy);
        cmsDeleteTransform(hTransform);
        image->normalizeFloatTo65535();
    } else {
        
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

        for (int i = cy; i < cy + ch; i++) {
            float R, G, B;
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];

            for (int j = cx; j < cx + cw; j++) {

                float fy = (Color::c1By116 * rL[j]) / 327.68f + Color::c16By116; // (L+16)/116
                float fx = (0.002f * ra[j]) / 327.68f + fy;
                float fz = fy - (0.005f * rb[j]) / 327.68f;
                float LL = rL[j] / 327.68f;

                float x_ = 65535.0f * Color::f2xyz(fx) * Color::D50x;
                //float y_ = 65535.0 * Color::f2xyz(fy);
                float z_ = 65535.0f * Color::f2xyz(fz) * Color::D50z;
                float y_ = (LL > (float)Color::epskap) ? 65535.0f * fy * fy * fy : 65535.0f * LL / (float)Color::kappa;

                Color::xyz2srgb(x_, y_, z_, R, G, B);

                image->r(i - cy, j - cx) = Color::gamma2curve[CLIP(R)];
                image->g(i - cy, j - cx) = Color::gamma2curve[CLIP(G)];
                image->b(i - cy, j - cx) = Color::gamma2curve[CLIP(B)];
            }
        }
    }

    return image;
}


void ImProcFunctions::workingtrc(const Imagefloat* src, Imagefloat* dst, int cw, int ch, int mul, const Glib::ustring &profile, double gampos, double slpos, cmsHTRANSFORM &transform, bool normalizeIn, bool normalizeOut, bool keepTransForm) const
{
    const TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);

    const float toxyz[3][3] = {
        {
            static_cast<float>(wprof[0][0] / ((normalizeIn ? 65535.0 : 1.0))), //I have suppressed / Color::D50x
            static_cast<float>(wprof[0][1] / ((normalizeIn ? 65535.0 : 1.0))),
            static_cast<float>(wprof[0][2] / ((normalizeIn ? 65535.0 : 1.0)))
        }, {
            static_cast<float>(wprof[1][0] / (normalizeIn ? 65535.0 : 1.0)),
            static_cast<float>(wprof[1][1] / (normalizeIn ? 65535.0 : 1.0)),
            static_cast<float>(wprof[1][2] / (normalizeIn ? 65535.0 : 1.0))
        }, {
            static_cast<float>(wprof[2][0] / ((normalizeIn ? 65535.0 : 1.0))), //I have suppressed / Color::D50z
            static_cast<float>(wprof[2][1] / ((normalizeIn ? 65535.0 : 1.0))),
            static_cast<float>(wprof[2][2] / ((normalizeIn ? 65535.0 : 1.0)))
        }
    };

    cmsHTRANSFORM hTransform = nullptr;
    if (transform) {
        hTransform = transform;
    } else {
        double pwr = 1.0 / gampos;
        double ts = slpos;
        int five = mul;


        if (gampos < 1.0) {
            pwr = gampos;
            gampos = 1. / gampos;
            five = -mul;
        }

        //  int select_temp = 1; //5003K
        constexpr double eps = 0.000000001; // not divide by zero

        enum class ColorTemp {
            D50 = 5003, // for Widegamut, ProPhoto Best, Beta -> D50
            D65 = 6504, // for sRGB, AdobeRGB, Bruce Rec2020  -> D65
            D60 = 6005  // for ACES AP0 and AP1

        };
        ColorTemp temp = ColorTemp::D50;

        float p[6]; //primaries

        //primaries for 10 working profiles ==> output profiles
        if (profile == "WideGamut") {
            p[0] = 0.7350;    //Widegamut primaries
            p[1] = 0.2650;
            p[2] = 0.1150;
            p[3] = 0.8260;
            p[4] = 0.1570;
            p[5] = 0.0180;
        } else if (profile == "Adobe RGB") {
            p[0] = 0.6400;    //Adobe primaries
            p[1] = 0.3300;
            p[2] = 0.2100;
            p[3] = 0.7100;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (profile == "sRGB") {
            p[0] = 0.6400;    // sRGB primaries
            p[1] = 0.3300;
            p[2] = 0.3000;
            p[3] = 0.6000;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (profile == "BruceRGB") {
            p[0] = 0.6400;    // Bruce primaries
            p[1] = 0.3300;
            p[2] = 0.2800;
            p[3] = 0.6500;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (profile == "Beta RGB") {
            p[0] = 0.6888;    // Beta primaries
            p[1] = 0.3112;
            p[2] = 0.1986;
            p[3] = 0.7551;
            p[4] = 0.1265;
            p[5] = 0.0352;
        } else if (profile == "BestRGB") {
            p[0] = 0.7347;    // Best primaries
            p[1] = 0.2653;
            p[2] = 0.2150;
            p[3] = 0.7750;
            p[4] = 0.1300;
            p[5] = 0.0350;
        } else if (profile == "Rec2020") {
            p[0] = 0.7080;    // Rec2020 primaries
            p[1] = 0.2920;
            p[2] = 0.1700;
            p[3] = 0.7970;
            p[4] = 0.1310;
            p[5] = 0.0460;
            temp = ColorTemp::D65;
        } else if (profile == "ACESp0") {
            p[0] = 0.7347;    // ACES P0 primaries
            p[1] = 0.2653;
            p[2] = 0.0000;
            p[3] = 1.0;
            p[4] = 0.0001;
            p[5] = -0.0770;
            temp = ColorTemp::D60;
        } else if (profile == "ACESp1") {
            p[0] = 0.713;    // ACES P1 primaries
            p[1] = 0.293;
            p[2] = 0.165;
            p[3] = 0.830;
            p[4] = 0.128;
            p[5] = 0.044;
            temp = ColorTemp::D60;
        } else if (profile == "ProPhoto") {
            p[0] = 0.7347;    //ProPhoto and default primaries
            p[1] = 0.2653;
            p[2] = 0.1596;
            p[3] = 0.8404;
            p[4] = 0.0366;
            p[5] = 0.0001;
        } else {
            p[0] = 0.7347;    //default primaries always unused
            p[1] = 0.2653;
            p[2] = 0.1596;
            p[3] = 0.8404;
            p[4] = 0.0366;
            p[5] = 0.0001;
        }

        if (slpos == 0) {
            slpos = eps;
        }

        GammaValues g_a; //gamma parameters
        Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope : return parameters for LCMS2


        cmsFloat64Number gammaParams[7];
        gammaParams[4] = g_a[3] * ts;
        gammaParams[0] = gampos;
        gammaParams[1] = 1. / (1.0 + g_a[4]);
        gammaParams[2] = g_a[4] / (1.0 + g_a[4]);
        gammaParams[3] = 1. / slpos;
        gammaParams[5] = 0.0;
        gammaParams[6] = 0.0;
       // printf("ga0=%f ga1=%f ga2=%f ga3=%f ga4=%f\n", ga0, ga1, ga2, ga3, ga4);

        // 7 parameters for smoother curves
        cmsCIExyY xyD;
        cmsWhitePointFromTemp(&xyD, (double)temp);
        if (profile == "ACESp0") {
            xyD = {0.32168, 0.33767, 1.0};//refine white point to avoid differences
        }

        cmsToneCurve* GammaTRC[3];
        GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(NULL, five, gammaParams);//5 = more smoother than 4

        const cmsCIExyYTRIPLE Primaries = {
            {p[0], p[1], 1.0}, // red
            {p[2], p[3], 1.0}, // green
            {p[4], p[5], 1.0}  // blue
        };
        const cmsHPROFILE oprofdef = cmsCreateRGBProfile(&xyD, &Primaries, GammaTRC);
        cmsFreeToneCurve(GammaTRC[0]);

        if (oprofdef) {
            constexpr cmsUInt32Number flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;
            const cmsHPROFILE iprof = ICCStore::getInstance()->getXYZProfile();
            lcmsMutex->lock();
            hTransform = cmsCreateTransform(iprof, TYPE_RGB_FLT, oprofdef, TYPE_RGB_FLT, params->icm.outputIntent, flags);
            lcmsMutex->unlock();
        }
    }
    if (hTransform) {
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            AlignedBuffer<float> pBuf(cw * 3);
            const float normalize = normalizeOut ? 65535.f : 1.f;

#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 16) nowait
#endif

            for (int i = 0; i < ch; ++i) {
                float *p = pBuf.data;
                for (int j = 0; j < cw; ++j) {
                    const float r = src->r(i, j);
                    const float g = src->g(i, j);
                    const float b = src->b(i, j);

                    *(p++) = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
                    *(p++) = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                    *(p++) = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;
                }
                p = pBuf.data;
                cmsDoTransform(hTransform, p, p, cw);
                for (int j = 0; j < cw; ++j) {
                    dst->r(i, j) = *(p++) * normalize;
                    dst->g(i, j) = *(p++) * normalize;
                    dst->b(i, j) = *(p++) * normalize;
                }
            }
        }
        if (!keepTransForm) {
            cmsDeleteTransform(hTransform);
            hTransform = nullptr;
        }
        transform = hTransform;
    }
}


}
