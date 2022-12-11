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
#include <glibmm/ustring.h>

#include "alignedbuffer.h"
#include "color.h"
#include "iccmatrices.h"
#include "iccstore.h"
#include "image8.h"
#include "imagefloat.h"
#include "improcfun.h"
#include "labimage.h"
#include "procparams.h"
#include "rtengine.h"
#include "settings.h"
#include "utils.h"

namespace rtengine
{

namespace {

inline void copyAndClampLine(const float *src, unsigned char *dst, const int W)
{
    for (int j = 0; j < W * 3; ++j) {
        dst[j] = uint16ToUint8Rounded(CLIP(src[j] * MAXVALF));
    }
}


inline void copyAndClamp(const LabImage *src, unsigned char *dst, const double rgb_xyz[3][3], bool multiThread)
{
    const int W = src->W;
    const int H = src->H;

    float rgb_xyzf[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rgb_xyzf[i][j] = rgb_xyz[i][j];
        }
    }

#ifdef __SSE2__
    vfloat rgb_xyzv[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rgb_xyzv[i][j] = F2V(rgb_xyzf[i][j]);
        }
    }
#endif
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
    for (int i = 0; i < H; ++i) {
        float* rL = src->L[i];
        float* ra = src->a[i];
        float* rb = src->b[i];
        int ix = i * 3 * W;

#ifdef __SSE2__
        float rbuffer[W] ALIGNED16;
        float gbuffer[W] ALIGNED16;
        float bbuffer[W] ALIGNED16;
        int j = 0;
        for (; j < W - 3; j += 4) {
            vfloat R, G, B;
            vfloat x_, y_, z_;
            Color::Lab2XYZ(LVFU(rL[j]), LVFU(ra[j]), LVFU(rb[j]), x_, y_, z_ );
            Color::xyz2rgb(x_, y_, z_, R, G, B, rgb_xyzv);
            STVF(rbuffer[j], Color::gamma2curve[R]);
            STVF(gbuffer[j], Color::gamma2curve[G]);
            STVF(bbuffer[j], Color::gamma2curve[B]);
        }
        for (; j < W; ++j) {
            float R, G, B;
            float x_, y_, z_;
            Color::Lab2XYZ(rL[j], ra[j], rb[j], x_, y_, z_ );
            Color::xyz2rgb(x_, y_, z_, R, G, B, rgb_xyzf);
            rbuffer[j] = Color::gamma2curve[R];
            gbuffer[j] = Color::gamma2curve[G];
            bbuffer[j] = Color::gamma2curve[B];
        }

        for (j = 0; j < W; ++j) {
            dst[ix++] = uint16ToUint8Rounded(rbuffer[j]);
            dst[ix++] = uint16ToUint8Rounded(gbuffer[j]);
            dst[ix++] = uint16ToUint8Rounded(bbuffer[j]);
        }

#else
        for (int j = 0; j < W; ++j) {
            float R, G, B;
            float x_, y_, z_;
            Color::Lab2XYZ(rL[j], ra[j], rb[j], x_, y_, z_ );
            Color::xyz2rgb(x_, y_, z_, R, G, B, rgb_xyzf);

            dst[ix++] = uint16ToUint8Rounded(Color::gamma2curve[R]);
            dst[ix++] = uint16ToUint8Rounded(Color::gamma2curve[G]);
            dst[ix++] = uint16ToUint8Rounded(Color::gamma2curve[B]);
        }
#endif
    }
}

} // namespace


float gammalog(float x, float p, float s, float g3, float g4)
{
    return x <= g3 ? x * s : (1.f + g4) * xexpf(xlogf(x) / p) - g4;//continuous
}

#ifdef __SSE2__
vfloat gammalog(vfloat x, vfloat p, vfloat s, vfloat g3, vfloat g4)
{
    return vself(vmaskf_le(x, g3), x * s, (F2V(1.f) + g4) * xexpf(xlogf(x) / p) - g4);//continuous
}
#endif

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

    cmsHPROFILE oprof = nullptr;

    if (settings->HistogramWorking && consider_histogram_settings) {
        profile = icm.workingProfile;
    } else {
        profile = icm.outputProfile;

        if (icm.outputProfile.empty() || icm.outputProfile == ColorManagementParams::NoICMString) {
            profile = "sRGB";
        }

        oprof = ICCStore::getInstance()->getProfile(profile);
    }

    if (oprof) {
        const cmsUInt32Number flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE | (icm.outputBPC ? cmsFLAGS_BLACKPOINTCOMPENSATION : 0); // NOCACHE is important for thread safety

        lcmsMutex->lock();
        cmsHPROFILE LabIProf  = cmsCreateLab4Profile(nullptr);
        cmsHTRANSFORM hTransform = cmsCreateTransform (LabIProf, TYPE_Lab_DBL, oprof, TYPE_RGB_FLT, icm.outputIntent, flags);
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

void ImProcFunctions::preserv(LabImage *nprevl, LabImage *provis, int cw, int ch)
{//avoid too strong in middle values chroma when changing primaries
  float pres = 0.01f * params->icm.preser;
  float neutral = 2000000000.f;//if a2 + b2 < 200000000 scale 0..100 a and b about : 140 > a & b > -140  decrease effect 
  float medneutral = 10000000.f;//plein effect 10 > a & b > -10
  float aaneu = 1.f / (medneutral - neutral);
  float bbneu = - aaneu * neutral;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 16) nowait
#endif
    for (int i = 0; i < ch; ++i) 
        for (int j = 0; j < cw; ++j) {
            float neu = SQR(provis->a[i][j]) + SQR(provis->b[i][j]);
            if (neu < medneutral) {//plein effect
                nprevl->a[i][j] = intp(pres, provis->a[i][j], nprevl->a[i][j]); 
                nprevl->b[i][j] = intp(pres, provis->b[i][j], nprevl->b[i][j]); 
            } else if (neu < neutral) {//decrease effect
                float presred = aaneu * neu + bbneu;
                nprevl->a[i][j] = intp(pres * presred, provis->a[i][j], nprevl->a[i][j]); 
                nprevl->b[i][j] = intp(pres * presred, provis->b[i][j], nprevl->b[i][j]); 
            } 
        }
}

void ImProcFunctions::workingtrc(const Imagefloat* src, Imagefloat* dst, int cw, int ch, int mul, Glib::ustring &profile, double gampos, double slpos, int &illum, int prim, cmsHTRANSFORM &transform, bool normalizeIn, bool normalizeOut, bool keepTransForm) const
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

    if (profile == "sRGB" || profile == "Adobe RGB" || profile == "ProPhoto" || profile == "WideGamut" || profile == "BruceRGB" || profile == "Beta RGB" || profile == "BestRGB" || profile == "Rec2020" || profile == "ACESp0" || profile == "ACESp1") {
        if (settings->verbose) {
            printf("Profile=%s\n", profile.c_str());
        }
    } else {
        if (settings->verbose) {
            printf("profile not accepted\n");
        }
        return;
    }

    if (mul == -5 &&  gampos == 2.4 && slpos == 12.92310) {//must be change if we change settings RT sRGB
        //only in this case we can shortcut..all process..no gamut control..because we reduce...leads to very small differences, but big speedup
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16) if (multiThread) 
#endif

            for (int i = 0; i < ch; ++i) 
                for (int j = 0; j < cw; ++j) {
                    float r = src->r(i, j);
                    float g = src->g(i, j);
                    float b = src->b(i, j);
                    r = (Color::igammatab_srgb[r]) / 65535.f;
                    g = (Color::igammatab_srgb[g]) / 65535.f;
                    b = (Color::igammatab_srgb[b]) / 65535.f;
                    dst->r(i, j) = r;
                    dst->g(i, j) = g;
                    dst->b(i, j) = b;
                }
       return;

    }
 
    if (mul == 1 ||(params->icm.wprim == ColorManagementParams::Primaries::DEFAULT && params->icm.will == ColorManagementParams::Illuminant::DEFAULT)) {//shortcut and speedup when no call primaries and illuminant - no gamut control...in this case be careful
        GammaValues g_a; //gamma parameters
        double pwr = 1.0 / static_cast<double>(gampos);
        Color::calcGamma(pwr, slpos, g_a); // call to calcGamma with selected gamma and slope

#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < ch; ++y) {
            int x = 0;
#ifdef __SSE2__
            for (; x < cw - 3; x += 4) {
                STVFU(dst->r(y,x), F2V(65536.f) * gammalog(LVFU(src->r(y,x)), F2V(gampos), F2V(slpos), F2V(g_a[3]), F2V(g_a[4])));
                STVFU(dst->g(y,x), F2V(65536.f) * gammalog(LVFU(src->g(y,x)), F2V(gampos), F2V(slpos), F2V(g_a[3]), F2V(g_a[4])));
                STVFU(dst->b(y,x), F2V(65536.f) * gammalog(LVFU(src->b(y,x)), F2V(gampos), F2V(slpos), F2V(g_a[3]), F2V(g_a[4])));
           }
#endif
            for (; x < cw; ++x) {
                dst->r(y,x) = 65536.f * gammalog(src->r(y,x), gampos, slpos, g_a[3], g_a[4]);
                dst->g(y,x) = 65536.f * gammalog(src->g(y,x), gampos, slpos, g_a[3], g_a[4]);
                dst->b(y,x) = 65536.f * gammalog(src->b(y,x), gampos, slpos, g_a[3], g_a[4]);
            }
        }
        return;
    }
        

    float redxx = params->icm.redx;
    float redyy = params->icm.redy;
    float bluxx = params->icm.blux;
    float bluyy = params->icm.bluy;
    float grexx = params->icm.grex;
    float greyy = params->icm.grey;

    if (prim == 12) {//convert datas area to xy
        float redgraphx =  params->icm.labgridcieALow;
        float redgraphy =  params->icm.labgridcieBLow;
        float blugraphx =  params->icm.labgridcieAHigh;
        float blugraphy =  params->icm.labgridcieBHigh;
        float gregraphx =  params->icm.labgridcieGx;
        float gregraphy =  params->icm.labgridcieGy;
        redxx = 0.55f * (redgraphx + 1.f) - 0.1f;
        redxx = rtengine::LIM(redxx, 0.41f, 1.f);//limit values for xy (arbitrary)
        redyy = 0.55f * (redgraphy + 1.f) - 0.1f;
        redyy = rtengine::LIM(redyy, 0.f, 0.7f);
        bluxx = 0.55f * (blugraphx + 1.f) - 0.1f;
        bluxx = rtengine::LIM(bluxx, -0.1f, 0.5f);
        bluyy = 0.55f * (blugraphy + 1.f) - 0.1f;
        bluyy = rtengine::LIM(bluyy, -0.1f, 0.49f);
        grexx = 0.55f * (gregraphx + 1.f) - 0.1f;
        grexx = rtengine::LIM(grexx, -0.1f, 0.4f);
        greyy = 0.55f * (gregraphy + 1.f) - 0.1f;
        greyy = rtengine::LIM(greyy, 0.5f, 1.f);
    }
    //fixed crash when there is no space or too small..just a line...Possible if bx, by aligned with Gx,Gy Rx,Ry
    float ac = (greyy - redyy) / (grexx - redxx);
    float bc = greyy - ac * grexx;
    float yc = ac * bluxx + bc;
    if ((bluyy < yc + 0.0004f) &&  (bluyy > yc - 0.0004f)) {//under 0.0004 in some case crash because space too small
        return;
    }


    switch (ColorManagementParams::Primaries(prim)) {
        case ColorManagementParams::Primaries::DEFAULT: {
            break;
        }

        case ColorManagementParams::Primaries::SRGB: {
            profile = "sRGB";
            break;
        }

        case ColorManagementParams::Primaries::ADOBE_RGB: {
            profile = "Adobe RGB";
            break;
        }

        case ColorManagementParams::Primaries::PRO_PHOTO: {
            profile = "ProPhoto";
            break;
        }

        case ColorManagementParams::Primaries::REC2020: {
            profile = "Rec2020";
            break;
        }

        case ColorManagementParams::Primaries::ACES_P1: {
            profile = "ACESp1";
            break;
        }

        case ColorManagementParams::Primaries::WIDE_GAMUT: {
            profile = "WideGamut";
            break;
        }

        case ColorManagementParams::Primaries::ACES_P0: {
            profile = "ACESp0";
            break;
        }

        case ColorManagementParams::Primaries::BRUCE_RGB: {
            profile = "BruceRGB";
            break;
        }

        case ColorManagementParams::Primaries::BETA_RGB: {
            profile = "Beta RGB";
            break;
        }

        case ColorManagementParams::Primaries::BEST_RGB: {
            profile = "BestRGB";
            break;
        }

        case ColorManagementParams::Primaries::CUSTOM: {
            profile = "Custom";
            break;
        }

        case ColorManagementParams::Primaries::CUSTOM_GRID: {
            profile = "Custom";
            break;
        }
    }
    
        if (settings->verbose  && prim != 0) {
            printf("prim=%i Profile Destination=%s\n", prim, profile.c_str());
        }
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
        double tempv4 = 5003.;
        float p[6]; //primaries

        //primaries for 10 working profiles ==> output profiles
        if (profile == "WideGamut") {
            p[0] = 0.7350;    //Widegamut primaries
            p[1] = 0.2650;
            p[2] = 0.1150;
            p[3] = 0.8260;
            p[4] = 0.1570;
            p[5] = 0.0180;
            illum = toUnderlying(ColorManagementParams::Illuminant::D50);
        } else if (profile == "Adobe RGB") {
            p[0] = 0.6400;    //Adobe primaries
            p[1] = 0.3300;
            p[2] = 0.2100;
            p[3] = 0.7100;
            p[4] = 0.1500;
            p[5] = 0.0600;
            tempv4 = 6504.;
            illum = toUnderlying(ColorManagementParams::Illuminant::D65);
        } else if (profile == "sRGB") {
            p[0] = 0.6400;    // sRGB primaries
            p[1] = 0.3300;
            p[2] = 0.3000;
            p[3] = 0.6000;
            p[4] = 0.1500;
            p[5] = 0.0600;
            tempv4 = 6504.;
            illum = toUnderlying(ColorManagementParams::Illuminant::D65);
        } else if (profile == "BruceRGB") {
            p[0] = 0.6400;    // Bruce primaries
            p[1] = 0.3300;
            p[2] = 0.2800;
            p[3] = 0.6500;
            p[4] = 0.1500;
            p[5] = 0.0600;
            tempv4 = 6504.;
            illum = toUnderlying(ColorManagementParams::Illuminant::D65);
       } else if (profile == "Beta RGB") {
            p[0] = 0.6888;    // Beta primaries
            p[1] = 0.3112;
            p[2] = 0.1986;
            p[3] = 0.7551;
            p[4] = 0.1265;
            p[5] = 0.0352;
            illum = toUnderlying(ColorManagementParams::Illuminant::D50);
        } else if (profile == "BestRGB") {
            p[0] = 0.7347;    // Best primaries
            p[1] = 0.2653;
            p[2] = 0.2150;
            p[3] = 0.7750;
            p[4] = 0.1300;
            p[5] = 0.0350;
            illum = toUnderlying(ColorManagementParams::Illuminant::D50);
        } else if (profile == "Rec2020") {
            p[0] = 0.7080;    // Rec2020 primaries
            p[1] = 0.2920;
            p[2] = 0.1700;
            p[3] = 0.7970;
            p[4] = 0.1310;
            p[5] = 0.0460;
            tempv4 = 6504.;
            illum = toUnderlying(ColorManagementParams::Illuminant::D65);
        } else if (profile == "ACESp0") {
            p[0] = 0.7347;    // ACES P0 primaries
            p[1] = 0.2653;
            p[2] = 0.0000;
            p[3] = 1.0;
            p[4] = 0.0001;
            p[5] = -0.0770;
            tempv4 = 6004.;
            illum = toUnderlying(ColorManagementParams::Illuminant::D60);
        } else if (profile == "ACESp1") {
            p[0] = 0.713;    // ACES P1 primaries
            p[1] = 0.293;
            p[2] = 0.165;
            p[3] = 0.830;
            p[4] = 0.128;
            p[5] = 0.044;
            tempv4 = 6004.;
            illum = toUnderlying(ColorManagementParams::Illuminant::D60);
        } else if (profile == "ProPhoto") {
            p[0] = 0.7347;    //ProPhoto and default primaries
            p[1] = 0.2653;
            p[2] = 0.1596;
            p[3] = 0.8404;
            p[4] = 0.0366;
            p[5] = 0.0001;
            illum = toUnderlying(ColorManagementParams::Illuminant::D50);
        } else if (profile == "Custom") {
            p[0] = redxx;   
            p[1] = redyy;
            p[2] = grexx;
            p[3] = greyy;
            p[4] = bluxx;
            p[5] = bluyy;
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
        Glib::ustring ills = "D50";
        switch (ColorManagementParams::Illuminant(illum)) {
            case ColorManagementParams::Illuminant::DEFAULT:
            case ColorManagementParams::Illuminant::STDA:
            case ColorManagementParams::Illuminant::TUNGSTEN_2000K:
            case ColorManagementParams::Illuminant::TUNGSTEN_1500K: {
                break;
            }

            case ColorManagementParams::Illuminant::D41: {
                tempv4 = 4100.;
                ills = "D41";
                break;
            }

            case ColorManagementParams::Illuminant::D50: {
                tempv4 = 5003.;
                ills = "D50";
                break;
            }

            case ColorManagementParams::Illuminant::D55: {
                tempv4 = 5500.;
                ills = "D55";
                break;
            }

            case ColorManagementParams::Illuminant::D60: {
                tempv4 = 6004.;
                ills = "D60";
                break;
            }

            case ColorManagementParams::Illuminant::D65: {
                tempv4 = 6504.;
                ills = "D65";
                break;
            }

            case ColorManagementParams::Illuminant::D80: {
                tempv4 = 8000.;
                ills = "D80";
                break;
            }

            case ColorManagementParams::Illuminant::D120: {
                tempv4 = 12000.;
                ills = "D120";
                break;
            }
        }

        cmsWhitePointFromTemp(&xyD, tempv4);

        switch (ColorManagementParams::Illuminant(illum)) {
            case ColorManagementParams::Illuminant::DEFAULT:
            case ColorManagementParams::Illuminant::D55:
            case ColorManagementParams::Illuminant::D80: {
                break;
            }

            case ColorManagementParams::Illuminant::D41: {
                break;
            }

            case ColorManagementParams::Illuminant::D50: {
                xyD = {0.3457, 0.3585, 1.0}; // near LCMS values but not perfect... it's a compromise!!
                break;
            }

            case ColorManagementParams::Illuminant::D60: {
                xyD = {0.32168, 0.33767, 1.0};
                break;
            }

            case ColorManagementParams::Illuminant::D65: {
                xyD = {0.312700492, 0.329000939, 1.0};
                break;
            }

            case ColorManagementParams::Illuminant::D120: {
                xyD = {0.269669, 0.28078, 1.0};
                break;
            }

            case ColorManagementParams::Illuminant::STDA: {
                xyD = {0.447573, 0.407440, 1.0};
                ills = "stdA 2875K";
                break;
            }

            case ColorManagementParams::Illuminant::TUNGSTEN_2000K: {
                xyD = {0.526591, 0.41331, 1.0};
                ills = "Tungsten 2000K";
                break;
            }

            case ColorManagementParams::Illuminant::TUNGSTEN_1500K: {
                xyD = {0.585703, 0.393157, 1.0};
                ills = "Tungsten 1500K";
                break;
            }
        }

        //D41  0.377984  0.381229
        //D55  0.332424  0.347426
        //D80  0.293755  0.309185
        //D75  0.299021  0.314852
        cmsToneCurve* GammaTRC[3];
        GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(NULL, five, gammaParams);//5 = more smoother than 4
        cmsHPROFILE oprofdef = nullptr;

        const cmsCIExyYTRIPLE Primaries = {
            {p[0], p[1], 1.0}, // red
            {p[2], p[3], 1.0}, // green
            {p[4], p[5], 1.0}  // blue
        };
        oprofdef = cmsCreateRGBProfile(&xyD, &Primaries, GammaTRC);
        cmsWriteTag(oprofdef, cmsSigRedTRCTag, GammaTRC[0]);
        cmsWriteTag(oprofdef, cmsSigGreenTRCTag, GammaTRC[1]);
        cmsWriteTag(oprofdef, cmsSigBlueTRCTag, GammaTRC[2]);

      //to read XYZ values and illuminant
        if (rtengine::settings->verbose) {
            cmsCIEXYZ *redT = static_cast<cmsCIEXYZ*>(cmsReadTag(oprofdef, cmsSigRedMatrixColumnTag));
            cmsCIEXYZ *greenT  = static_cast<cmsCIEXYZ*>(cmsReadTag(oprofdef, cmsSigGreenMatrixColumnTag));
            cmsCIEXYZ *blueT  = static_cast<cmsCIEXYZ*>(cmsReadTag(oprofdef, cmsSigBlueMatrixColumnTag));
            printf("Illuminant=%s\n", ills.c_str());
            printf("rX=%f gX=%f bX=%f\n", redT->X, greenT->X, blueT->X);
            printf("rY=%f gY=%f bY=%f\n", redT->Y, greenT->Y, blueT->Y);
            printf("rZ=%f gZ=%f bZ=%f\n", redT->Z, greenT->Z, blueT->Z);
        }

        cmsFreeToneCurve(GammaTRC[0]);
        if (oprofdef) {
            constexpr cmsUInt32Number flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE | cmsFLAGS_BLACKPOINTCOMPENSATION | cmsFLAGS_GAMUTCHECK;
            const cmsHPROFILE iprof = ICCStore::getInstance()->getXYZProfile();
            lcmsMutex->lock();
            hTransform = cmsCreateTransform(iprof, TYPE_RGB_FLT, oprofdef, TYPE_RGB_FLT, params->icm.aRendIntent, flags);
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
