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

extern void filmlike_clip(float *r, float *g, float *b);

namespace {

inline void clipLAB(float iL, float ia, float ib, float &oL, float &oa, float &ob, const float scale, const float wp[3][3], const float wip[3][3])
{
    if (iL < 0.f) {
        oL = oa = ob = 0.f;
    } else if (iL > 32768.f) {
        
        float X, Y, Z;
        float r, g, b;
        Color::Lab2XYZ(iL, ia, ib, X, Y, Z);
        Color::xyz2rgb(X, Y, Z, r, g, b, wip);
        filmlike_clip(&r, &g, &b);
        Color::rgbxyz(r, g, b, X, Y, Z, wp);
        Color::XYZ2Lab(X, Y, Z, oL, oa, ob);
        oL /= scale;
        oa /= scale;
        ob /= scale;
        
        // oL = 32768.f / scale;
        // oa = ob = 0.f;
    } else {
        oL = iL / scale;
        oa = ia / scale;
        ob = ib / scale;
    }
}


inline void clipLAB(float iL, float ia, float ib, double &oL, double &oa, double &ob, const float scale, const float wp[3][3], const float wip[3][3])
{
    float tL, ta, tb;
    clipLAB(iL, ia, ib, tL, ta, tb, scale, wp, wip);
    oL = tL;
    oa = ta;
    ob = tb;
}

} // namespace

extern const Settings* settings;

#define DECLARE_WORKING_MATRICES_(space) \
    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix ( space ); \
    const float wp[3][3] = {                                            \
        {static_cast<float> (wprof[0][0]), static_cast<float> (wprof[0][1]), static_cast<float> (wprof[0][2])}, \
        {static_cast<float> (wprof[1][0]), static_cast<float> (wprof[1][1]), static_cast<float> (wprof[1][2])}, \
        {static_cast<float> (wprof[2][0]), static_cast<float> (wprof[2][1]), static_cast<float> (wprof[2][2])} \
    };                                                                  \
                                                                        \
    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix ( space ); \
    const float wip[3][3] = {                                           \
        {static_cast<float> (wiprof[0][0]), static_cast<float> (wiprof[0][1]), static_cast<float> (wiprof[0][2])}, \
        {static_cast<float> (wiprof[1][0]), static_cast<float> (wiprof[1][1]), static_cast<float> (wiprof[1][2])}, \
        {static_cast<float> (wiprof[2][0]), static_cast<float> (wiprof[2][1]), static_cast<float> (wiprof[2][2])} \
    }
    

// Used in ImProcCoordinator::updatePreviewImage  (rtengine/improccoordinator.cc)
//         Crop::update                           (rtengine/dcrop.cc)
//         Thumbnail::processImage                (rtengine/rtthumbnail.cc)
//
// If monitorTransform, divide by 327.68 then apply monitorTransform (which can integrate soft-proofing)
// otherwise divide by 327.68, convert to xyz and apply the sRGB transform, before converting with gamma2curve
void ImProcFunctions::lab2monitorRgb (LabImage* lab, Image8* image)
{
    DECLARE_WORKING_MATRICES_(params->icm.working);
    
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

                for (int j = 0; j < W; j++) {
                    clipLAB(rL[j], ra[j], rb[j], buffer[iy], buffer[iy+1], buffer[iy+2], 327.68f, wp, wip);
                    iy += 3;
                }

                cmsDoTransform (monitorTransform, buffer, data + ix, W);
            }
        } // End of parallelization
    } else {

        int W = lab->W;
        int H = lab->H;
        unsigned char * data = image->data;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

        for (int i = 0; i < H; ++i) {
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
            int ix = i * 3 * W;

            float R, G, B;
            float x_, y_, z_;
            float L, a, b;

            for (int j = 0; j < W; ++j) {

                //float L1=rL[j],a1=ra[j],b1=rb[j];//for testing
                clipLAB(rL[j], ra[j], rb[j], L, a, b, 1.f, wp, wip);

                Color::Lab2XYZ(L, a, b, x_, y_, z_ );

                Color::xyz2srgb(x_, y_, z_, R, G, B);

                /* copy RGB */
                //int R1=((int)gamma2curve[(R)])
                data[ix++] = uint16ToUint8Rounded(Color::gamma2curve[R]);
                data[ix++] = uint16ToUint8Rounded(Color::gamma2curve[G]);
                data[ix++] = uint16ToUint8Rounded(Color::gamma2curve[B]);
            }
        }
    }
}



// Used in ImProcCoordinator::updatePreviewImage  (rtengine/improccoordinator.cc)
//         Crop::update                           (rtengine/dcrop.cc)
//
// Generate an Image8
//
// If output profile used, divide by 327.68 then apply the "profile" profile (eventually with a standard gamma)
// otherwise divide by 327.68, convert to xyz and apply the RGB transform, before converting with gamma2curve
Image8* ImProcFunctions::lab2rgb (LabImage* lab, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm, bool consider_histogram_settings)
{
    DECLARE_WORKING_MATRICES_(icm.working);
    
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
    Glib::ustring profile;

    bool standard_gamma;

    if(settings->HistogramWorking && consider_histogram_settings) {
        profile = icm.working;
        standard_gamma = true;
    } else {
        profile = icm.output;
        if (icm.output.empty() || icm.output == ColorManagementParams::NoICMString) {
            profile = "sRGB";
        }
        standard_gamma = false;
    }

    cmsHPROFILE oprof = ICCStore::getInstance()->getProfile (profile);

    if (oprof) {
        cmsHPROFILE oprofG = oprof;

        if (standard_gamma) {
            oprofG = ICCStore::makeStdGammaProfile(oprof);
        }

        cmsUInt32Number flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;
        if (icm.outputBPC) {
            flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
        }
        lcmsMutex->lock ();
        cmsHPROFILE LabIProf  = cmsCreateLab4Profile(nullptr);
        cmsHTRANSFORM hTransform = cmsCreateTransform (LabIProf, TYPE_Lab_DBL, oprofG, TYPE_RGB_8, icm.outputIntent, flags);  // NOCACHE is important for thread safety
        cmsCloseProfile(LabIProf);
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
                    clipLAB(rL[j], ra[j], rb[j], buffer[iy], buffer[iy+1], buffer[iy+2], 327.68f, wp, wip);
                    iy += 3;
                }

                cmsDoTransform (hTransform, buffer, data + ix, cw);
            }
        } // End of parallelization

        cmsDeleteTransform(hTransform);

        if (oprofG != oprof) {
            cmsCloseProfile(oprofG);
        }
    } else {

        const auto xyz_rgb = ICCStore::getInstance()->workingSpaceInverseMatrix (profile);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

        for (int i = cy; i < cy + ch; ++i) {
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
            int ix = 3 * i * cw;

            float R, G, B;
            float x_, y_, z_;
            float L, a, b;

            for (int j = cx; j < cx + cw; ++j) {
                clipLAB(rL[j], ra[j], rb[j], L, a, b, 1.f, wp, wip);
                Color::Lab2XYZ(rL[j], ra[j], rb[j], x_, y_, z_);

                Color::xyz2rgb(x_, y_, z_, R, G, B, xyz_rgb);

                image->data[ix++] = uint16ToUint8Rounded(Color::gamma2curve[R]);
                image->data[ix++] = uint16ToUint8Rounded(Color::gamma2curve[G]);
                image->data[ix++] = uint16ToUint8Rounded(Color::gamma2curve[B]);
            }
        }
    }

    return image;
}


/** @brief Convert the final Lab image to the output RGB color space
 *
 * Used in processImage   (rtengine/simpleprocess.cc)
 *
 * Provide a pointer to a 7 floats array for "ga" (uninitialized ; this array will be filled with the gamma values) if you want
 * to use the custom gamma scenario. Thoses gamma values will correspond to the ones of the chosen standard output profile
 * (Prophoto if non standard output profile given)
 *
 * If "ga" is NULL, then we're considering standard gamma with the chosen output profile.
 *
 * Generate an Image16
 *
 * If a custom gamma profile can be created, divide by 327.68, convert to xyz and apply the custom gamma transform
 * otherwise divide by 327.68, convert to xyz and apply the sRGB transform, before converting with gamma2curve
 */
Imagefloat* ImProcFunctions::lab2rgbOut (LabImage* lab, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm, GammaValues *ga)
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

    Imagefloat* image = new Imagefloat (cw, ch);

    cmsHPROFILE oprof = nullptr;
    if (ga) {
        lcmsMutex->lock ();
        ICCStore::getInstance()->getGammaArray(icm, *ga);
        oprof = ICCStore::getInstance()->createGammaProfile(icm, *ga);
        lcmsMutex->unlock ();
    } else {
        oprof = ICCStore::getInstance()->getProfile (icm.output);
    }

    if (oprof) {
        cmsUInt32Number flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;
        if (icm.outputBPC) {
            flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
        }
        lcmsMutex->lock ();
        cmsHPROFILE iprof = cmsCreateLab4Profile(nullptr);
        cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_Lab_FLT, oprof, TYPE_RGB_FLT, icm.outputIntent, flags);
        lcmsMutex->unlock ();

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

                setUnlessOOG(image->r(i - cy, j - cx), Color::gamma2curve[CLIP(R)]);
                setUnlessOOG(image->g(i - cy, j - cx), Color::gamma2curve[CLIP(G)]);
                setUnlessOOG(image->b(i - cy, j - cx), Color::gamma2curve[CLIP(B)]);
            }
        }
    }

    return image;
}

}
