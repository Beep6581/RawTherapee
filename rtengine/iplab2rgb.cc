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

extern const Settings* settings;

namespace {

inline void copyAndClampLine(const float *src, unsigned char *dst, const int W)
{
    for (int j = 0, iy = 0; j < W; ++j) {
        float r = src[iy] * MAXVALF;
        float g = src[iy+1] * MAXVALF;
        float b = src[iy+2] * MAXVALF;
        if (r > MAXVALF || g > MAXVALF || b > MAXVALF) {
            filmlike_clip(&r, &g, &b);
        }
        dst[iy] = uint16ToUint8Rounded(CLIP(r));
        dst[iy+1] = uint16ToUint8Rounded(CLIP(g));
        dst[iy+2] = uint16ToUint8Rounded(CLIP(b));
        iy += 3;
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

            if (R > MAXVALF || G > MAXVALF || B > MAXVALF) {
                filmlike_clip(&R, &G, &B);
            }

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

        int W = lab->W;
        int H = lab->H;
        unsigned char * data = image->data;

        // cmsDoTransform is relatively expensive
#ifdef _OPENMP
        #pragma omp parallel firstprivate(lab, data, W, H)
#endif
        {
            AlignedBuffer<float> pBuf(3 * lab->W);
            AlignedBuffer<float> mBuf(3 * lab->W);

            AlignedBuffer<float> gwBuf1;
            AlignedBuffer<float> gwBuf2;

            if (gamutWarning) {
                gwBuf1.resize(3 * lab->W);
                gwBuf2.resize(3 * lab->W);
            }

            float *buffer = pBuf.data;
            float *outbuffer = mBuf.data;

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

                cmsDoTransform (monitorTransform, buffer, outbuffer, W);
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
        profile = icm.working;
        standard_gamma = true;
    } else {
        profile = icm.output;

        if (icm.output.empty() || icm.output == ColorManagementParams::NoICMString) {
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
Imagefloat* ImProcFunctions::lab2rgbOut(LabImage* lab, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm, GammaValues *ga)
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

    cmsHPROFILE oprof = nullptr;

    if (ga) {
        lcmsMutex->lock();
        ICCStore::getInstance()->getGammaArray(icm, *ga);
        oprof = ICCStore::getInstance()->createGammaProfile(icm, *ga);
        lcmsMutex->unlock();
    } else {
        oprof = ICCStore::getInstance()->getProfile(icm.output);
		Glib::ustring outtest = icm.output;
        std::string fileis_RTv2 = outtest.substr(0, 4);
		//printf("IsRTv2=%s\n", fileis_RTv2.c_str());
		if(fileis_RTv2 == "RTv2") {//Only fot ICC v2 : read tag from desc to retrieve gamma and slope save before in generate ICC v2
		//due to bug in LCMS in CmsToneCurve
		//printf("icmout=%s \n",icm.output.c_str());
			GammaValues g_b; //gamma parameters
			GammaValues gb; //gamma parameters
			const double eps = 0.000000001; // not divide by zero
			double gammatag = 2.4;
			double slopetag = 12.92;
			cmsMLU *modelDescMLU = (cmsMLU*) (cmsReadTag(oprof, cmsSigDeviceModelDescTag));
			if (modelDescMLU) {
				cmsUInt32Number count = cmsMLUgetWide(modelDescMLU, "eng", "USA", nullptr, 0);  // get buffer length first
				if (count) {
					wchar_t *buffer = new wchar_t[count];
					count = cmsMLUgetWide(modelDescMLU, "eng", "USA", buffer, count); // now put the string in the buffer
					char* cModelDesc = g_utf16_to_utf8((unsigned short int*)buffer, -1, nullptr, nullptr, nullptr); // convert to utf-8 in a buffer allocated by glib
					delete [] buffer;
					if (cModelDesc) {
						Glib::ustring modelDesc(cModelDesc);
						g_free(cModelDesc);
					//    printf("dmdd=%s\n", modelDesc.c_str());
					
						std::size_t pos = modelDesc.find("g");
						std::size_t posmid = modelDesc.find("s");
						std::size_t posend = modelDesc.find("!"); 
						std::string strgamma = modelDesc.substr(pos + 1, (posmid - pos));
						gammatag = std::stod(strgamma.c_str());
						std::string strslope = modelDesc.substr(posmid + 1, (posend - posmid));
						slopetag = std::stod(strslope.c_str());
					//	printf("gam=%f slo=%f\n", gammatag, slopetag);
					}
				} else {
					printf("Error: lab2rgbOut  /  String length is null!\n");
				}
			} else {
				printf("Error: lab2rgbOut  /  cmsReadTag/cmsSigDeviceModelDescTag failed!\n");
			}

			double pwr = 1.0 / gammatag;
			double ts = slopetag;
			double slope = slopetag == 0 ? eps : slopetag;

			int mode = 0;
			Color::calcGamma(pwr, ts, mode, g_b); // call to calcGamma with selected gamma and slope : return parameters for LCMS2
			gb[4] = g_b[3] * ts;
			gb[0] = gammatag;
			gb[1] = 1. / (1.0 + g_b[4]);
			gb[2] = g_b[4] / (1.0 + g_b[4]);
			gb[3] = 1. / slope;
			gb[5] = 0.0;
			gb[6] = 0.0;
     
			cmsToneCurve* GammaTRC[3];	
			cmsFloat64Number Parameters[7] = { gb[0],  gb[1], gb[2], gb[3], gb[4], gb[5], gb[6] } ;
		
			GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(nullptr, 5, Parameters); //5 = smoother than 4
			cmsWriteTag(oprof, cmsSigRedTRCTag, GammaTRC[0]);
			cmsWriteTag(oprof, cmsSigGreenTRCTag, GammaTRC[1]);
			cmsWriteTag(oprof, cmsSigBlueTRCTag, GammaTRC[2]);
			cmsFreeToneCurve(GammaTRC[0]);
		}  
		
		
    }

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


Imagefloat* ImProcFunctions::workingtrc(Imagefloat* working, int cw, int ch, int mul, Glib::ustring profile, double gampos, double slpos, double &ga0, double &ga1, double &ga2, double &ga3, double &ga4, double &ga5, double &ga6)
{
    TMatrix wprof;

        wprof = ICCStore::getInstance()->workingSpaceMatrix(params->icm.working);

    double dx = Color::D50x;
    double dz = Color::D50z;
    {
        dx = dz = 1.0;
    }
    double toxyz[3][3] = {
        {
            (wprof[0][0] / dx), //I have suppressed / Color::D50x
            (wprof[0][1] / dx),
            (wprof[0][2] / dx)
        }, {
            (wprof[1][0]),
            (wprof[1][1]),
            (wprof[1][2])
        }, {
            (wprof[2][0] / dz), //I have suppressed / Color::D50z
            (wprof[2][1] / dz),
            (wprof[2][2] / dz)
        }
    };

    Imagefloat* image = new  Imagefloat(cw, ch);

    double pwr;
    double ts;
    ts = slpos;

    int five = mul;

    ga6 = 0.0;
    pwr = 1.0 / gampos;

    if (gampos < 1.0) {
        pwr = gampos;
        gampos = 1. / gampos;
        five = -mul;
    }

    //  int select_temp = 1; //5003K
    const double eps = 0.000000001; // not divide by zero

    enum class ColorTemp {
        D50 = 5003,  // for Widegamut, Prophoto Best, Beta -> D50
        D65 = 6504,   // for sRGB, AdobeRGB, Bruce Rec2020  -> D65
        D60 = 6005        //for ACESP0 and AcesP1

    };
    ColorTemp temp = ColorTemp::D50;

    cmsHPROFILE oprofdef;
    float p[6]; //primaries

    if (true) {
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
        int mode = 0;
        Color::calcGamma(pwr, ts, mode, g_a); // call to calcGamma with selected gamma and slope : return parameters for LCMS2

        ga4 = g_a[3] * ts;
        ga0 = gampos;
        ga1 = 1. / (1.0 + g_a[4]);
        ga2 = g_a[4] / (1.0 + g_a[4]);
        ga3 = 1. / slpos;
        ga5 = 0.0;
		ga6 = 0.0;
       // printf("ga0=%f ga1=%f ga2=%f ga3=%f ga4=%f\n", ga0, ga1, ga2, ga3, ga4);

        cmsCIExyY       xyD;

        cmsCIExyYTRIPLE Primaries = {
            {p[0], p[1], 1.0}, // red
            {p[2], p[3], 1.0}, // green
            {p[4], p[5], 1.0}  // blue
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
        cmsWhitePointFromTemp(&xyD, (double)temp);
        GammaTRC[0] = GammaTRC[1] = GammaTRC[2] =   cmsBuildParametricToneCurve(NULL, five, Parameters);//5 = more smoother than 4
        oprofdef = cmsCreateRGBProfile(&xyD, &Primaries, GammaTRC);
        cmsFreeToneCurve(GammaTRC[0]);
    }

    if (oprofdef) {
        #pragma omp parallel for if (multiThread)

        for (int i = 0; i < ch; i++) {
            float* rr = working->r(i);
            float* rg = working->g(i);
            float* rb = working->b(i);
			
            float* xa = (float*)image->r(i);
            float* ya = (float*)image->g(i);
            float* za = (float*)image->b(i);
			
			
            for (int j = 0; j < cw; j++) {
                float r1 = rr[j];
                float g1 = rg[j];
                float b1 = rb[j];


                float x_ = toxyz[0][0] * r1 + toxyz[0][1] * g1 + toxyz[0][2] * b1;
                float y_ = toxyz[1][0] * r1 + toxyz[1][1] * g1 + toxyz[1][2] * b1;
                float z_ = toxyz[2][0] * r1 + toxyz[2][1] * g1 + toxyz[2][2] * b1;
	
			
                xa[j] = ( x_) ;
                ya[j] = ( y_);
                za[j] = ( z_);

            }
        }

        cmsUInt32Number flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;


        lcmsMutex->lock();
        cmsHPROFILE iprof = ICCStore::getInstance()->getXYZProfile();
		//   cmsHTRANSFORM hTransform = cmsCreateTransform(iprof, TYPE_RGB_16, oprofdef, TYPE_RGB_16, params->icm.outputIntent,  cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
		cmsHTRANSFORM hTransform = cmsCreateTransform(iprof, TYPE_RGB_FLT, oprofdef, TYPE_RGB_FLT, params->icm.outputIntent, flags);
		lcmsMutex->unlock();

        image->ExecCMSTransform2(hTransform);
		
        cmsDeleteTransform(hTransform);
        image->normalizeFloatTo65535();

    }


    return image;

}


}
