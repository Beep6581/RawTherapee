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
#include <rtengine.h>
#include <improcfun.h>
#include <glibmm.h>
#include <iccstore.h>
#include <iccmatrices.h>
#include <mytime.h>
//#include <sRGBgamutbdy.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {

#undef CLIP
#undef CLIPTO
#undef CMAXVAL

#define CMAXVAL 0xffff
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIP01(a) ((a)>0?((a)<1?(a):1):0)
	
	
#define D50x 0.96422
#define D50z 0.82521

extern const Settings* settings;
	
const double (*wprof[])[3]  = {xyz_sRGB, xyz_adobe, xyz_prophoto, xyz_widegamut, xyz_bruce, xyz_beta, xyz_best};
const double (*iwprof[])[3] = {sRGB_xyz, adobe_xyz, prophoto_xyz, widegamut_xyz, bruce_xyz, beta_xyz, best_xyz};
const char* wprofnames[] = {"sRGB", "Adobe RGB", "ProPhoto", "WideGamut", "BruceRGB", "Beta RGB", "BestRGB"};
const int numprof = 7;

void ImProcFunctions::lab2rgb (LabImage* lab, Image8* image) {
	//MyTime tBeg,tEnd;
 //   tBeg.set();
	//gamutmap(lab);

	if (monitorTransform) {
        
        // cmsDoTransform is relatively expensive
        #pragma omp parallel for if (multiThread)
		for (int i=0; i<lab->H; i++) {
            float buffer[3*lab->W];
            float g;

            const int ix = i * 3 * lab->W;
            int iy = 0;

			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];

            float fy,fx,fz,x_,y_,z_;

			for (int j=0; j<lab->W; j++) {
								
				fy = (0.00862069 * rL[j]) / 327.68 + 0.137932; // (L+16)/116
				fx = (0.002 * ra[j]) / 327.68 + fy;
				fz = fy - (0.005 * rb[j]) / 327.68;
				
				x_ = f2xyz(fx)*D50x;//should this be 32767???  buffer is short int !!!
				y_ = f2xyz(fy);
				z_ = f2xyz(fz)*D50z;

                buffer[iy++] = CLIP01(x_);
                buffer[iy++] = CLIP01(y_);
                buffer[iy++] = CLIP01(z_);
			}

            if (settings->LCMSSafeMode) lcmsMutex->lock ();
            cmsDoTransform (monitorTransform, buffer, image->data + ix, lab->W);
            if (settings->LCMSSafeMode) lcmsMutex->unlock ();
		}
        
	} else {

		#pragma omp parallel for if (multiThread)
		for (int i=0; i<lab->H; i++) {
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			int ix = i * 3 * lab->W;
			float g;

			float R,G,B;
            float fy,fx,fz,x_,y_,z_;

			for (int j=0; j<lab->W; j++) {
			
				//float L1=rL[j],a1=ra[j],b1=rb[j];//for testing
				
				fy = (0.00862069 * rL[j]) / 327.68 + 0.137932; // (L+16)/116
				fx = (0.002 * ra[j]) / 327.68 + fy;
				fz = fy - (0.005 * rb[j]) / 327.68;
				
				x_ = 65535.0 * f2xyz(fx)*D50x;
				y_ = 65535.0 * f2xyz(fy);
				z_ = 65535.0 * f2xyz(fz)*D50z;

				xyz2srgb(x_,y_,z_,R,G,B);

				/* copy RGB */
				image->data[ix++] = (int)gamma2curve[(R)] >> 8;
				image->data[ix++] = (int)gamma2curve[(G)] >> 8;
				image->data[ix++] = (int)gamma2curve[(B)] >> 8;
			}
		}
	}

    //tEnd.set();
    //printf("lab2rgb %i %d\n", lab->W, tEnd.etime(tBeg));
}

Image8* ImProcFunctions::lab2rgb (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile) {

	//gamutmap(lab);
	
    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw>lab->W) cw = lab->W-cx;
    if (cy+ch>lab->H) ch = lab->H-cy;

    Image8* image = new Image8 (cw, ch);

    cmsHPROFILE oprof = iccStore->getProfile (profile);

    if (oprof) {
        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16, oprof, TYPE_RGB_8, settings->colorimetricIntent,
            settings->LCMSSafeMode ? cmsFLAGS_NOOPTIMIZE : cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety
        lcmsMutex->unlock ();

        // cmsDoTransform is relatively expensive
		#pragma omp parallel for if (multiThread)
        for (int i=cy; i<cy+ch; i++) {
            short buffer [3*cw];
            float g;

            const int ix = i * 3 * cw;
            int iy = 0;

            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];

            for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

                buffer[iy++] = CLIP((int)x_);
                buffer[iy++] = CLIP((int)y_);
                buffer[iy++] = CLIP((int)z_);
            }

            if (settings->LCMSSafeMode) lcmsMutex->lock ();
            cmsDoTransform (hTransform, buffer, image->data + ix, cw);
            if (settings->LCMSSafeMode) lcmsMutex->unlock ();
        }

        cmsDeleteTransform(hTransform);
    } else {
		
		float rgb_xyz[3][3];
		
		for (int i=0; i<numprof; i++) {
			if (profile==wprofnames[i]) {
				for (int m=0; m<3; m++) 
					for (int n=0; n<3; n++) {
						rgb_xyz[m][n] = iwprof[i][m][n];
					}
				break;
			}
		}
		
		#pragma omp parallel for if (multiThread)
        for (int i=cy; i<cy+ch; i++) {
			float g;
			float R,G,B;
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
			int ix = 3*i*cw;
            for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

				xyz2rgb(x_,y_,z_,R,G,B,rgb_xyz);

                image->data[ix++] = (int)gamma2curve[(R)] >> 8;
                image->data[ix++] = (int)gamma2curve[(G)] >> 8;
                image->data[ix++] = (int)gamma2curve[(B)] >> 8;
            }
        }
    }
    return image;
}

Image16* ImProcFunctions::lab2rgb16 (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile) {
	
	//gamutmap(lab);

    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw>lab->W) cw = lab->W-cx;
    if (cy+ch>lab->H) ch = lab->H-cy;

    Image16* image = new Image16 (cw, ch);

    cmsHPROFILE oprof = iccStore->getProfile (profile);

    if (oprof) {
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
			float g;
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			short* xa = (short*)image->r[i-cy];
			short* ya = (short*)image->g[i-cy];
			short* za = (short*)image->b[i-cy];
			for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

				xa[j-cx] = CLIP((int)x_);
				ya[j-cx] = CLIP((int)y_);
				za[j-cx] = CLIP((int)z_);
			}
		}

        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
		cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16_PLANAR, oprof, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, cmsFLAGS_NOOPTIMIZE);
        lcmsMutex->unlock ();

		cmsDoTransform (hTransform, image->data, image->data, image->planestride);
		cmsDeleteTransform(hTransform);
	} else {
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
			float g;
			float R,G,B;
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

				xyz2srgb(x_,y_,z_,R,G,B);

				image->r[i-cy][j-cx] = (int)gamma2curve[(R)];
				image->g[i-cy][j-cx] = (int)gamma2curve[(G)];
				image->b[i-cy][j-cx] = (int)gamma2curve[(B)];
			}
		}
	}
    return image;
}
	
//#include "sRGBgamutbdy.cc"

}
