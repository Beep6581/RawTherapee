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
	
	
#define epsilon 0.00885645 //216/24389
#define kappa 903.2963 //24389/27
#define kappainv 0.00110706 //inverse of kappa
#define kapeps 8 // kappa*epsilon
#define Lab2xyz(f) (( (g=f*f*f) > epsilon) ? g : (116*f-16)*kappainv)
	
#define D50x 0.96422
#define D50z 0.82521

extern const Settings* settings;

void ImProcFunctions::lab2rgb (LabImage* lab, Image8* image) {

	if (monitorTransform) {
	    int ix = 0;
		float g;
        short* buffer = new short [3*lab->W];
		for (int i=0; i<lab->H; i++) {
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			int iy = 0;
			for (int j=0; j<lab->W; j++) {
								
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535*Lab2xyz(fx)*D50x;
				float y_ = 65535*Lab2xyz(fy);
				float z_ = 65535*Lab2xyz(fz)*D50z;

                buffer[iy++] = CLIP(x_);
                buffer[iy++] = CLIP(y_);
                buffer[iy++] = CLIP(z_);
			}
            cmsDoTransform (monitorTransform, buffer, image->data + ix, lab->W);
            ix += 3*lab->W;
		}
        delete [] buffer;
	}
	else {
		#pragma omp parallel for if (multiThread)
		for (int i=0; i<lab->H; i++) {
			float g;
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			int ix = 3*i*lab->W;
			for (int j=0; j<lab->W; j++) {
				
				float L1=rL[j],a1=ra[j],b1=rb[j];//for testing
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535*Lab2xyz(fx)*D50x;
				float y_ = 65535*Lab2xyz(fy);
				float z_ = 65535*Lab2xyz(fz)*D50z;

				/* XYZ-D50 to RGB */
				int R = (int)( 3.1338561*x_ - 1.6168667*y_ - 0.4906146*z_);
				int G = (int)(-0.9787684*x_ + 1.9161415*y_ + 0.0334540*z_);
				int B = (int)( 0.0719453*x_ - 0.2289914*y_ + 1.4052427*z_);
				
				// XYZ-D65 to RGB 
				//3.2404542 -1.5371385 -0.4985314
				//-0.9692660  1.8760108  0.0415560
				//0.0556434 -0.2040259  1.0572252

				/* copy RGB */
				image->data[ix++] = (int)gamma2curve[CLIP(R)] >> 8;
				image->data[ix++] = (int)gamma2curve[CLIP(G)] >> 8;
				image->data[ix++] = (int)gamma2curve[CLIP(B)] >> 8;
			}
		}
	}
}

Image8* ImProcFunctions::lab2rgb (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile) {

    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw>lab->W) cw = lab->W-cx;
    if (cy+ch>lab->H) ch = lab->H-cy;

    Image8* image = new Image8 (cw, ch);

    cmsHPROFILE oprof = iccStore->getProfile (profile);

    if (oprof) {
        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16, oprof, TYPE_RGB_8, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
        int ix = 0;
		float g;
        short* buffer = new short [3*cw];
        for (int i=cy; i<cy+ch; i++) {
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
            int iy = 0;
            for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535*Lab2xyz(fx)*D50x;
				float y_ = 65535*Lab2xyz(fy);
				float z_ = 65535*Lab2xyz(fz)*D50z;

                buffer[iy++] = CLIP((int)x_);
                buffer[iy++] = CLIP((int)y_);
                buffer[iy++] = CLIP((int)z_);
            }
            cmsDoTransform (hTransform, buffer, image->data + ix, cw);
            ix += 3*cw;
        }
        delete [] buffer;
        cmsDeleteTransform(hTransform);
    }
    else {
		#pragma omp parallel for if (multiThread)
        for (int i=cy; i<cy+ch; i++) {
			float g;
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
			int ix = 3*i*cw;
            for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535*Lab2xyz(fx)*D50x;
				float y_ = 65535*Lab2xyz(fy);
				float z_ = 65535*Lab2xyz(fz)*D50z;

				/* XYZ-D50 to RGB */
				int R = (int)(3.1338561*x_ - 1.6168667*y_ - 0.4906146*z_);
				int G = (int)(-0.9787684*x_ + 1.9161415*y_ + 0.0334540*z_);
				int B = (int)(0.0719453*x_ -0.2289914*y_ + 1.4052427*z_);

                image->data[ix++] = (int)gamma2curve[CLIP(R)] >> 8;
                image->data[ix++] = (int)gamma2curve[CLIP(G)] >> 8;
                image->data[ix++] = (int)gamma2curve[CLIP(B)] >> 8;
            }
        }
    }
    return image;
}

Image16* ImProcFunctions::lab2rgb16 (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile) {

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
				
				float x_ = 65535*Lab2xyz(fx)*D50x;
				float y_ = 65535*Lab2xyz(fy);
				float z_ = 65535*Lab2xyz(fz)*D50z;

				xa[j-cx] = CLIP((int)x_);
				ya[j-cx] = CLIP((int)y_);
				za[j-cx] = CLIP((int)z_);
			}
		}
        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
		cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16_PLANAR, oprof, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
		cmsDoTransform (hTransform, image->data, image->data, image->planestride/2);
		cmsDeleteTransform(hTransform);
	}
	else {
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
			float g;
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535*Lab2xyz(fx)*D50x;
				float y_ = 65535*Lab2xyz(fy);
				float z_ = 65535*Lab2xyz(fz)*D50z;

				/* XYZ-D50 to RGB */
				int R = (int)(3.1338561*x_ - 1.6168667*y_ - 0.4906146*z_);
				int G = (int)(-0.9787684*x_ + 1.9161415*y_ + 0.0334540*z_);
				int B = (int)(0.0719453*x_ -0.2289914*y_ + 1.4052427*z_);

				image->r[i-cy][j-cx] = gamma2curve[CLIP(R)];
				image->g[i-cy][j-cx] = gamma2curve[CLIP(G)];
				image->b[i-cy][j-cx] = gamma2curve[CLIP(B)];
			}
		}
	}
    return image;
}

}
