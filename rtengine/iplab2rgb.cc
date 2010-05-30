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

extern const Settings* settings;

void ImProcFunctions::lab2rgb (LabImage* lab, Image8* image) {

	if (monitorTransform) {
	    int ix = 0;
        short* buffer = new short [3*lab->W];
		for (int i=0; i<lab->H; i++) {
			unsigned short* rL = lab->L[i];
			short* ra = lab->a[i];
			short* rb = lab->b[i];
			int iy = 0;
			for (int j=0; j<lab->W; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);

                y_ = ycache[y_];
				x_ = xcache[x_];
				z_ = zcache[z_];

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
			unsigned short* rL = lab->L[i];
			short* ra = lab->a[i];
			short* rb = lab->b[i];
			int ix = 3*i*lab->W;
			for (int j=0; j<lab->W; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);

                y_ = ycache[y_];
				x_ = xcache[x_];
				z_ = zcache[z_];

				/* XYZ-D50 to RGB */
				int R = (25689*x_-13261*y_-4022*z_) >> 13;
				int G = (-8017*x_+15697*y_+274*z_) >> 13;
				int B = (590*x_-1877*y_+11517*z_) >> 13;

				/* copy RGB */
				image->data[ix++] = gamma2curve[CLIP(R)] >> 8;
				image->data[ix++] = gamma2curve[CLIP(G)] >> 8;
				image->data[ix++] = gamma2curve[CLIP(B)] >> 8;
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

    cmsHPROFILE oprof = iccStore.getProfile (profile);

    if (oprof) {
        cmsHPROFILE iprof = iccStore.getXYZProfile ();
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16, oprof, TYPE_RGB_8, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
        int ix = 0;
        short* buffer = new short [3*cw];
        for (int i=cy; i<cy+ch; i++) {
            unsigned short* rL = lab->L[i];
            short* ra = lab->a[i];
            short* rb = lab->b[i];
            int iy = 0;
            for (int j=cx; j<cx+cw; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);

                y_ = ycache[y_];
				x_ = xcache[x_];
				z_ = zcache[z_];

                buffer[iy++] = CLIP(x_);
                buffer[iy++] = CLIP(y_);
                buffer[iy++] = CLIP(z_);
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
            unsigned short* rL = lab->L[i];
            short* ra = lab->a[i];
            short* rb = lab->b[i];
			int ix = 3*i*cw;
            for (int j=cx; j<cx+cw; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);

                y_ = ycache[y_];
				x_ = xcache[x_];
				z_ = zcache[z_];

                int R = (25689*x_-13261*y_-4022*z_) >> 13;
                int G = (-8017*x_+15697*y_+274*z_) >> 13;
                int B = (590*x_-1877*y_+11517*z_) >> 13;

                image->data[ix++] = gamma2curve[CLIP(R)] >> 8;
                image->data[ix++] = gamma2curve[CLIP(G)] >> 8;
                image->data[ix++] = gamma2curve[CLIP(B)] >> 8;
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

    cmsHPROFILE oprof = iccStore.getProfile (profile);

    if (oprof) {
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
			unsigned short* rL = lab->L[i];
			short* ra = lab->a[i];
			short* rb = lab->b[i];
			short* xa = (short*)image->r[i-cy];
			short* ya = (short*)image->g[i-cy];
			short* za = (short*)image->b[i-cy];
			for (int j=cx; j<cx+cw; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);

                y_ = ycache[y_];
				x_ = xcache[x_];
				z_ = zcache[z_];

				xa[j-cx] = CLIP(x_);
				ya[j-cx] = CLIP(y_);
				za[j-cx] = CLIP(z_);
			}
		}
        cmsHPROFILE iprof = iccStore.getXYZProfile ();
        lcmsMutex->lock ();
		cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16_PLANAR, oprof, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
		cmsDoTransform (hTransform, image->data, image->data, image->planestride/2);
		cmsDeleteTransform(hTransform);
	}
	else {
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
			unsigned short* rL = lab->L[i];
			short* ra = lab->a[i];
			short* rb = lab->b[i];
			for (int j=cx; j<cx+cw; j++) {

                int y_ = rL[j];
                int x_ = rL[j]+10486+ra[j]*152/chroma_scale+141556;
                int z_ = rL[j]+10486-rb[j]*380/chroma_scale+369619;

                x_ = CLIPTO(x_,0,369820);
                y_ = CLIPTO(y_,0,825745);

                y_ = ycache[y_];
				x_ = xcache[x_];
				z_ = zcache[z_];

				int R = (25689*x_-13261*y_-4022*z_) >> 13;
				int G = (-8017*x_+15697*y_+274*z_) >> 13;
				int B = (590*x_-1877*y_+11517*z_) >> 13;

				image->r[i-cy][j-cx] = gamma2curve[CLIP(R)];
				image->g[i-cy][j-cx] = gamma2curve[CLIP(G)];
				image->b[i-cy][j-cx] = gamma2curve[CLIP(B)];
			}
		}
	}
    return image;
}

}
