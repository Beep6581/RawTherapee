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
#ifdef _OPENMP
#include <omp.h>
#endif
#include <minmax.h>
#include <gauss.h>
#include <bilateral2.h>

namespace rtengine {

#undef CLIP
#undef CMAXVAL
#undef ABS

#define CMAXVAL 0xffff
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define ABS(a) ((a)<0?-(a):(a))

void ImProcFunctions::dcdamping (float** aI, unsigned short** aO, float damping, int W, int H) {

#ifdef _OPENMP
#pragma omp for
#endif
	for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
            float I = aI[i][j];
            float O = (float)aO[i][j];
            if (O==0.0 || I==0.0) {
                aI[i][j] = 0.0;
                continue;
            }
            float U = -(O * log(I/O) - I + O) * 2.0 / (damping*damping);
            U = MIN(U,1.0);
            U = U*U*U*U*(5.0-U*4.0);
            aI[i][j] = (O - I) / I * U + 1.0;
        }
}

void ImProcFunctions::deconvsharpening (LabImage* lab, unsigned short** b2) {

    if (params->sharpening.enabled==false || params->sharpening.deconvamount<1)
        return;

    int W = lab->W, H = lab->H;

    float** tmpI = new float*[H];
    for (int i=0; i<H; i++) {
        tmpI[i] = new float[W];
        for (int j=0; j<W; j++)
            tmpI[i][j] = (float)lab->L[i][j];
    }

    float** tmp = (float**)b2;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(W,H));
    float damping = params->sharpening.deconvdamping / 5.0;
    bool needdamp = params->sharpening.deconvdamping > 0;
    for (int k=0; k<params->sharpening.deconviter; k++) {

    	// apply blur function (gaussian blur)
        gaussHorizontal<float> (tmpI, tmp, buffer, W, H, params->sharpening.deconvradius / scale, multiThread);
        gaussVertical<float>   (tmp, tmp,  buffer, W, H, params->sharpening.deconvradius / scale, multiThread);

    	if (!needdamp) {
#ifdef _OPENMP
#pragma omp for
#endif
            for (int i=0; i<H; i++)
                for (int j=0; j<W; j++)
                    if (tmp[i][j]>0)
                        tmp[i][j] = (float)lab->L[i][j] / tmp[i][j];
        }
        else
			dcdamping (tmp, lab->L, damping, W, H);

        gaussHorizontal<float> (tmp, tmp, buffer, W, H, params->sharpening.deconvradius / scale, multiThread);
        gaussVertical<float>   (tmp, tmp, buffer, W, H, params->sharpening.deconvradius / scale, multiThread);

#ifdef _OPENMP
#pragma omp for
#endif
        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++)
                tmpI[i][j] = tmpI[i][j] * tmp[i][j];
		} // end for
    delete buffer;
    } // end parallel

#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++)
            lab->L[i][j] = lab->L[i][j]*(100-params->sharpening.deconvamount) / 100 + (int)CLIP(tmpI[i][j])*params->sharpening.deconvamount / 100;



    for (int i=0; i<H; i++)
        delete [] tmpI[i];
    delete [] tmpI;
}

void ImProcFunctions::sharpening (LabImage* lab, unsigned short** b2) {

    if (params->sharpening.method=="rld") {
        deconvsharpening (lab, b2);
        return;
    }

    if (params->sharpening.enabled==false || params->sharpening.amount<1 || lab->W<8 || lab->H<8)
        return;

    int W = lab->W, H = lab->H;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    unsigned short** b3;

    AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(W,H));
    if (params->sharpening.edgesonly==false) {

        gaussHorizontal<unsigned short> (lab->L, b2, buffer, W, H, params->sharpening.radius / scale, multiThread);
        gaussVertical<unsigned short>   (b2,     b2, buffer, W, H, params->sharpening.radius / scale, multiThread);
    }
    else {
        b3 = new unsigned short*[H];
        for (int i=0; i<H; i++)
            b3[i] = new unsigned short[W];

		bilateral<unsigned short, unsigned int> (lab->L, (unsigned short**)b3, b2, W, H, params->sharpening.edges_radius / scale, params->sharpening.edges_tolerance, multiThread);
		gaussHorizontal<unsigned short> (b3, b2, buffer, W, H, params->sharpening.radius / scale, multiThread);
		gaussVertical<unsigned short>   (b2, b2, buffer, W, H, params->sharpening.radius / scale, multiThread);
    }
    delete buffer;

    unsigned short** base = lab->L;
    if (params->sharpening.edgesonly)
        base = b3;

    if (params->sharpening.halocontrol==false) {
		#pragma omp for
    	for (int i=0; i<H; i++)
            for (int j=0; j<W; j++) {
                int diff = base[i][j] - b2[i][j];
                if (ABS(diff)>params->sharpening.threshold) {
                    int val = lab->L[i][j] + params->sharpening.amount * diff / 100;
                    lab->L[i][j] = CLIP(val);
                }
            }
    }
    else
		sharpenHaloCtrl (lab, b2, base, W, H);

    if (params->sharpening.edgesonly) {
        for (int i=0; i<H; i++)
            delete [] b3[i];
        delete [] b3;
    }
    }
}

void ImProcFunctions::sharpenHaloCtrl (LabImage* lab, unsigned short** blurmap, unsigned short** base, int W, int H) {

    int scale = 100 * (100-params->sharpening.halocontrol_amount);
    unsigned short** nL = base;
	#pragma omp parallel for if (multiThread)
    for (int i=2; i<H-2; i++) {
        int max1 = 0, max2 = 0, min1 = 0, min2 = 0, maxn, minn, np1, np2, np3, min, max;
        for (int j=2; j<W-2; j++) {
            int diff = base[i][j] - blurmap[i][j];
            if (ABS(diff) > params->sharpening.threshold) {
                // compute maximum/minimum in a delta environment
                np1 = 2*(nL[i-2][j] + nL[i-2][j+1] + nL[i-2][j+2] + nL[i-1][j] + nL[i-1][j+1] + nL[i-1][j+2] + nL[i][j] + nL[i][j+1] + nL[i][j+2]) / 27 + nL[i-1][j+1] / 3;
                np2 = 2*(nL[i-1][j] + nL[i-1][j+1] + nL[i-1][j+2] + nL[i][j] + nL[i][j+1] + nL[i][j+2] + nL[i+1][j] + nL[i+1][j+1] + nL[i+1][j+2]) / 27 + nL[i][j+1] / 3;
                np3 = 2*(nL[i][j] + nL[i][j+1] + nL[i][j+2] + nL[i+1][j] + nL[i+1][j+1] + nL[i+1][j+2] + nL[i+2][j] + nL[i+2][j+1] + nL[i+2][j+2]) / 27 + nL[i+1][j+1] / 3;
                MINMAX3(np1,np2,np3,maxn,minn);
                MAX3(max1,max2,maxn,max);
                MIN3(min1,min2,minn,min);
                max1 = max2; max2 = maxn;
                min1 = min2; min2 = minn;
                if (max < lab->L[i][j])
                    max = lab->L[i][j];
                if (min > lab->L[i][j])
                    min = lab->L[i][j];
                int val = lab->L[i][j] + params->sharpening.amount * diff / 100;
                int newL = CLIP(val);
                // applying halo control
                if (newL > max)
                    newL = max + (newL-max) * scale / 10000;
                else if (newL<min)
                    newL = min - (min-newL) * scale / 10000;
                lab->L[i][j] = newL;
            }
        }
    }
}

}
