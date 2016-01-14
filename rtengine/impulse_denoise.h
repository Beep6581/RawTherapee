/*
 *  This file is part of RawTherapee.
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but widthITheightOUT ANY widthARRANTY; without even the implied warranty of
 *  MERCheightANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  2010 Emil Martinec <ejmartin@uchicago.edu>
 *
 */
#include <cstddef>
#include "rt_math.h"
#include "labimage.h"
#include "improcfun.h"
#include "cieimage.h"
#include "sleef.c"
#include "opthelper.h"
#define BENCHMARK
#include "StopWatch.h"

using namespace std;

namespace rtengine
{

SSEFUNCTION void ImProcFunctions::impulse_nr (LabImage* lab, double thresh)
{
    BENCHFUN
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // impulse noise removal
    // local variables

    int width = lab->W;
    int height = lab->H;

    // buffer for the lowpass image
    float ** lpf = new float *[height];
    lpf[0] = new float [width * height];
    // buffer for the highpass image
    char ** impish = new char *[height];
    impish[0] = new char [width * height];

    for (int i = 1; i < height; i++) {
        lpf[i] = lpf[i - 1] + width;
        impish[i] = impish[i - 1] + width;
    }


    //The cleaning algorithm starts here

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // modified bilateral filter for lowpass image, omitting input pixel; or Gaussian blur

    const float eps = 1.0;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur (lab->L, lpf, width, height, max(2.0, thresh - 1.0));
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    float impthr = max(1.0, 5.5 - thresh);
    float impthrDiv24 = impthr / 24.0f;         //Issue 1671: moved the Division outside the loop, impthr can be optimized out too, but I let in the code at the moment


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int i1, j1, j;
        float hpfabs, hfnbrave;
#ifdef __SSE2__
        vfloat hfnbravev, hpfabsv;
        vfloat impthrDiv24v = F2V( impthrDiv24 );
        vfloat onev = F2V( 1.0f );
#endif
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                hpfabs = fabs(lab->L[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++) {
                        hfnbrave += fabs(lab->L[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }

#ifdef __SSE2__

            for (; j < width - 5; j += 4) {
                hfnbravev = ZEROV;
                hpfabsv = vabsf(LVFU(lab->L[i][j]) - LVFU(lpf[i][j]));

                //block average of high pass data
                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ ) {
                    for (j1 = j - 2; j1 <= j + 2; j1++) {
                        hfnbravev += vabsf(LVFU(lab->L[i1][j1]) - LVFU(lpf[i1][j1]));
                    }
                }

                int mask = _mm_movemask_ps((hfnbravev - hpfabsv) * impthrDiv24v - hpfabsv);
                impish[i][j] = (mask & 1);
                impish[i][j + 1] = ((mask & 2) >> 1);
                impish[i][j + 2] = ((mask & 4) >> 2);
                impish[i][j + 3] = ((mask & 8) >> 3);
            }

#endif

            for (; j < width - 2; j++) {
                hpfabs = fabs(lab->L[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++) {
                        hfnbrave += fabs(lab->L[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }

            for (; j < width; j++) {
                hpfabs = fabs(lab->L[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++) {
                        hfnbrave += fabs(lab->L[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }
        }
    }

//now impulsive values have been identified

// Issue 1671:
// often, noise isn't evenly distributed, e.g. only a few noisy pixels in the bright sky, but many in the dark foreground,
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
// race conditions are avoided by the array impish
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int i1, j1, j;
        float wtdsum[3], dirwt, norm;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1 / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * lab->L[i1][j1];
                        wtdsum[1] += dirwt * lab->a[i1][j1];
                        wtdsum[2] += dirwt * lab->b[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    lab->L[i][j] = wtdsum[0] / norm; //low pass filter
                    lab->a[i][j] = wtdsum[1] / norm; //low pass filter
                    lab->b[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }

            for (; j < width - 2; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1 / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * lab->L[i1][j1];
                        wtdsum[1] += dirwt * lab->a[i1][j1];
                        wtdsum[2] += dirwt * lab->b[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    lab->L[i][j] = wtdsum[0] / norm; //low pass filter
                    lab->a[i][j] = wtdsum[1] / norm; //low pass filter
                    lab->b[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }

            for (; j < width; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1 / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * lab->L[i1][j1];
                        wtdsum[1] += dirwt * lab->a[i1][j1];
                        wtdsum[2] += dirwt * lab->b[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    lab->L[i][j] = wtdsum[0] / norm; //low pass filter
                    lab->a[i][j] = wtdsum[1] / norm; //low pass filter
                    lab->b[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }
        }
    }
//now impulsive values have been corrected

    delete [] lpf[0];
    delete [] impish[0];

}


SSEFUNCTION void ImProcFunctions::impulse_nrcam (CieImage* ncie, double thresh, float **buffers[3])
{
    BENCHFUN
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // impulse noise removal
    // local variables

    int width = ncie->W;
    int height = ncie->H;


    float piid = 3.14159265f / 180.f;

    // buffer for the lowpass image
    float ** lpf = buffers[0];
    // buffer for the highpass image
    float ** impish = buffers[1];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // modified bilateral filter for lowpass image, omitting input pixel; or Gaussian blur



    //The cleaning algorithm starts here

    //rangeblur<unsigned short, unsigned int> (lab->L, lpf, impish /*used as buffer here*/, width, height, thresh, false);
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur (ncie->sh_p, lpf, width, height, max(2.0, thresh - 1.0));
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    float impthr = max(1.0f, 5.0f - (float)thresh);
    float impthrDiv24 = impthr / 24.0f;         //Issue 1671: moved the Division outside the loop, impthr can be optimized out too, but I let in the code at the moment

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int i1, j1, j;
        float hpfabs, hfnbrave;
#ifdef __SSE2__
        vfloat hfnbravev, hpfabsv;
        vfloat impthrDiv24v = F2V( impthrDiv24 );
        vfloat onev = F2V( 1.0f );
#endif
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                hpfabs = fabs(ncie->sh_p[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++) {
                        hfnbrave += fabs(ncie->sh_p[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }

#ifdef __SSE2__

            for (; j < width - 5; j += 4) {
                hpfabsv = vabsf(LVFU(ncie->sh_p[i][j]) - LVFU(lpf[i][j]));
                hfnbravev = ZEROV;

                //block average of high pass data
                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ ) {
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        hfnbravev += vabsf(LVFU(ncie->sh_p[i1][j1]) - LVFU(lpf[i1][j1]));
                    }

                }

                STVFU(impish[i][j], vselfzero(vmaskf_gt(hpfabsv, (hfnbravev - hpfabsv)*impthrDiv24v), onev));
            }

#endif

            for (; j < width - 2; j++) {
                hpfabs = fabs(ncie->sh_p[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        hfnbrave += fabs(ncie->sh_p[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }

            for (; j < width; j++) {
                hpfabs = fabs(ncie->sh_p[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++ ) {
                        hfnbrave += fabs(ncie->sh_p[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }
        }
    }

//now impulsive values have been identified

    const float eps = 1.0f;

    float** sraa = buffers[0]; // we can reuse buffers[0] because lpf is not needed anymore at this point
    float** srbb = buffers[2];

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

#ifdef __SSE2__
        vfloat2 sincosvalv;
        vfloat piidv = F2V( piid );
        vfloat tempv;
#endif
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
            int j = 0;
#ifdef __SSE2__

            for (; j < width - 3; j += 4) {
                sincosvalv = xsincosf(piidv * LVFU(ncie->h_p[i][j]));
                tempv = LVFU(ncie->C_p[i][j]);
                STVFU(sraa[i][j], tempv * sincosvalv.y);
                STVFU(srbb[i][j], tempv * sincosvalv.x);
            }

#endif

            for (; j < width; j++) {
                float2 sincosval = xsincosf(piid * ncie->h_p[i][j]);
                sraa[i][j] = ncie->C_p[i][j] * sincosval.y;
                srbb[i][j] = ncie->C_p[i][j] * sincosval.x;
            }
        }
    }

// Issue 1671:
// often, noise isn't evenly distributed, e.g. only a few noisy pixels in the bright sky, but many in the dark foreground,
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
// race conditions are avoided by the array impish
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int i1, j1, j;
        float wtdsum[3], dirwt, norm;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0f;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1.f / (SQR(ncie->sh_p[i1][j1] - ncie->sh_p[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * ncie->sh_p[i1][j1];
                        wtdsum[1] += dirwt * sraa[i1][j1];
                        wtdsum[2] += dirwt * srbb[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    ncie->sh_p[i][j] = wtdsum[0] / norm; //low pass filter
                    sraa[i][j] = wtdsum[1] / norm; //low pass filter
                    srbb[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }

            for (; j < width - 2; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0f;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1.f / (SQR(ncie->sh_p[i1][j1] - ncie->sh_p[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * ncie->sh_p[i1][j1];
                        wtdsum[1] += dirwt * sraa[i1][j1];
                        wtdsum[2] += dirwt * srbb[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    ncie->sh_p[i][j] = wtdsum[0] / norm; //low pass filter
                    sraa[i][j] = wtdsum[1] / norm; //low pass filter
                    srbb[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }

            for (; j < width; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0f;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1.f / (SQR(ncie->sh_p[i1][j1] - ncie->sh_p[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * ncie->sh_p[i1][j1];
                        wtdsum[1] += dirwt * sraa[i1][j1];
                        wtdsum[2] += dirwt * srbb[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    ncie->sh_p[i][j] = wtdsum[0] / norm; //low pass filter
                    sraa[i][j] = wtdsum[1] / norm; //low pass filter
                    srbb[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }
        }
    }

//now impulsive values have been corrected

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef __SSE2__
        vfloat interav, interbv;
        vfloat piidv = F2V(piid);
#endif // __SSE2__
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int i = 0; i < height; i++ ) {
            int j = 0;
#ifdef __SSE2__

            for(; j < width - 3; j += 4) {
                interav = LVFU(sraa[i][j]);
                interbv = LVFU(srbb[i][j]);
                STVFU(ncie->h_p[i][j], (xatan2f(interbv, interav)) / piidv);
                STVFU(ncie->C_p[i][j], vsqrtf(SQRV(interbv) + SQRV(interav)));
            }

#endif

            for(; j < width; j++) {
                float intera = sraa[i][j];
                float interb = srbb[i][j];
                ncie->h_p[i][j] = (xatan2f(interb, intera)) / piid;
                ncie->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
            }
        }
    }

}


}
