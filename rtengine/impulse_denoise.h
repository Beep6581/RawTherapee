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

using namespace std;

namespace rtengine
{

SSEFUNCTION void ImProcFunctions::impulse_nr (LabImage* lab, double thresh)
{

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // impulse noise removal
    // local variables

    int width = lab->W;
    int height = lab->H;

    // buffer for the lowpass image
    float ** lpf = new float *[height];
    // buffer for the highpass image
    float ** impish = new float *[height];

    for (int i = 0; i < height; i++) {
        lpf[i] = new float [width];
        //memset (lpf[i], 0, width*sizeof(float));
        impish[i] = new float [width];
        //memset (impish[i], 0, width*sizeof(unsigned short));
    }


    //The cleaning algorithm starts here

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // modified bilateral filter for lowpass image, omitting input pixel; or Gaussian blur

    const float eps = 1.0;

    //rangeblur<unsigned short, unsigned int> (lab->L, lpf, impish /*used as buffer here*/, width, height, thresh, false);
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        AlignedBufferMP<double> buffer(max(width, height));

        gaussHorizontal<float> (lab->L, lpf, buffer, width, height, max(2.0, thresh - 1.0));
        gaussVertical<float>   (lpf, lpf, buffer, width, height, max(2.0, thresh - 1.0));
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
        __m128 hfnbravev, hpfabsv;
        __m128 impthrDiv24v = _mm_set1_ps( impthrDiv24 );
        __m128 onev = _mm_set1_ps( 1.0f );
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
                hfnbravev = _mm_setzero_ps( );
                hpfabsv = vabsf(LVFU(lab->L[i][j]) - LVFU(lpf[i][j]));

                //block average of high pass data
                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++) {
                        hfnbravev += vabsf(LVFU(lab->L[i1][j1]) - LVFU(lpf[i1][j1]));
                    }

                _mm_storeu_ps(&impish[i][j], vself(vmaskf_gt(hpfabsv, (hfnbravev - hpfabsv)*impthrDiv24v), onev, _mm_setzero_ps()));
            }

            for (; j < width - 2; j++) {
                hpfabs = fabs(lab->L[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++) {
                        hfnbrave += fabs(lab->L[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }

#else

            for (; j < width - 2; j++) {
                hpfabs = fabs(lab->L[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++) {
                        hfnbrave += fabs(lab->L[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }

#endif

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
                        if (i1 == i && j1 == j) {
                            continue;
                        }

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
                        if (i1 == i && j1 == j) {
                            continue;
                        }

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
                        if (i1 == i && j1 == j) {
                            continue;
                        }

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

    for (int i = 0; i < height; i++) {
        delete [] lpf[i];
        delete [] impish[i];
    }

    delete [] lpf;
    delete [] impish;

}


SSEFUNCTION void ImProcFunctions::impulse_nrcam (CieImage* ncie, double thresh, float **buffers[3])
{
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
        AlignedBufferMP<double> buffer(max(width, height));
        gaussHorizontal<float> (ncie->sh_p, lpf, buffer, width, height, max(2.0, thresh - 1.0));
        gaussVertical<float>   (lpf, lpf, buffer, width, height, max(2.0, thresh - 1.0));
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
        __m128 hfnbravev, hpfabsv;
        __m128 impthrDiv24v = _mm_set1_ps( impthrDiv24 );
        __m128 onev = _mm_set1_ps( 1.0f );
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
                hfnbravev = _mm_setzero_ps();

                //block average of high pass data
                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ ) {
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        hfnbravev += vabsf(LVFU(ncie->sh_p[i1][j1]) - LVFU(lpf[i1][j1]));
                    }

                    _mm_storeu_ps(&impish[i][j], vself(vmaskf_gt(hpfabsv, (hfnbravev - hpfabsv)*impthrDiv24v), onev, _mm_setzero_ps()));
                }
            }

            for (; j < width - 2; j++) {
                hpfabs = fabs(ncie->sh_p[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        hfnbrave += fabs(ncie->sh_p[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }

#else

            for (; j < width - 2; j++) {
                hpfabs = fabs(ncie->sh_p[i][j] - lpf[i][j]);

                //block average of high pass data
                for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        hfnbrave += fabs(ncie->sh_p[i1][j1] - lpf[i1][j1]);
                    }

                impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
            }

#endif

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
        int j;
        float2 sincosval;
#ifdef __SSE2__
        vfloat2 sincosvalv;
        __m128 piidv = _mm_set1_ps( piid );
        __m128 tempv;
#endif
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
#ifdef __SSE2__

            for (j = 0; j < width - 3; j += 4) {
                sincosvalv = xsincosf(piidv * LVFU(ncie->h_p[i][j]));
                tempv = LVFU(ncie->C_p[i][j]);
                _mm_storeu_ps(&sraa[i][j], tempv * sincosvalv.y);
                _mm_storeu_ps(&srbb[i][j], tempv * sincosvalv.x);
            }

            for (; j < width; j++) {
                sincosval = xsincosf(piid * ncie->h_p[i][j]);
                sraa[i][j] = ncie->C_p[i][j] * sincosval.y;
                srbb[i][j] = ncie->C_p[i][j] * sincosval.x;
            }

#else

            for (j = 0; j < width; j++) {
                sincosval = xsincosf(piid * ncie->h_p[i][j]);
                sraa[i][j] = ncie->C_p[i][j] * sincosval.y;
                srbb[i][j] = ncie->C_p[i][j] * sincosval.x;
            }

#endif
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
                        if (i1 == i && j1 == j) {
                            continue;
                        }

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
                        if (i1 == i && j1 == j) {
                            continue;
                        }

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
                        if (i1 == i && j1 == j) {
                            continue;
                        }

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
        __m128  interav, interbv;
        __m128 piidv = _mm_set1_ps(piid);
#endif // __SSE2__
        int j;
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int i = 0; i < height; i++ ) {
#ifdef __SSE2__

            for(j = 0; j < width - 3; j += 4) {
                interav = LVFU(sraa[i][j]);
                interbv = LVFU(srbb[i][j]);
                _mm_storeu_ps(&ncie->h_p[i][j], (xatan2f(interbv, interav)) / piidv);
                _mm_storeu_ps(&ncie->C_p[i][j], _mm_sqrt_ps(SQRV(interbv) + SQRV(interav)));
            }

            for(; j < width; j++) {
                float intera = sraa[i][j];
                float interb = srbb[i][j];
                ncie->h_p[i][j] = (xatan2f(interb, intera)) / piid;
                ncie->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
            }

#else

            for(j = 0; j < width; j++) {
                float intera = sraa[i][j];
                float interb = srbb[i][j];
                ncie->h_p[i][j] = (xatan2f(interb, intera)) / piid;
                ncie->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
            }

#endif
        }
    }

}


}
