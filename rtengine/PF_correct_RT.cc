////////////////////////////////////////////////////////////////
//
//      Chromatic Aberration Auto-correction
//
//      copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: November 24, 2010
// optimized: September 2013, Ingo Weyrich
//
//  PF_correct_RT.cc is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "gauss.h"
#include "improcfun.h"
#include "sleef.c"
#include "mytime.h"
#include "../rtgui/myflatcurve.h"
#include "rt_math.h"
#include "opthelper.h"
#include "median.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace rtengine
{
extern const Settings* settings;

SSEFUNCTION void ImProcFunctions::PF_correct_RT(LabImage * src, LabImage * dst, double radius, int thresh)
{
    const int halfwin = ceil(2 * radius) + 1;

    FlatCurve* chCurve = nullptr;

    if (params->defringe.huecurve.size() && FlatCurveType(params->defringe.huecurve.at(0)) > FCT_Linear) {
        chCurve = new FlatCurve(params->defringe.huecurve);
    }

    // local variables
    const int width = src->W, height = src->H;
    //temporary array to store chromaticity
    float (*fringe);
    fringe = (float (*)) malloc (height * width * sizeof(*fringe));

    LabImage * tmp1;
    tmp1 = new LabImage(width, height);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur (src->a, tmp1->a, src->W, src->H, radius);
        gaussianBlur (src->b, tmp1->b, src->W, src->H, radius);
    }

    float chromave = 0.0f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float chromaChfactor = 1.0f;
#ifdef _OPENMP
        #pragma omp for reduction(+:chromave)
#endif

        for(int i = 0; i < height; i++ ) {
#ifdef __SSE2__

            // vectorized per row precalculation of the atan2 values
            if (chCurve) {
                int k = 0;

                for(; k < width - 3; k += 4) {
                    STVFU(fringe[i * width + k], xatan2f(LVFU(src->b[i][k]), LVFU(src->a[i][k])));
                }

                for(; k < width; k++) {
                    fringe[i * width + k] = xatan2f(src->b[i][k], src->a[i][k]);
                }
            }

#endif // __SSE2__

            for(int j = 0; j < width; j++) {
                if (chCurve) {
#ifdef __SSE2__
                    // use the precalculated atan values
                    float HH = fringe[i * width + j];
#else
                    // no precalculated values without SSE => calculate
                    float HH = xatan2f(src->b[i][j], src->a[i][j]);
#endif
                    float chparam = float((chCurve->getVal((Color::huelab_to_huehsv2(HH))) - 0.5f) * 2.0f); //get C=f(H)

                    if(chparam > 0.f) {
                        chparam /= 2.f;    // reduced action if chparam > 0
                    }

                    chromaChfactor = 1.0f + chparam;
                }

                float chroma = SQR(chromaChfactor * (src->a[i][j] - tmp1->a[i][j])) + SQR(chromaChfactor * (src->b[i][j] - tmp1->b[i][j])); //modulate chroma function hue
                chromave += chroma;
                fringe[i * width + j] = chroma;
            }
        }
    }

    chromave /= (height * width);
    float threshfactor = SQR(thresh / 33.f) * chromave * 5.0f;


// now chromave is calculated, so we postprocess fringe to reduce the number of divisions in future
#ifdef __SSE2__
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        __m128 sumv = F2V( chromave );
        __m128 onev = F2V( 1.0f );
#ifdef _OPENMP
        #pragma omp for nowait
#endif

        for(int j = 0; j < width * height - 3; j += 4) {
            STVFU(fringe[j], onev / (LVFU(fringe[j]) + sumv));
        }

        #pragma omp single

        for(int j = width * height - (width * height) % 4; j < width * height; j++) {
            fringe[j] = 1.f / (fringe[j] + chromave);
        }
    }

#else
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(int j = 0; j < width * height; j++) {
        fringe[j] = 1.f / (fringe[j] + chromave);
    }

#endif

    // because we changed the values of fringe we also have to recalculate threshfactor
    threshfactor = 1.0f / (threshfactor + chromave);

// Issue 1674:
// often, CA isn't evenly distributed, e.g. a lot in contrasty regions and none in the sky.
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
// Issue 1972: Split this loop in three parts to avoid most of the min and max-operations
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for(int i = 0; i < height; i++ ) {
        int j;

        for(j = 0; j < halfwin - 1; j++) {
            tmp1->a[i][j] = src->a[i][j];
            tmp1->b[i][j] = src->b[i][j];

            //test for pixel darker than some fraction of neighborhood ave, near an edge, more saturated than average
            if (fringe[i * width + j] < threshfactor) {
                float atot = 0.f;
                float btot = 0.f;
                float norm = 0.f;
                float wt;

                for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                    for (int j1 = 0; j1 < j + halfwin; j1++) {
                        //neighborhood average of pixels weighted by chrominance
                        wt = fringe[i1 * width + j1];
                        atot += wt * src->a[i1][j1];
                        btot += wt * src->b[i1][j1];
                        norm += wt;
                    }

                tmp1->a[i][j] = atot / norm;
                tmp1->b[i][j] = btot / norm;
            }
        }

        for(; j < width - halfwin + 1; j++) {
            tmp1->a[i][j] = src->a[i][j];
            tmp1->b[i][j] = src->b[i][j];

            //test for pixel darker than some fraction of neighborhood ave, near an edge, more saturated than average
            if (fringe[i * width + j] < threshfactor) {
                float atot = 0.f;
                float btot = 0.f;
                float norm = 0.f;
                float wt;

                for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                    for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                        //neighborhood average of pixels weighted by chrominance
                        wt = fringe[i1 * width + j1];
                        atot += wt * src->a[i1][j1];
                        btot += wt * src->b[i1][j1];
                        norm += wt;
                    }

                tmp1->a[i][j] = atot / norm;
                tmp1->b[i][j] = btot / norm;
            }
        }

        for(; j < width; j++) {
            tmp1->a[i][j] = src->a[i][j];
            tmp1->b[i][j] = src->b[i][j];

            //test for pixel darker than some fraction of neighborhood ave, near an edge, more saturated than average
            if (fringe[i * width + j] < threshfactor) {
                float atot = 0.f;
                float btot = 0.f;
                float norm = 0.f;
                float wt;

                for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                    for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                        //neighborhood average of pixels weighted by chrominance
                        wt = fringe[i1 * width + j1];
                        atot += wt * src->a[i1][j1];
                        btot += wt * src->b[i1][j1];
                        norm += wt;
                    }

                tmp1->a[i][j] = atot / norm;
                tmp1->b[i][j] = btot / norm;
            }
        }
    }//end of ab channel averaging

    if(src != dst)
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int i = 0; i < height; i++ ) {
            for(int j = 0; j < width; j++) {
                dst->L[i][j] = src->L[i][j];
            }
        }

#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(int i = 0; i < height; i++ ) {
        for(int j = 0; j < width; j++) {
            dst->a[i][j] = tmp1->a[i][j];
            dst->b[i][j] = tmp1->b[i][j];
        }
    }


    delete tmp1;

    if(chCurve) {
        delete chCurve;
    }

    free(fringe);
}

SSEFUNCTION void ImProcFunctions::PF_correct_RTcam(CieImage * src, CieImage * dst, double radius, int thresh)
{
    const int halfwin = ceil(2 * radius) + 1;

    FlatCurve* chCurve = nullptr;

    if (params->defringe.huecurve.size() && FlatCurveType(params->defringe.huecurve.at(0)) > FCT_Linear) {
        chCurve = new FlatCurve(params->defringe.huecurve);
    }

    // local variables
    const int width = src->W, height = src->H;
    const float piid = 3.14159265f / 180.f;
    const float eps2 = 0.01f;

    //temporary array to store chromaticity
    float (*fringe);
    fringe = (float (*)) malloc (height * width * sizeof(*fringe));

    float** sraa;
    sraa = new float*[height];

    for (int i = 0; i < height; i++) {
        sraa[i] = new float[width];
    }

    float** srbb;
    srbb = new float*[height];

    for (int i = 0; i < height; i++) {
        srbb[i] = new float[width];
    }

    float** tmaa;
    tmaa = new float*[height];

    for (int i = 0; i < height; i++) {
        tmaa[i] = new float[width];
    }

    float** tmbb;
    tmbb = new float*[height];

    for (int i = 0; i < height; i++) {
        tmbb[i] = new float[width];
    }


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float2 sincosval;
#ifdef __SSE2__
        int j;
        vfloat2 sincosvalv;
        __m128 piidv = F2V(piid);
#endif // __SSE2__
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
#ifdef __SSE2__

            for (j = 0; j < width - 3; j += 4) {
                sincosvalv = xsincosf(piidv * LVFU(src->h_p[i][j]));
                STVFU(sraa[i][j], LVFU(src->C_p[i][j])*sincosvalv.y);
                STVFU(srbb[i][j], LVFU(src->C_p[i][j])*sincosvalv.x);
            }

            for (; j < width; j++) {
                sincosval = xsincosf(piid * src->h_p[i][j]);
                sraa[i][j] = src->C_p[i][j] * sincosval.y;
                srbb[i][j] = src->C_p[i][j] * sincosval.x;
            }

#else

            for (int j = 0; j < width; j++) {
                sincosval = xsincosf(piid * src->h_p[i][j]);
                sraa[i][j] = src->C_p[i][j] * sincosval.y;
                srbb[i][j] = src->C_p[i][j] * sincosval.x;
            }

#endif
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur (sraa, tmaa, src->W, src->H, radius);
        gaussianBlur (srbb, tmbb, src->W, src->H, radius);
    }

    float chromave = 0.0f;

#ifdef __SSE2__

    if( chCurve ) {
// vectorized precalculation of the atan2 values
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            int j;
#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height; i++ )
            {
                for(j = 0; j < width - 3; j += 4) {
                    STVFU(fringe[i * width + j], xatan2f(LVFU(srbb[i][j]), LVFU(sraa[i][j])));
                }

                for(; j < width; j++) {
                    fringe[i * width + j] = xatan2f(srbb[i][j], sraa[i][j]);
                }
            }
        }
    }

#endif

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float chromaChfactor = 1.0f;
#ifdef _OPENMP
        #pragma omp for reduction(+:chromave)
#endif

        for(int i = 0; i < height; i++ ) {
            for(int j = 0; j < width; j++) {
                if (chCurve) {
#ifdef __SSE2__
                    // use the precalculated atan values
                    float HH = fringe[i * width + j];
#else
                    // no precalculated values without SSE => calculate
                    float HH = xatan2f(srbb[i][j], sraa[i][j]);
#endif
                    float chparam = float((chCurve->getVal((Color::huelab_to_huehsv2(HH))) - 0.5f) * 2.0f); //get C=f(H)

                    if(chparam > 0.f) {
                        chparam /= 2.f;    // reduced action if chparam > 0
                    }

                    chromaChfactor = 1.0f + chparam;
                }

                float chroma = SQR(chromaChfactor * (sraa[i][j] - tmaa[i][j])) + SQR(chromaChfactor * (srbb[i][j] - tmbb[i][j])); //modulate chroma function hue
                chromave += chroma;
                fringe[i * width + j] = chroma;
            }
        }
    }

    chromave /= (height * width);
    float threshfactor = SQR(thresh / 33.f) * chromave * 5.0f; // Calculated once to eliminate mult inside the next loop

// now chromave is calculated, so we postprocess fringe to reduce the number of divisions in future
#ifdef __SSE2__
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        __m128 sumv = F2V( chromave + eps2 );
        __m128 onev = F2V( 1.0f );
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int j = 0; j < width * height - 3; j += 4) {
            STVFU(fringe[j], onev / (LVFU(fringe[j]) + sumv));
        }
    }

    for(int j = width * height - (width * height) % 4; j < width * height; j++) {
        fringe[j] = 1.f / (fringe[j] + chromave + eps2);
    }

#else
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(int j = 0; j < width * height; j++) {
        fringe[j] = 1.f / (fringe[j] + chromave + eps2);
    }

#endif

    // because we changed the values of fringe we also have to recalculate threshfactor
    threshfactor = 1.0f / (threshfactor + chromave + eps2);

// Issue 1674:
// often, CA isn't evenly distributed, e.g. a lot in contrasty regions and none in the sky.
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
// Issue 1972: Split this loop in three parts to avoid most of the min and max-operations
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for(int i = 0; i < height; i++ ) {
        int j;

        for(j = 0; j < halfwin - 1; j++) {
            tmaa[i][j] = sraa[i][j];
            tmbb[i][j] = srbb[i][j];

            if (fringe[i * width + j] < threshfactor) {
                float atot = 0.f;
                float btot = 0.f;
                float norm = 0.f;
                float wt;

                for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                    for (int j1 = 0; j1 < j + halfwin; j1++) {
                        //neighborhood average of pixels weighted by chrominance
                        wt = fringe[i1 * width + j1];
                        atot += wt * sraa[i1][j1];
                        btot += wt * srbb[i1][j1];
                        norm += wt;
                    }

                if(norm > 0.f) {
                    tmaa[i][j] = (atot / norm);
                    tmbb[i][j] = (btot / norm);
                }
            }
        }

        for(; j < width - halfwin + 1; j++) {
            tmaa[i][j] = sraa[i][j];
            tmbb[i][j] = srbb[i][j];

            if (fringe[i * width + j] < threshfactor) {
                float atot = 0.f;
                float btot = 0.f;
                float norm = 0.f;
                float wt;

                for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                    for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                        //neighborhood average of pixels weighted by chrominance
                        wt = fringe[i1 * width + j1];
                        atot += wt * sraa[i1][j1];
                        btot += wt * srbb[i1][j1];
                        norm += wt;
                    }

                if(norm > 0.f) {
                    tmaa[i][j] = (atot / norm);
                    tmbb[i][j] = (btot / norm);
                }
            }
        }

        for(; j < width; j++) {
            tmaa[i][j] = sraa[i][j];
            tmbb[i][j] = srbb[i][j];

            if (fringe[i * width + j] < threshfactor) {
                float atot = 0.f;
                float btot = 0.f;
                float norm = 0.f;
                float wt;

                for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                    for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                        //neighborhood average of pixels weighted by chrominance
                        wt = fringe[i1 * width + j1];
                        atot += wt * sraa[i1][j1];
                        btot += wt * srbb[i1][j1];
                        norm += wt;
                    }

                if(norm > 0.f) {
                    tmaa[i][j] = (atot / norm);
                    tmbb[i][j] = (btot / norm);
                }
            }
        }
    } //end of ab channel averaging


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef __SSE2__
        int j;
        __m128 interav, interbv;
        __m128 piidv = F2V(piid);
#endif
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int i = 0; i < height; i++ ) {
#ifdef __SSE2__

            for(j = 0; j < width - 3; j += 4) {
                STVFU(dst->sh_p[i][j], LVFU(src->sh_p[i][j]));
                interav = LVFU(tmaa[i][j]);
                interbv = LVFU(tmbb[i][j]);
                STVFU(dst->h_p[i][j], (xatan2f(interbv, interav)) / piidv);
                STVFU(dst->C_p[i][j], vsqrtf(SQRV(interbv) + SQRV(interav)));
            }

            for(; j < width; j++) {
                dst->sh_p[i][j] = src->sh_p[i][j];
                float intera = tmaa[i][j];
                float interb = tmbb[i][j];
                dst->h_p[i][j] = (xatan2f(interb, intera)) / piid;
                dst->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
            }

#else

            for(int j = 0; j < width; j++) {
                dst->sh_p[i][j] = src->sh_p[i][j];
                float intera = tmaa[i][j];
                float interb = tmbb[i][j];
                dst->h_p[i][j] = (xatan2f(interb, intera)) / piid;
                dst->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
            }

#endif
        }
    }

    for (int i = 0; i < height; i++) {
        delete [] sraa[i];
    }

    delete [] sraa;

    for (int i = 0; i < height; i++) {
        delete [] srbb[i];
    }

    delete [] srbb;

    for (int i = 0; i < height; i++) {
        delete [] tmaa[i];
    }

    delete [] tmaa;

    for (int i = 0; i < height; i++) {
        delete [] tmbb[i];
    }

    delete [] tmbb;

    if(chCurve) {
        delete chCurve;
    }

    free(fringe);
}

SSEFUNCTION void ImProcFunctions::Badpixelscam(CieImage * src, CieImage * dst, double radius, int thresh, int mode, float skinprot, float chrom, int hotbad)
{
    const int halfwin = ceil(2 * radius) + 1;
    MyTime t1, t2;
    t1.set();

    const int width = src->W, height = src->H;
    const float piid = 3.14159265f / 180.f;

    int i1, j1;
    const float eps = 1.0f;
    const float eps2 = 0.01f;

    float** sraa;
    sraa = new float*[height];

    for (int i = 0; i < height; i++) {
        sraa[i] = new float[width];
    }

    float** srbb;
    srbb = new float*[height];

    for (int i = 0; i < height; i++) {
        srbb[i] = new float[width];
    }

    float** tmaa;
    tmaa = new float*[height];

    for (int i = 0; i < height; i++) {
        tmaa[i] = new float[width];
    }

    float** tmbb;
    tmbb = new float*[height];

    for (int i = 0; i < height; i++) {
        tmbb[i] = new float[width];
    }

    float* badpix = (float*)malloc(width * height * sizeof(float));

    float** tmL;
    tmL = new float*[height];

    for (int i = 0; i < height; i++) {
        tmL[i] = new float[width];
    }


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float2 sincosval;
#ifdef __SSE2__
        int j;
        vfloat2 sincosvalv;
        __m128 piidv = F2V(piid);
#endif // __SSE2__
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
#ifdef __SSE2__

            for (j = 0; j < width - 3; j += 4) {
                sincosvalv = xsincosf(piidv * LVFU(src->h_p[i][j]));
                STVFU(sraa[i][j], LVFU(src->C_p[i][j])*sincosvalv.y);
                STVFU(srbb[i][j], LVFU(src->C_p[i][j])*sincosvalv.x);
            }

            for (; j < width; j++) {
                sincosval = xsincosf(piid * src->h_p[i][j]);
                sraa[i][j] = src->C_p[i][j] * sincosval.y;
                srbb[i][j] = src->C_p[i][j] * sincosval.x;
            }

#else

            for (int j = 0; j < width; j++) {
                sincosval = xsincosf(piid * src->h_p[i][j]);
                sraa[i][j] = src->C_p[i][j] * sincosval.y;
                srbb[i][j] = src->C_p[i][j] * sincosval.x;
            }

#endif
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        //chroma a and b
        if(mode == 2) { //choice of gaussian blur
            gaussianBlur (sraa, tmaa, src->W, src->H, radius);
            gaussianBlur (srbb, tmbb, src->W, src->H, radius);
        }

        //luma sh_p
        gaussianBlur (src->sh_p, tmL, src->W, src->H, 2.0);//low value to avoid artifacts
    }

    if(mode == 1) { //choice of median
        #pragma omp parallel
        {
            int ip, in, jp, jn;
            #pragma omp for nowait  //nowait because next loop inside this parallel region is independent on this one

            for (int i = 0; i < height; i++) {
                if (i < 2) {
                    ip = i + 2;
                } else {
                    ip = i - 2;
                }

                if (i > height - 3) {
                    in = i - 2;
                } else {
                    in = i + 2;
                }

                for (int j = 0; j < width; j++) {
                    if (j < 2) {
                        jp = j + 2;
                    } else {
                        jp = j - 2;
                    }

                    if (j > width - 3) {
                        jn = j - 2;
                    } else {
                        jn = j + 2;
                    }

                    tmaa[i][j] = median(sraa[ip][jp], sraa[ip][j], sraa[ip][jn], sraa[i][jp], sraa[i][j], sraa[i][jn], sraa[in][jp], sraa[in][j], sraa[in][jn]);
                }
            }

            #pragma omp for

            for (int i = 0; i < height; i++) {
                if (i < 2) {
                    ip = i + 2;
                } else {
                    ip = i - 2;
                }

                if (i > height - 3) {
                    in = i - 2;
                } else {
                    in = i + 2;
                }

                for (int j = 0; j < width; j++) {
                    if (j < 2) {
                        jp = j + 2;
                    } else {
                        jp = j - 2;
                    }

                    if (j > width - 3) {
                        jn = j - 2;
                    } else {
                        jn = j + 2;
                    }

                    tmbb[i][j] = median(srbb[ip][jp], srbb[ip][j], srbb[ip][jn], srbb[i][jp], srbb[i][j], srbb[i][jn], srbb[in][jp], srbb[in][j], srbb[in][jn]);
                }
            }
        }
    }

//luma badpixels
    const float sh_thr = 4.5f;//low value for luma sh_p to avoid artifacts
    const float shthr = sh_thr / 24.0f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int j;
#ifdef __SSE2__
        __m128 shfabsv, shmedv;
        __m128 shthrv = F2V(shthr);
        __m128 onev = F2V(1.0f);
#endif // __SSE2__
#ifdef _OPENMP
        #pragma omp for private(i1,j1)
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                float shfabs = fabs(src->sh_p[i][j] - tmL[i][j]);
                float shmed = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++ ) {
                        shmed += fabs(src->sh_p[i1][j1] - tmL[i1][j1]);
                    }

                badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
            }

#ifdef __SSE2__

            for (; j < width - 5; j += 4) {
                shfabsv = vabsf(LVFU(src->sh_p[i][j]) - LVFU(tmL[i][j]));
                shmedv = ZEROV;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        shmedv += vabsf(LVFU(src->sh_p[i1][j1]) - LVFU(tmL[i1][j1]));
                    }

                STVFU(badpix[i * width + j], vself(vmaskf_gt(shfabsv, (shmedv - shfabsv)*shthrv), onev, ZEROV));
            }

            for (; j < width - 2; j++) {
                float shfabs = fabs(src->sh_p[i][j] - tmL[i][j]);
                float shmed = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        shmed += fabs(src->sh_p[i1][j1] - tmL[i1][j1]);
                    }

                badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
            }

#else

            for (; j < width - 2; j++) {
                float shfabs = fabs(src->sh_p[i][j] - tmL[i][j]);
                float shmed = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        shmed += fabs(src->sh_p[i1][j1] - tmL[i1][j1]);
                    }

                badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
            }

#endif

            for (; j < width; j++) {
                float shfabs = fabs(src->sh_p[i][j] - tmL[i][j]);
                float shmed = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++ ) {
                        shmed += fabs(src->sh_p[i1][j1] - tmL[i1][j1]);
                    }

                badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
            }
        }
    }


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int j;
#ifdef _OPENMP
        #pragma omp for private(i1,j1) schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                if (!badpix[i * width + j]) {
                    continue;
                }

                float norm = 0.0f;
                float shsum = 0.0f;
                float sum = 0.0f;
                int tot = 0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++ ) {
                        if (i1 == i && j1 == j) {
                            continue;
                        }

                        if (badpix[i1 * width + j1]) {
                            continue;
                        }

                        sum += src->sh_p[i1][j1];
                        tot++;
                        float dirsh = 1.f / (SQR(src->sh_p[i1][j1] - src->sh_p[i][j]) + eps);
                        shsum += dirsh * src->sh_p[i1][j1];
                        norm += dirsh;
                    }

                if (norm > 0.f) {
                    src->sh_p[i][j] = shsum / norm;
                } else {
                    if(tot > 0) {
                        src->sh_p[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width - 2; j++) {
                if (!badpix[i * width + j]) {
                    continue;
                }

                float norm = 0.0f;
                float shsum = 0.0f;
                float sum = 0.0f;
                int tot = 0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        if (i1 == i && j1 == j) {
                            continue;
                        }

                        if (badpix[i1 * width + j1]) {
                            continue;
                        }

                        sum += src->sh_p[i1][j1];
                        tot++;
                        float dirsh = 1.f / (SQR(src->sh_p[i1][j1] - src->sh_p[i][j]) + eps);
                        shsum += dirsh * src->sh_p[i1][j1];
                        norm += dirsh;
                    }

                if (norm > 0.f) {
                    src->sh_p[i][j] = shsum / norm;
                } else {
                    if(tot > 0) {
                        src->sh_p[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width; j++) {
                if (!badpix[i * width + j]) {
                    continue;
                }

                float norm = 0.0f;
                float shsum = 0.0f;
                float sum = 0.0f;
                int tot = 0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++ ) {
                        if (i1 == i && j1 == j) {
                            continue;
                        }

                        if (badpix[i1 * width + j1]) {
                            continue;
                        }

                        sum += src->sh_p[i1][j1];
                        tot++;
                        float dirsh = 1.f / (SQR(src->sh_p[i1][j1] - src->sh_p[i][j]) + eps);
                        shsum += dirsh * src->sh_p[i1][j1];
                        norm += dirsh;
                    }

                if (norm > 0.f) {
                    src->sh_p[i][j] = shsum / norm;
                } else {
                    if(tot > 0) {
                        src->sh_p[i][j] = sum / tot;
                    }
                }
            }
        }
    }
// end luma badpixels


// begin chroma badpixels
    float chrommed = 0.f;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:chrommed)
#endif

    for(int i = 0; i < height; i++ ) {
        for(int j = 0; j < width; j++) {
            float chroma = SQR(sraa[i][j] - tmaa[i][j]) + SQR(srbb[i][j] - tmbb[i][j]);
            chrommed += chroma;
            badpix[i * width + j] = chroma;
        }
    }

    chrommed /= (height * width);
    float threshfactor = (thresh * chrommed) / 33.f;

// now chrommed is calculated, so we postprocess badpix to reduce the number of divisions in future
#ifdef __SSE2__
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int j;
        __m128 sumv = F2V( chrommed + eps2 );
        __m128 onev = F2V( 1.0f );
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int i = 0; i < height; i++) {
            for(j = 0; j < width - 3; j += 4) {
                STVFU(badpix[i * width + j], onev / (LVFU(badpix[i * width + j]) + sumv));
            }

            for(; j < width; j++) {
                badpix[i * width + j] = 1.f / (badpix[i * width + j] + chrommed + eps2);
            }
        }
    }
#else
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(int i = 0; i < height; i++)
        for(int j = 0; j < width; j++) {
            badpix[i * width + j] = 1.f / (badpix[i * width + j] + chrommed + eps2);
        }

#endif

    // because we changed the values of badpix we also have to recalculate threshfactor
    threshfactor = 1.0f / (threshfactor + chrommed + eps2);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int j;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for(int i = 0; i < height; i++ ) {
            for(j = 0; j < halfwin; j++) {
                tmaa[i][j] = sraa[i][j];
                tmbb[i][j] = srbb[i][j];

                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f;
                    float btot = 0.f;
                    float norm = 0.f;
                    float wt;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = 0; j1 < j + halfwin; j1++) {
                            wt = badpix[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }

                    if(norm > 0.f) {
                        tmaa[i][j] = (atot / norm);
                        tmbb[i][j] = (btot / norm);
                    }
                }
            }

            for(; j < width - halfwin; j++) {
                tmaa[i][j] = sraa[i][j];
                tmbb[i][j] = srbb[i][j];

                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f;
                    float btot = 0.f;
                    float norm = 0.f;
                    float wt;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            wt = badpix[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }

                    if(norm > 0.f) {
                        tmaa[i][j] = (atot / norm);
                        tmbb[i][j] = (btot / norm);
                    }
                }
            }

            for(; j < width; j++) {
                tmaa[i][j] = sraa[i][j];
                tmbb[i][j] = srbb[i][j];

                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f;
                    float btot = 0.f;
                    float norm = 0.f;
                    float wt;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                            wt = badpix[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }

                    if(norm > 0.f) {
                        tmaa[i][j] = (atot / norm);
                        tmbb[i][j] = (btot / norm);
                    }
                }
            }
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int i = 0; i < height; i++ ) {
            for(int j = 0; j < width; j++) {
                float intera = tmaa[i][j];
                float interb = tmbb[i][j];
                float CC = sqrt(SQR(interb) + SQR(intera));

                if(hotbad == 0) {
                    if(CC < chrom && skinprot != 0.f) {
                        dst->h_p[i][j] = (xatan2f(interb, intera)) / piid;
                        dst->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
                    }
                } else {
                    dst->h_p[i][j] = (xatan2f(interb, intera)) / piid;
                    dst->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
                }
            }
        }
    }

    if(src != dst) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for(int i = 0; i < height; i++ )
            for(int j = 0; j < width; j++) {
                dst->sh_p[i][j] = src->sh_p[i][j];
            }
    }


    for (int i = 0; i < height; i++) {
        delete [] sraa[i];
    }

    delete [] sraa;

    for (int i = 0; i < height; i++) {
        delete [] srbb[i];
    }

    delete [] srbb;

    for (int i = 0; i < height; i++) {
        delete [] tmaa[i];
    }

    delete [] tmaa;

    for (int i = 0; i < height; i++) {
        delete [] tmbb[i];
    }

    delete [] tmbb;

    for (int i = 0; i < height; i++) {
        delete [] tmL[i];
    }

    delete [] tmL;

    free(badpix);

    t2.set();

    if( settings->verbose ) {
        printf("Ciecam badpixels:- %d usec\n", t2.etime(t1));
    }


}

SSEFUNCTION void ImProcFunctions::BadpixelsLab(LabImage * src, LabImage * dst, double radius, int thresh, int mode, float skinprot, float chrom)
{
    const int halfwin = ceil(2 * radius) + 1;
    MyTime t1, t2;
    t1.set();

    const int width = src->W, height = src->H;

    int i1, j1;
    const float eps = 1.0f;
    const float eps2 = 0.01f;

    float** sraa;
    sraa = new float*[height];

    for (int i = 0; i < height; i++) {
        sraa[i] = new float[width];
    }

    float** srbb;
    srbb = new float*[height];

    for (int i = 0; i < height; i++) {
        srbb[i] = new float[width];
    }

    float** tmaa;
    tmaa = new float*[height];

    for (int i = 0; i < height; i++) {
        tmaa[i] = new float[width];
    }

    float** tmbb;
    tmbb = new float*[height];

    for (int i = 0; i < height; i++) {
        tmbb[i] = new float[width];
    }

    float* badpix = (float*)malloc(width * height * sizeof(float));

    float** tmL;
    tmL = new float*[height];

    for (int i = 0; i < height; i++) {
        tmL[i] = new float[width];
    }


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
//  float2 sincosval;
#ifdef __SSE2__
        int j;
//  vfloat2 sincosvalv;
//  __m128 piidv = F2V(piid);
#endif // __SSE2__
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
#ifdef __SSE2__

            for (j = 0; j < width - 3; j += 4) {
                STVFU(sraa[i][j], LVFU(src->a[i][j]));
                STVFU(srbb[i][j], LVFU(src->b[i][j]));
            }

            for (; j < width; j++) {
                sraa[i][j] = src->a[i][j];
                srbb[i][j] = src->b[i][j];
            }

#else

            for (int j = 0; j < width; j++) {
                sraa[i][j] = src->a[i][j];
                srbb[i][j] = src->b[i][j];
            }

#endif
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        //chroma a and b
        if(mode >= 2) { //choice of gaussian blur
            gaussianBlur (sraa, tmaa, src->W, src->H, radius);
            gaussianBlur (srbb, tmbb, src->W, src->H, radius);
        }

        //luma sh_p
        gaussianBlur (src->L, tmL, src->W, src->H, 2.0);//low value to avoid artifacts
    }

    if(mode == 1) { //choice of median
        #pragma omp parallel
        {
            int ip, in, jp, jn;
            #pragma omp for nowait  //nowait because next loop inside this parallel region is independent on this one

            for (int i = 0; i < height; i++) {
                if (i < 2) {
                    ip = i + 2;
                } else {
                    ip = i - 2;
                }

                if (i > height - 3) {
                    in = i - 2;
                } else {
                    in = i + 2;
                }

                for (int j = 0; j < width; j++) {
                    if (j < 2) {
                        jp = j + 2;
                    } else {
                        jp = j - 2;
                    }

                    if (j > width - 3) {
                        jn = j - 2;
                    } else {
                        jn = j + 2;
                    }

                    tmaa[i][j] = median(sraa[ip][jp], sraa[ip][j], sraa[ip][jn], sraa[i][jp], sraa[i][j], sraa[i][jn], sraa[in][jp], sraa[in][j], sraa[in][jn]);
                }
            }

            #pragma omp for

            for (int i = 0; i < height; i++) {
                if (i < 2) {
                    ip = i + 2;
                } else {
                    ip = i - 2;
                }

                if (i > height - 3) {
                    in = i - 2;
                } else {
                    in = i + 2;
                }

                for (int j = 0; j < width; j++) {
                    if (j < 2) {
                        jp = j + 2;
                    } else {
                        jp = j - 2;
                    }

                    if (j > width - 3) {
                        jn = j - 2;
                    } else {
                        jn = j + 2;
                    }

                    tmbb[i][j] = median(srbb[ip][jp], srbb[ip][j], srbb[ip][jn], srbb[i][jp], srbb[i][j], srbb[i][jn], srbb[in][jp], srbb[in][j], srbb[in][jn]);
                }
            }
        }
    }

//luma badpixels
    const float sh_thr = 4.5f;//low value for luma sh_p to avoid artifacts
    const float shthr = sh_thr / 24.0f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int j;
#ifdef __SSE2__
        __m128 shfabsv, shmedv;
        __m128 shthrv = F2V(shthr);
        __m128 onev = F2V(1.0f);
#endif // __SSE2__
#ifdef _OPENMP
        #pragma omp for private(i1,j1)
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                float shfabs = fabs(src->L[i][j] - tmL[i][j]);
                float shmed = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++ ) {
                        shmed += fabs(src->L[i1][j1] - tmL[i1][j1]);
                    }

                badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
            }

#ifdef __SSE2__

            for (; j < width - 5; j += 4) {
                shfabsv = vabsf(LVFU(src->L[i][j]) - LVFU(tmL[i][j]));
                shmedv = ZEROV;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        shmedv += vabsf(LVFU(src->L[i1][j1]) - LVFU(tmL[i1][j1]));
                    }

                STVFU(badpix[i * width + j], vself(vmaskf_gt(shfabsv, (shmedv - shfabsv)*shthrv), onev, ZEROV));
            }

            for (; j < width - 2; j++) {
                float shfabs = fabs(src->L[i][j] - tmL[i][j]);
                float shmed = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        shmed += fabs(src->L[i1][j1] - tmL[i1][j1]);
                    }

                badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
            }

#else

            for (; j < width - 2; j++) {
                float shfabs = fabs(src->L[i][j] - tmL[i][j]);
                float shmed = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        shmed += fabs(src->L[i1][j1] - tmL[i1][j1]);
                    }

                badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
            }

#endif

            for (; j < width; j++) {
                float shfabs = fabs(src->L[i][j] - tmL[i][j]);
                float shmed = 0.0f;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++ ) {
                        shmed += fabs(src->L[i1][j1] - tmL[i1][j1]);
                    }

                badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
            }
        }
    }


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int j;
#ifdef _OPENMP
        #pragma omp for private(i1,j1) schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                if (!badpix[i * width + j]) {
                    continue;
                }

                float norm = 0.0f;
                float shsum = 0.0f;
                float sum = 0.0f;
                int tot = 0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++ ) {
                        if (i1 == i && j1 == j) {
                            continue;
                        }

                        if (badpix[i1 * width + j1]) {
                            continue;
                        }

                        sum += src->L[i1][j1];
                        tot++;
                        float dirsh = 1.f / (SQR(src->L[i1][j1] - src->L[i][j]) + eps);
                        shsum += dirsh * src->L[i1][j1];
                        norm += dirsh;
                    }

                if (norm > 0.f) {
                    src->L[i][j] = shsum / norm;
                } else {
                    if(tot > 0) {
                        src->L[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width - 2; j++) {
                if (!badpix[i * width + j]) {
                    continue;
                }

                float norm = 0.0f;
                float shsum = 0.0f;
                float sum = 0.0f;
                int tot = 0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        if (i1 == i && j1 == j) {
                            continue;
                        }

                        if (badpix[i1 * width + j1]) {
                            continue;
                        }

                        sum += src->L[i1][j1];
                        tot++;
                        float dirsh = 1.f / (SQR(src->L[i1][j1] - src->L[i][j]) + eps);
                        shsum += dirsh * src->L[i1][j1];
                        norm += dirsh;
                    }

                if (norm > 0.f) {
                    src->L[i][j] = shsum / norm;
                } else {
                    if(tot > 0) {
                        src->L[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width; j++) {
                if (!badpix[i * width + j]) {
                    continue;
                }

                float norm = 0.0f;
                float shsum = 0.0f;
                float sum = 0.0f;
                int tot = 0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++ ) {
                        if (i1 == i && j1 == j) {
                            continue;
                        }

                        if (badpix[i1 * width + j1]) {
                            continue;
                        }

                        sum += src->L[i1][j1];
                        tot++;
                        float dirsh = 1.f / (SQR(src->L[i1][j1] - src->L[i][j]) + eps);
                        shsum += dirsh * src->L[i1][j1];
                        norm += dirsh;
                    }

                if (norm > 0.f) {
                    src->L[i][j] = shsum / norm;
                } else {
                    if(tot > 0) {
                        src->L[i][j] = sum / tot;
                    }
                }
            }
        }
    }
// end luma badpixels

    if(mode == 3) {
// begin chroma badpixels
        float chrommed = 0.f;
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:chrommed)
#endif

        for(int i = 0; i < height; i++ ) {
            for(int j = 0; j < width; j++) {
                float chroma = SQR(sraa[i][j] - tmaa[i][j]) + SQR(srbb[i][j] - tmbb[i][j]);
                chrommed += chroma;
                badpix[i * width + j] = chroma;
            }
        }

        chrommed /= (height * width);
        float threshfactor = (thresh * chrommed) / 33.f;

// now chrommed is calculated, so we postprocess badpix to reduce the number of divisions in future
#ifdef __SSE2__
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            int j;
            __m128 sumv = F2V( chrommed + eps2 );
            __m128 onev = F2V( 1.0f );
#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height; i++) {
                for(j = 0; j < width - 3; j += 4) {
                    STVFU(badpix[i * width + j], onev / (LVFU(badpix[i * width + j]) + sumv));
                }

                for(; j < width; j++) {
                    badpix[i * width + j] = 1.f / (badpix[i * width + j] + chrommed + eps2);
                }
            }
        }
#else
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for(int i = 0; i < height; i++)
            for(int j = 0; j < width; j++) {
                badpix[i * width + j] = 1.f / (badpix[i * width + j] + chrommed + eps2);
            }

#endif

        // because we changed the values of badpix we also have to recalculate threshfactor
        threshfactor = 1.0f / (threshfactor + chrommed + eps2);

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            int j;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for(int i = 0; i < height; i++ ) {
                for(j = 0; j < halfwin; j++) {
                    tmaa[i][j] = sraa[i][j];
                    tmbb[i][j] = srbb[i][j];

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f;
                        float btot = 0.f;
                        float norm = 0.f;
                        float wt;

                        for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                            for (int j1 = 0; j1 < j + halfwin; j1++) {
                                wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if(norm > 0.f) {
                            tmaa[i][j] = (atot / norm);
                            tmbb[i][j] = (btot / norm);
                        }
                    }
                }

                for(; j < width - halfwin; j++) {
                    tmaa[i][j] = sraa[i][j];
                    tmbb[i][j] = srbb[i][j];

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f;
                        float btot = 0.f;
                        float norm = 0.f;
                        float wt;

                        for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                            for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                                wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if(norm > 0.f) {
                            tmaa[i][j] = (atot / norm);
                            tmbb[i][j] = (btot / norm);
                        }
                    }
                }

                for(; j < width; j++) {
                    tmaa[i][j] = sraa[i][j];
                    tmbb[i][j] = srbb[i][j];

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f;
                        float btot = 0.f;
                        float norm = 0.f;
                        float wt;

                        for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                            for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                                wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if(norm > 0.f) {
                            tmaa[i][j] = (atot / norm);
                            tmbb[i][j] = (btot / norm);
                        }
                    }
                }
            }
        }

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height; i++ ) {
                for(int j = 0; j < width; j++) {
                    float intera = tmaa[i][j];
                    float interb = tmbb[i][j];
                    float CC = sqrt(SQR(interb / 327.68) + SQR(intera / 327.68f));

                    if(CC < chrom && skinprot != 0.f) {
                        dst->a[i][j] = intera;
                        dst->b[i][j] = interb;
                    }
                }
            }
        }
    }

    if(src != dst) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for(int i = 0; i < height; i++ )
            for(int j = 0; j < width; j++) {
                dst->L[i][j] = src->L[i][j];
            }
    }


    for (int i = 0; i < height; i++) {
        delete [] sraa[i];
    }

    delete [] sraa;

    for (int i = 0; i < height; i++) {
        delete [] srbb[i];
    }

    delete [] srbb;

    for (int i = 0; i < height; i++) {
        delete [] tmaa[i];
    }

    delete [] tmaa;

    for (int i = 0; i < height; i++) {
        delete [] tmbb[i];
    }

    delete [] tmbb;

    for (int i = 0; i < height; i++) {
        delete [] tmL[i];
    }

    delete [] tmL;

    free(badpix);

    t2.set();

    if( settings->verbose ) {
        printf("Lab artifacts:- %d usec\n", t2.etime(t1));
    }


}

}
