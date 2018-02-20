////////////////////////////////////////////////////////////////
//
//      Chromatic Aberration Auto-correction
//
//      copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: November 24, 2010
// optimized: September 2013, Ingo Weyrich
// further optimized: February 2018, Ingo Weyrich
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
#include "../rtgui/myflatcurve.h"
#include "rt_math.h"
#include "opthelper.h"
#include "median.h"
#include "jaggedarray.h"
#define BENCHMARK
#include "StopWatch.h"

using namespace std;

namespace rtengine
{

void ImProcFunctions::PF_correct_RT(LabImage * src, double radius, int thresh)
{
    BENCHFUN
    const int halfwin = ceil(2 * radius) + 1;

    std::unique_ptr<FlatCurve> chCurve;
    if (params->defringe.huecurve.size() && FlatCurveType(params->defringe.huecurve.at(0)) > FCT_Linear) {
        chCurve.reset(new FlatCurve(params->defringe.huecurve));
    }

    // local variables
    const int width = src->W, height = src->H;

    //temporary array to store chromaticity
    std::unique_ptr<float[]> fringe(new float[width * height]);

    const JaggedArray<float> tmpa(width, height);
    const JaggedArray<float> tmpb(width, height);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur(src->a, tmpa, width, height, radius);
        gaussianBlur(src->b, tmpb, width, height, radius);
    }

    double chromave = 0.0; // use double precision for large summations

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float chromaChfactor = 1.f;
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

#endif

            for(int j = 0; j < width; j++) {
                if (chCurve) {
#ifdef __SSE2__
                    // use the precalculated atan values
                    float HH = fringe[i * width + j];
#else
                    // no precalculated values without SSE => calculate
                    float HH = xatan2f(src->b[i][j], src->a[i][j]);
#endif
                    float chparam = chCurve->getVal((Color::huelab_to_huehsv2(HH))) - 0.5f; //get C=f(H)

                    if(chparam < 0.f) {
                        chparam *= 2.f;    // increased action if chparam < 0
                    }

                    chromaChfactor = SQR(1.f + chparam);
                }

                float chroma = chromaChfactor * (SQR(src->a[i][j] - tmpa[i][j]) + SQR(src->b[i][j] - tmpb[i][j])); //modulate chroma function hue
                chromave += chroma;
                fringe[i * width + j] = chroma;
            }
        }
    }

    chromave /= (height * width);

    if(chromave > 0.0) {
        // now as chromave is calculated, we postprocess fringe to reduce the number of divisions in future
#ifdef _OPENMP
        #pragma omp parallel for simd
#endif

        for(int j = 0; j < width * height; j++) {
            fringe[j] = 1.f / (fringe[j] + chromave);
        }

        const float threshfactor = 1.f / (SQR(thresh / 33.f) * chromave * 5.0f + chromave);

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
            int j = 0;
            for(; j < halfwin - 1; j++) {

                //test for pixel darker than some fraction of neighbourhood ave, near an edge, more saturated than average
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = 0; j1 < j + halfwin; j1++) {
                            //neighbourhood average of pixels weighted by chrominance
                            float wt = fringe[i1 * width + j1];
                            atot += wt * src->a[i1][j1];
                            btot += wt * src->b[i1][j1];
                            norm += wt;
                        }

                    src->a[i][j] = atot / norm;
                    src->b[i][j] = btot / norm;
                }
            }

            for(; j < width - halfwin + 1; j++) {

                //test for pixel darker than some fraction of neighbourhood ave, near an edge, more saturated than average
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            //neighbourhood average of pixels weighted by chrominance
                            float wt = fringe[i1 * width + j1];
                            atot += wt * src->a[i1][j1];
                            btot += wt * src->b[i1][j1];
                            norm += wt;
                        }

                    src->a[i][j] = atot / norm;
                    src->b[i][j] = btot / norm;
                }
            }

            for(; j < width; j++) {

                //test for pixel darker than some fraction of neighbourhood ave, near an edge, more saturated than average
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                            //neighbourhood average of pixels weighted by chrominance
                            float wt = fringe[i1 * width + j1];
                            atot += wt * src->a[i1][j1];
                            btot += wt * src->b[i1][j1];
                            norm += wt;
                        }

                    src->a[i][j] = atot / norm;
                    src->b[i][j] = btot / norm;
                }
            }
        }//end of ab channel averaging
    }
}

void ImProcFunctions::PF_correct_RTcam(CieImage * src, double radius, int thresh)
{
    BENCHFUN
    const int halfwin = ceil(2 * radius) + 1;

    std::unique_ptr<FlatCurve> chCurve;

    if (params->defringe.huecurve.size() && FlatCurveType(params->defringe.huecurve.at(0)) > FCT_Linear) {
        chCurve.reset(new FlatCurve(params->defringe.huecurve));
    }

    // local variables
    const int width = src->W, height = src->H;

    //temporary array to store chromaticity
    std::unique_ptr<float[]> fringe(new float[width * height]);

    float **sraa = src->h_p; // we use the src->h_p buffer to avoid memory allocation/deallocation and reduce memory pressure
    float **srbb = src->C_p; // we use the src->C_p buffer to avoid memory allocation/deallocation and reduce memory pressure
    const JaggedArray<float> tmaa(width, height);
    const JaggedArray<float> tmbb(width, height);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef __SSE2__
        vfloat piDiv180v = F2V(RT_PI_F_180);
#endif
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
            int j = 0;
#ifdef __SSE2__

            for (; j < width - 3; j += 4) {
                vfloat2 sincosvalv = xsincosf(piDiv180v * LVFU(src->h_p[i][j]));
                STVFU(sraa[i][j], LVFU(src->C_p[i][j]) * sincosvalv.y);
                STVFU(srbb[i][j], LVFU(src->C_p[i][j]) * sincosvalv.x);
            }
#endif
            for (; j < width; j++) {
                float2 sincosval = xsincosf(RT_PI_F_180 * src->h_p[i][j]);
                sraa[i][j] = src->C_p[i][j] * sincosval.y;
                srbb[i][j] = src->C_p[i][j] * sincosval.x;
            }
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur(sraa, tmaa, width, height, radius);
        gaussianBlur(srbb, tmbb, width, height, radius);
    }

#ifdef __SSE2__

    if(chCurve) {
        // vectorized precalculation of the atan2 values
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height; i++ ) {
                int j = 0;
                for(; j < width - 3; j += 4) {
                    STVFU(fringe[i * width + j], xatan2f(LVFU(srbb[i][j]), LVFU(sraa[i][j])));
                }

                for(; j < width; j++) {
                    fringe[i * width + j] = xatan2f(srbb[i][j], sraa[i][j]);
                }
            }
        }
    }

#endif

    double chromave = 0.0; // use double precision for large summations

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float chromaChfactor = 1.f;
#ifdef _OPENMP
        #pragma omp for reduction(+:chromave)
#endif

        for(int i = 0; i < height; i++ ) {
            for(int j = 0; j < width; j++) {
                if (chCurve) {
#ifdef __SSE2__
                    // use the precalculated atan2 values
                    float HH = fringe[i * width + j];
#else
                    // no precalculated values without SSE => calculate
                    float HH = xatan2f(srbb[i][j], sraa[i][j]);
#endif
                    float chparam = chCurve->getVal(Color::huelab_to_huehsv2(HH)) - 0.5f; //get C=f(H)

                    if(chparam < 0.f) {
                        chparam *= 2.f;    // increase action if chparam < 0
                    }

                    chromaChfactor = SQR(1.f + chparam);
                }

                float chroma = chromaChfactor * (SQR(sraa[i][j] - tmaa[i][j]) + SQR(srbb[i][j] - tmbb[i][j])); //modulate chroma function hue
                chromave += chroma;
                fringe[i * width + j] = chroma;
            }
        }
    }

    chromave /= (height * width);

    if(chromave > 0.0) {
        // now as chromave is calculated, we postprocess fringe to reduce the number of divisions in future
#ifdef _OPENMP
        #pragma omp parallel for simd
#endif

        for(int j = 0; j < width * height; j++) {
            fringe[j] = 1.f / (fringe[j] + chromave);
        }

        const float threshfactor = 1.f / (SQR(thresh / 33.f) * chromave * 5.0f + chromave);

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
            int j = 0;
            for(; j < halfwin - 1; j++) {
                tmaa[i][j] = sraa[i][j];
                tmbb[i][j] = srbb[i][j];

                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = 0; j1 < j + halfwin; j1++) {
                            //neighbourhood average of pixels weighted by chrominance
                            float wt = fringe[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }

                    if(norm > 0.f) {
                        tmaa[i][j] = atot / norm;
                        tmbb[i][j] = btot / norm;
                    }
                }
            }

            for(; j < width - halfwin + 1; j++) {
                tmaa[i][j] = sraa[i][j];
                tmbb[i][j] = srbb[i][j];

                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            //neighbourhood average of pixels weighted by chrominance
                            float wt = fringe[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }

                    if(norm > 0.f) {
                        tmaa[i][j] = atot / norm;
                        tmbb[i][j] = btot / norm;
                    }
                }
            }

            for(; j < width; j++) {
                tmaa[i][j] = sraa[i][j];
                tmbb[i][j] = srbb[i][j];

                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f,  norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                        for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                            //neighbourhood average of pixels weighted by chrominance
                            float wt = fringe[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }

                    if(norm > 0.f) {
                        tmaa[i][j] = atot / norm;
                        tmbb[i][j] = btot / norm;
                    }
                }
            }
        } //end of ab channel averaging

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for(int i = 0; i < height; i++ ) {
            int j = 0;
#ifdef __SSE2__

            for(; j < width - 3; j += 4) {
                vfloat interav = LVFU(tmaa[i][j]);
                vfloat interbv = LVFU(tmbb[i][j]);
                STVFU(src->h_p[i][j], xatan2f(interbv, interav) / F2V(RT_PI_F_180));
                STVFU(src->C_p[i][j], vsqrtf(SQRV(interbv) + SQRV(interav)));
            }
#endif
            for(; j < width; j++) {
                float intera = tmaa[i][j];
                float interb = tmbb[i][j];
                src->h_p[i][j] = xatan2f(interb, intera) / RT_PI_F_180;
                src->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
            }
        }
    }
}

void ImProcFunctions::Badpixelscam(CieImage * src, double radius, int thresh, int mode, float chrom, bool hotbad)
{
    BENCHFUN
    if(mode == 2 && radius < 0.25) { // for gauss sigma less than 0.25 gaussianblur() just calls memcpy => nothing to do here
        return;
    }

    const int width = src->W, height = src->H;

    constexpr float eps = 1.f;

    const JaggedArray<float> tmL(width, height);

    std::unique_ptr<float[]> badpix(new float[width * height]);

    if(radius >= 0.5) { // for gauss sigma less than 0.25 gaussianblur() just calls memcpy => nothing to do here
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            //luma sh_p
            gaussianBlur(src->sh_p, tmL, width, height, radius / 2.0);//low value to avoid artifacts
        }

        //luma badpixels
        constexpr float sh_thr = 4.5f; //low value for luma sh_p to avoid artifacts
        constexpr float shthr = sh_thr / 24.0f; // divide by 24 because we are using a 5x5 grid and centre point is excluded from summation

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            vfloat shthrv = F2V(shthr);
            vfloat onev = F2V(1.f);
#endif // __SSE2__
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < height; i++) {
                int j = 0;
                for (; j < 2; j++) {
                    float shfabs = fabs(src->sh_p[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                        for (int j1 = 0; j1 <= j + 2; j1++ ) {
                            shmed += fabs(src->sh_p[i1][j1] - tmL[i1][j1]);
                        }

                    badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
                }

#ifdef __SSE2__

                for (; j < width - 5; j += 4) {
                    vfloat shfabsv = vabsf(LVFU(src->sh_p[i][j]) - LVFU(tmL[i][j]));
                    vfloat shmedv = ZEROV;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                        for (int j1 = j - 2; j1 <= j + 2; j1++ ) {
                            shmedv += vabsf(LVFU(src->sh_p[i1][j1]) - LVFU(tmL[i1][j1]));
                        }

                    STVFU(badpix[i * width + j], vselfzero(vmaskf_gt(shfabsv, (shmedv - shfabsv) * shthrv), onev));
                }
#endif
                for (; j < width - 2; j++) {
                    float shfabs = fabs(src->sh_p[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                        for (int j1 = j - 2; j1 <= j + 2; j1++ ) {
                            shmed += fabs(src->sh_p[i1][j1] - tmL[i1][j1]);
                        }

                    badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
                }

                for (; j < width; j++) {
                    float shfabs = fabs(src->sh_p[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                        for (int j1 = j - 2; j1 < width; j1++ ) {
                            shmed += fabs(src->sh_p[i1][j1] - tmL[i1][j1]);
                        }

                    badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
                }
            }
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            int j = 0;
            for (; j < 2; j++) {
                if (badpix[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ ) {
                        for (int j1 = 0; j1 <= j + 2; j1++ ) {
                            if (!badpix[i1 * width + j1]) {
                                sum += src->sh_p[i1][j1];
                                tot += 1.f;
                                float dirsh = 1.f / (SQR(src->sh_p[i1][j1] - src->sh_p[i][j]) + eps);
                                shsum += dirsh * src->sh_p[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        src->sh_p[i][j] = shsum / norm;
                    } else if (tot > 0.f) {
                        src->sh_p[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width - 2; j++) {
                if (badpix[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ ) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++ ) {
                            if (!badpix[i1 * width + j1]) {
                                sum += src->sh_p[i1][j1];
                                tot += 1.f;
                                float dirsh = 1.f / (SQR(src->sh_p[i1][j1] - src->sh_p[i][j]) + eps);
                                shsum += dirsh * src->sh_p[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        src->sh_p[i][j] = shsum / norm;
                    } else if(tot > 0.f) {
                        src->sh_p[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width; j++) {
                if (badpix[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ ) {
                        for (int j1 = j - 2; j1 < width; j1++ ) {
                            if (!badpix[i1 * width + j1]) {
                                sum += src->sh_p[i1][j1];
                                tot += 1.f;
                                float dirsh = 1.f / (SQR(src->sh_p[i1][j1] - src->sh_p[i][j]) + eps);
                                shsum += dirsh * src->sh_p[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        src->sh_p[i][j] = shsum / norm;
                    } else if(tot > 0.f) {
                        src->sh_p[i][j] = sum / tot;
                    }
                }
            }
        }
    }

// end luma badpixels

    if(hotbad) {

        const JaggedArray<float> sraa(width, height);
        const JaggedArray<float> srbb(width, height);

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {

#ifdef __SSE2__
            vfloat piDiv180v = F2V(RT_PI_F_180);
#endif // __SSE2__
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < height; i++) {
                int j = 0;
#ifdef __SSE2__

                for (; j < width - 3; j += 4) {
                    vfloat2 sincosvalv = xsincosf(piDiv180v * LVFU(src->h_p[i][j]));
                    STVFU(sraa[i][j], LVFU(src->C_p[i][j])*sincosvalv.y);
                    STVFU(srbb[i][j], LVFU(src->C_p[i][j])*sincosvalv.x);
                }
#endif
                for (; j < width; j++) {
                    float2 sincosval = xsincosf(RT_PI_F_180 * src->h_p[i][j]);
                    sraa[i][j] = src->C_p[i][j] * sincosval.y;
                    srbb[i][j] = src->C_p[i][j] * sincosval.x;
                }
            }
        }

        float ** tmaa = tmL; // reuse tmL buffer
        const JaggedArray<float> tmbb(width, height);

        if(mode == 2) { //choice of gaussian blur

#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                //chroma a and b
                gaussianBlur(sraa, tmaa, width, height, radius);
                gaussianBlur(srbb, tmbb, width, height, radius);
            }

        } else if(mode == 1) { //choice of median
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
#ifdef _OPENMP
                #pragma omp for nowait  //nowait because next loop inside this parallel region is independent on this one
#endif

                for (int i = 0; i < height; i++) {
                    int ip = i < 2 ? i + 2 : i -2;
                    int in = i > height - 3 ? i - 2 : i + 2;

                    for (int j = 0; j < width; j++) {
                        int jp = j < 2 ? j + 2 : j -2;
                        int jn = j > width - 3 ? j - 2 : j + 2;

                        tmaa[i][j] = median(sraa[ip][jp], sraa[ip][j], sraa[ip][jn], sraa[i][jp], sraa[i][j], sraa[i][jn], sraa[in][jp], sraa[in][j], sraa[in][jn]);
                    }
                }

#ifdef _OPENMP
                #pragma omp for
#endif
                for (int i = 0; i < height; i++) {
                    int ip = i < 2 ? i + 2 : i -2;
                    int in = i > height - 3 ? i - 2 : i + 2;

                    for (int j = 0; j < width; j++) {
                        int jp = j < 2 ? j + 2 : j -2;
                        int jn = j > width - 3 ? j - 2 : j + 2;

                        tmbb[i][j] = median(srbb[ip][jp], srbb[ip][j], srbb[ip][jn], srbb[i][jp], srbb[i][j], srbb[i][jn], srbb[in][jp], srbb[in][j], srbb[in][jn]);
                    }
                }
            }
        }

        // begin chroma badpixels
        double chrommed = 0.0; // use double precision for large summations
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

        if(chrommed > 0.0) {
        // now as chrommed is calculated, we postprocess badpix to reduce the number of divisions in future
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
#ifdef __SSE2__
                vfloat chrommedv = F2V(chrommed);
                vfloat onev = F2V(1.f);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

                for(int i = 0; i < height; i++) {
                    int j = 0;
#ifdef __SSE2__
                    for(; j < width - 3; j += 4) {
                        STVFU(badpix[i * width + j], onev / (LVFU(badpix[i * width + j]) + chrommedv));
                    }
#endif
                    for(; j < width; j++) {
                        badpix[i * width + j] = 1.f / (badpix[i * width + j] + chrommed);
                    }
                }
            }

            const float threshfactor = 1.f / ((thresh * chrommed) / 33.f + chrommed);
            const int halfwin = ceil(2 * radius) + 1;

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for(int i = 0; i < height; i++ ) {
                int j = 0;
                for(; j < halfwin; j++) {

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f, btot = 0.f, norm = 0.f;

                        for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                            for (int j1 = 0; j1 < j + halfwin; j1++) {
                                float wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if(norm > 0.f) {
                            const float intera = atot / norm;
                            const float interb = btot / norm;
                            const float CC = sqrt(SQR(interb) + SQR(intera));

                            if(CC < chrom) {
                                src->h_p[i][j] = xatan2f(interb, intera) / RT_PI_F_180;
                                src->C_p[i][j] = CC;
                            }
                        }
                    }
                }

#ifdef __SSE2__
                vfloat threshfactorv = F2V(threshfactor);
                vfloat chromv = F2V(chrom);
                vfloat piDiv180v = F2V(RT_PI_F_180);
                for(; j < width - halfwin - 3; j+=4) {

                    vmask selMask = vmaskf_lt(LVFU(badpix[i * width + j]), threshfactorv);
                    if(_mm_movemask_ps((vfloat)selMask)) {
                        vfloat atotv = ZEROV, btotv = ZEROV, normv = ZEROV;

                        for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                            for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                                vfloat wtv = LVFU(badpix[i1 * width + j1]);
                                atotv += wtv * LVFU(sraa[i1][j1]);
                                btotv += wtv * LVFU(srbb[i1][j1]);
                                normv += wtv;
                            }

                        selMask = vandm(selMask, vmaskf_gt(normv, ZEROV));
                        if(_mm_movemask_ps((vfloat)selMask)) {
                            vfloat interav = atotv / normv;
                            vfloat interbv = btotv / normv;
                            vfloat CCv = vsqrtf(SQRV(interbv) + SQRV(interav));

                            selMask = vandm(selMask, vmaskf_lt(CCv, chromv));
                            if(_mm_movemask_ps((vfloat)selMask)) {
                                STVFU(src->h_p[i][j], vself(selMask, xatan2f(interbv, interav) / piDiv180v, LVFU(src->h_p[i][j])));
                                STVFU(src->C_p[i][j], vself(selMask, CCv, LVFU(src->C_p[i][j])));
                            }
                        }
                    }
                }
#endif
                for(; j < width - halfwin; j++) {

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f, btot = 0.f, norm = 0.f;

                        for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                            for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                                float wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if(norm > 0.f) {
                            const float intera = atot / norm;
                            const float interb = btot / norm;
                            const float CC = sqrt(SQR(interb) + SQR(intera));

                            if(CC < chrom) {
                                src->h_p[i][j] = xatan2f(interb, intera) / RT_PI_F_180;
                                src->C_p[i][j] = CC;
                            }
                        }
                    }
                }

                for(; j < width; j++) {

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f, btot = 0.f, norm = 0.f;

                        for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++)
                            for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                                float wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if(norm > 0.f) {
                            const float intera = atot / norm;
                            const float interb = btot / norm;
                            const float CC = sqrt(SQR(interb) + SQR(intera));

                            if(CC < chrom) {
                                src->h_p[i][j] = xatan2f(interb, intera) / RT_PI_F_180;
                                src->C_p[i][j] = CC;
                            }
                        }
                    }
                }
            }
        }
    }
}

void ImProcFunctions::BadpixelsLab(LabImage * src, double radius, int thresh, float chrom)
{
    BENCHFUN

    if(radius < 0.25) { // for gauss sigma less than 0.25 gaussianblur() just calls memcpy => nothing to do here
        return;
    }

    const int halfwin = ceil(2 * radius) + 1;

    const int width = src->W, height = src->H;

    constexpr float eps = 1.f;

    const JaggedArray<float> tmL(width, height);

    std::unique_ptr<float[]> badpix(new float[width * height]);

    if(radius >= 0.5) { // for gauss sigma less than 0.25 gaussianblur() just calls memcpy => nothing to do here

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            // blur L channel
            gaussianBlur(src->L, tmL, width, height, radius / 2.0);//low value to avoid artifacts
        }

        //luma badpixels
        constexpr float sh_thr = 4.5f; //low value for luma sh_p to avoid artifacts
        constexpr float shthr = sh_thr / 24.0f; // divide by 24 because we are using a 5x5 grid and centre point is excluded from summation

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            vfloat shthrv = F2V(shthr);
            vfloat onev = F2V(1.f);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < height; i++) {
                int j = 0;
                for (; j < 2; j++) {
                    float shfabs = fabs(src->L[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++) {
                        for (int j1 = 0; j1 <= j + 2; j1++) {
                            shmed += fabs(src->L[i1][j1] - tmL[i1][j1]);
                        }
                    }
                    badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
                }

#ifdef __SSE2__

                for (; j < width - 5; j += 4) {
                    vfloat shfabsv = vabsf(LVFU(src->L[i][j]) - LVFU(tmL[i][j]));
                    vfloat shmedv = ZEROV;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            shmedv += vabsf(LVFU(src->L[i1][j1]) - LVFU(tmL[i1][j1]));
                        }
                    }
                    STVFU(badpix[i * width + j], vselfzero(vmaskf_gt(shfabsv, (shmedv - shfabsv) * shthrv), onev));
                }
#endif
                for (; j < width - 2; j++) {
                    float shfabs = fabs(src->L[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            shmed += fabs(src->L[i1][j1] - tmL[i1][j1]);
                        }
                    }
                    badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
                }

                for (; j < width; j++) {
                    float shfabs = fabs(src->L[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 < width; j1++) {
                            shmed += fabs(src->L[i1][j1] - tmL[i1][j1]);
                        }
                    }
                    badpix[i * width + j] = (shfabs > ((shmed - shfabs) * shthr));
                }
            }
        }

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            int j = 0;
            for (; j < 2; j++) {
                if (badpix[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++) {
                        for (int j1 = 0; j1 <= j + 2; j1++) {
                            if (!badpix[i1 * width + j1]) {
                                sum += src->L[i1][j1];
                                tot += 1.f;
                                float dirsh = 1.f / (SQR(src->L[i1][j1] - src->L[i][j]) + eps);
                                shsum += dirsh * src->L[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        src->L[i][j] = shsum / norm;
                    } else if(tot > 0.f) {
                        src->L[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width - 2; j++) {
                if (badpix[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            if (!badpix[i1 * width + j1]) {
                                sum += src->L[i1][j1];
                                tot += 1.f;
                                float dirsh = 1.f / (SQR(src->L[i1][j1] - src->L[i][j]) + eps);
                                shsum += dirsh * src->L[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        src->L[i][j] = shsum / norm;
                    } else if(tot > 0.f) {
                        src->L[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width; j++) {
                if (badpix[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 < width; j1++) {
                            if (!badpix[i1 * width + j1]) {
                                sum += src->L[i1][j1];
                                tot += 1.f;
                                float dirsh = 1.f / (SQR(src->L[i1][j1] - src->L[i][j]) + eps);
                                shsum += dirsh * src->L[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        src->L[i][j] = shsum / norm;
                    } else if(tot > 0.f) {
                        src->L[i][j] = sum / tot;
                    }
                }
            }
        }
    }

    // end luma badpixels

    float ** tmaa = tmL; // reuse tmL buffer
    const JaggedArray<float> tmbb(width, height);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        // blur chroma a and b
        gaussianBlur(src->a, tmaa, width, height, radius);
        gaussianBlur(src->b, tmbb, width, height, radius);
    }

// begin chroma badpixels
    double chrommed = 0.0; // use double precision for large summations

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:chrommed)
#endif

    for(int i = 0; i < height; i++ ) {
        for(int j = 0; j < width; j++) {
            float chroma = SQR(src->a[i][j] - tmaa[i][j]) + SQR(src->b[i][j] - tmbb[i][j]);
            chrommed += chroma;
            badpix[i * width + j] = chroma;
        }
    }

    chrommed /= (height * width);
    if(chrommed > 0.0) {

    // now as chrommed is calculated, we postprocess badpix to reduce the number of divisions in future

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            vfloat chrommedv = F2V(chrommed);
            vfloat onev = F2V(1.f);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height; i++) {
                int j = 0;
#ifdef __SSE2__
                for(; j < width - 3; j += 4) {
                    STVFU(badpix[i * width + j], onev / (LVFU(badpix[i * width + j]) + chrommedv));
                }
#endif
                for(; j < width; j++) {
                    badpix[i * width + j] = 1.f / (badpix[i * width + j] + chrommed);
                }
            }
        }

        const float threshfactor = 1.f / ((thresh * chrommed) / 33.f + chrommed);

        chrom *= 327.68f;
        chrom *= chrom;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for(int i = 0; i < height; i++ ) {
            int j = 0;
            for(; j < halfwin; j++) {

                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++) {
                        for (int j1 = 0; j1 < j + halfwin; j1++) {
                            float wt = badpix[i1 * width + j1];
                            atot += wt * src->a[i1][j1];
                            btot += wt * src->b[i1][j1];
                            norm += wt;
                        }
                    }
                    if(SQR(atot) + SQR(btot) < chrom * SQR(norm)) {
                        src->a[i][j] = atot / norm;
                        src->b[i][j] = btot / norm;
                    }
                }
            }

#ifdef __SSE2__
            vfloat chromv = F2V(chrom);
            vfloat threshfactorv = F2V(threshfactor);
            for(; j < width - halfwin - 3; j+=4) {

                vmask selMask = vmaskf_lt(LVFU(badpix[i * width + j]), threshfactorv);
                if (_mm_movemask_ps((vfloat)selMask)) {
                    vfloat atotv = ZEROV, btotv = ZEROV, normv = ZEROV;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++) {
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            vfloat wtv = LVFU(badpix[i1 * width + j1]);
                            atotv += wtv * LVFU(src->a[i1][j1]);
                            btotv += wtv * LVFU(src->b[i1][j1]);
                            normv += wtv;
                        }
                    }
                    selMask = vandm(selMask, vmaskf_lt(SQRV(atotv) + SQR(btotv), chromv * SQRV(normv)));
                    if(_mm_movemask_ps((vfloat)selMask)) {
                        vfloat aOrig = LVFU(src->a[i][j]);
                        vfloat bOrig = LVFU(src->b[i][j]);
                        STVFU(src->a[i][j], vself(selMask, atotv / normv, aOrig));
                        STVFU(src->b[i][j], vself(selMask, btotv / normv, bOrig));
                    }
                }
            }
#endif
            for(; j < width - halfwin; j++) {

                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++) {
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            float wt = badpix[i1 * width + j1];
                            atot += wt * src->a[i1][j1];
                            btot += wt * src->b[i1][j1];
                            norm += wt;
                        }
                    }
                    if(SQR(atot) + SQR(btot) < chrom * SQR(norm)) {
                        src->a[i][j] = atot / norm;
                        src->b[i][j] = btot / norm;
                    }
                }
            }

            for(; j < width; j++) {

                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = max(0, i - halfwin + 1); i1 < min(height, i + halfwin); i1++) {
                        for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                            float wt = badpix[i1 * width + j1];
                            atot += wt * src->a[i1][j1];
                            btot += wt * src->b[i1][j1];
                            norm += wt;
                        }
                    }
                    if(SQR(atot) + SQR(btot) < chrom * SQR(norm)) {
                        src->a[i][j] = atot / norm;
                        src->b[i][j] = btot / norm;
                    }
                }
            }
        }
    }
}

}
