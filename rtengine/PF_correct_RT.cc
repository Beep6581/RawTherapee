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

namespace rtengine
{

// Defringe in Lab mode
void ImProcFunctions::PF_correct_RT(LabImage * lab, double radius, int thresh)
{
    BENCHFUN
    std::unique_ptr<FlatCurve> chCurve;
    if (params->defringe.huecurve.size() && FlatCurveType(params->defringe.huecurve.at(0)) > FCT_Linear) {
        chCurve.reset(new FlatCurve(params->defringe.huecurve));
    }

    const int width = lab->W, height = lab->H;

    // temporary array to store chromaticity
    const std::unique_ptr<float[]> fringe(new float[width * height]);

    JaggedArray<float> tmpa(width, height);
    JaggedArray<float> tmpb(width, height);

    double chromave = 0.0; // use double precision for large summations

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur(lab->a, tmpa, width, height, radius);
        gaussianBlur(lab->b, tmpb, width, height, radius);

#ifdef _OPENMP
        #pragma omp for reduction(+:chromave) schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
#ifdef __SSE2__

            // vectorized per row precalculation of the atan2 values
            if (chCurve) {
                int k = 0;

                for (; k < width - 3; k += 4) {
                    STVFU(fringe[i * width + k], xatan2f(LVFU(lab->b[i][k]), LVFU(lab->a[i][k])));
                }

                for (; k < width; k++) {
                    fringe[i * width + k] = xatan2f(lab->b[i][k], lab->a[i][k]);
                }
            }

#endif

            for (int j = 0; j < width; j++) {
                float chromaChfactor = 1.f;
                if (chCurve) {
#ifdef __SSE2__
                    // use the precalculated atan values
                    const float HH = fringe[i * width + j];
#else
                    // no precalculated values without SSE => calculate
                    const float HH = xatan2f(lab->b[i][j], lab->a[i][j]);
#endif
                    float chparam = chCurve->getVal((Color::huelab_to_huehsv2(HH))) - 0.5f; // get C=f(H)

                    if (chparam < 0.f) {
                        chparam *= 2.f; // increased action if chparam < 0
                    }

                    chromaChfactor = SQR(1.f + chparam);
                }

                const float chroma = chromaChfactor * (SQR(lab->a[i][j] - tmpa[i][j]) + SQR(lab->b[i][j] - tmpb[i][j])); // modulate chroma function hue
                chromave += chroma;
                fringe[i * width + j] = chroma;
            }
        }
    }

    chromave /= height * width;

    if (chromave > 0.0) {
        // now as chromave is calculated, we postprocess fringe to reduce the number of divisions in future
#ifdef _OPENMP
        #pragma omp parallel for simd
#endif

        for (int j = 0; j < width * height; j++) {
            fringe[j] = 1.f / (fringe[j] + chromave);
        }

        const float threshfactor = 1.f / (SQR(thresh / 33.f) * chromave * 5.0f + chromave);
        const int halfwin = std::ceil(2 * radius) + 1;

// Issue 1674:
// often, colour fringe is not evenly distributed, e.g. a lot in contrasty regions and none in the sky.
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
// Issue 1972: Split this loop in three parts to avoid most of the min and max-operations
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            int j = 0;
            for (; j < halfwin - 1; j++) {

                // test for pixel darker than some fraction of neighbourhood ave, near an edge, more saturated than average
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++)
                        for (int j1 = 0; j1 < j + halfwin; j1++) {
                            // neighbourhood average of pixels weighted by chrominance
                            const float wt = fringe[i1 * width + j1];
                            atot += wt * lab->a[i1][j1];
                            btot += wt * lab->b[i1][j1];
                            norm += wt;
                        }

                    lab->a[i][j] = atot / norm;
                    lab->b[i][j] = btot / norm;
                }
            }

            for (; j < width - halfwin + 1; j++) {

                // test for pixel darker than some fraction of neighbourhood ave, near an edge, more saturated than average
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++)
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            // neighbourhood average of pixels weighted by chrominance
                            const float wt = fringe[i1 * width + j1];
                            atot += wt * lab->a[i1][j1];
                            btot += wt * lab->b[i1][j1];
                            norm += wt;
                        }

                    lab->a[i][j] = atot / norm;
                    lab->b[i][j] = btot / norm;
                }
            }

            for (; j < width; j++) {

                // test for pixel darker than some fraction of neighbourhood ave, near an edge, more saturated than average
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++)
                        for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                            // neighbourhood average of pixels weighted by chrominance
                            const float wt = fringe[i1 * width + j1];
                            atot += wt * lab->a[i1][j1];
                            btot += wt * lab->b[i1][j1];
                            norm += wt;
                        }

                    lab->a[i][j] = atot / norm;
                    lab->b[i][j] = btot / norm;
                }
            }
        } // end of ab channel averaging
    }
}

// Defringe in CIECAM02 mode
void ImProcFunctions::PF_correct_RTcam(CieImage * ncie, double radius, int thresh)
{
    BENCHFUN

    std::unique_ptr<FlatCurve> chCurve;

    if (params->defringe.huecurve.size() && FlatCurveType(params->defringe.huecurve.at(0)) > FCT_Linear) {
        chCurve.reset(new FlatCurve(params->defringe.huecurve));
    }

    const int width = ncie->W, height = ncie->H;

    // temporary array to store chromaticity
    const std::unique_ptr<float[]> fringe(new float[width * height]);

    float** const sraa = ncie->h_p; // we use the ncie->h_p buffer to avoid memory allocation/deallocation and reduce memory pressure
    float** const srbb = ncie->C_p; // we use the ncie->C_p buffer to avoid memory allocation/deallocation and reduce memory pressure
    JaggedArray<float> tmaa(width, height);
    JaggedArray<float> tmbb(width, height);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef __SSE2__
        const vfloat piDiv180v = F2V(RT_PI_F_180);
#endif
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < height; i++) {
            int j = 0;
#ifdef __SSE2__

            for (; j < width - 3; j += 4) {
                const vfloat2 sincosvalv = xsincosf(piDiv180v * LVFU(ncie->h_p[i][j]));
                STVFU(sraa[i][j], LVFU(ncie->C_p[i][j]) * sincosvalv.y);
                STVFU(srbb[i][j], LVFU(ncie->C_p[i][j]) * sincosvalv.x);
            }
#endif
            for (; j < width; j++) {
                const float2 sincosval = xsincosf(RT_PI_F_180 * ncie->h_p[i][j]);
                sraa[i][j] = ncie->C_p[i][j] * sincosval.y;
                srbb[i][j] = ncie->C_p[i][j] * sincosval.x;
            }
        }
    }

    double chromave = 0.0; // use double precision for large summations

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur(sraa, tmaa, width, height, radius);
        gaussianBlur(srbb, tmbb, width, height, radius);

        float chromaChfactor = 1.f;
#ifdef _OPENMP
        #pragma omp for reduction(+:chromave) schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
#ifdef __SSE2__
            // vectorized per row precalculation of the atan2 values
            if (chCurve) {
                int j = 0;
                for (; j < width - 3; j += 4) {
                    STVFU(fringe[i * width + j], xatan2f(LVFU(srbb[i][j]), LVFU(sraa[i][j])));
                }

                for (; j < width; j++) {
                    fringe[i * width + j] = xatan2f(srbb[i][j], sraa[i][j]);
                }
            }
#endif

            for (int j = 0; j < width; j++) {
                if (chCurve) {
#ifdef __SSE2__
                    // use the precalculated atan2 values
                    const float HH = fringe[i * width + j];
#else
                    // no precalculated values without SSE => calculate
                    const float HH = xatan2f(srbb[i][j], sraa[i][j]);
#endif
                    float chparam = chCurve->getVal(Color::huelab_to_huehsv2(HH)) - 0.5f; //get C=f(H)

                    if (chparam < 0.f) {
                        chparam *= 2.f;    // increase action if chparam < 0
                    }

                    chromaChfactor = SQR(1.f + chparam);
                }

                const float chroma = chromaChfactor * (SQR(sraa[i][j] - tmaa[i][j]) + SQR(srbb[i][j] - tmbb[i][j])); //modulate chroma function hue
                chromave += chroma;
                fringe[i * width + j] = chroma;
            }
        }
    }

    chromave /= height * width;

    if (chromave > 0.0) {
        // now as chromave is calculated, we postprocess fringe to reduce the number of divisions in future
#ifdef _OPENMP
        #pragma omp parallel for simd
#endif

        for (int j = 0; j < width * height; j++) {
            fringe[j] = 1.f / (fringe[j] + chromave);
        }

        const float threshfactor = 1.f / (SQR(thresh / 33.f) * chromave * 5.0f + chromave);
        const int halfwin = std::ceil(2 * radius) + 1;

// Issue 1674:
// often, colour fringe is not evenly distributed, e.g. a lot in contrasty regions and none in the sky.
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
// Issue 1972: Split this loop in three parts to avoid most of the min and max-operations

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            int j = 0;
            for (; j < halfwin - 1; j++) {
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;
                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++) {
                        for (int j1 = 0; j1 < j + halfwin; j1++) {
                            // neighbourhood average of pixels weighted by chrominance
                            const float wt = fringe[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }
                    }
                    tmaa[i][j] = atot / norm;
                    tmbb[i][j] = btot / norm;
                } else {
                    tmaa[i][j] = sraa[i][j];
                    tmbb[i][j] = srbb[i][j];
                }
            }

            for (; j < width - halfwin + 1; j++) {
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;
                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++) {
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            // neighbourhood average of pixels weighted by chrominance
                            const float wt = fringe[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }
                    }
                    tmaa[i][j] = atot / norm;
                    tmbb[i][j] = btot / norm;
                } else {
                    tmaa[i][j] = sraa[i][j];
                    tmbb[i][j] = srbb[i][j];
               }
            }

            for (; j < width; j++) {
                if (fringe[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f,  norm = 0.f;
                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++) {
                        for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                            // neighbourhood average of pixels weighted by chrominance
                            const float wt = fringe[i1 * width + j1];
                            atot += wt * sraa[i1][j1];
                            btot += wt * srbb[i1][j1];
                            norm += wt;
                        }
                    }
                    tmaa[i][j] = atot / norm;
                    tmbb[i][j] = btot / norm;
                } else {
                    tmaa[i][j] = sraa[i][j];
                    tmbb[i][j] = srbb[i][j];
                }
            }
            j = 0;
#ifdef __SSE2__

            for (; j < width - 3; j += 4) {
                const vfloat interav = LVFU(tmaa[i][j]);
                const vfloat interbv = LVFU(tmbb[i][j]);
                STVFU(ncie->h_p[i][j], xatan2f(interbv, interav) / F2V(RT_PI_F_180));
                STVFU(ncie->C_p[i][j], vsqrtf(SQRV(interbv) + SQRV(interav)));
            }
#endif
            for (; j < width; j++) {
                const float intera = tmaa[i][j];
                const float interb = tmbb[i][j];
                ncie->h_p[i][j] = xatan2f(interb, intera) / RT_PI_F_180;
                ncie->C_p[i][j] = sqrt(SQR(interb) + SQR(intera));
            }
        } // end of ab channel averaging
    }
}

// CIECAM02 hot/bad pixel filter
void ImProcFunctions::Badpixelscam(CieImage * ncie, double radius, int thresh, int mode, float chrom, bool hotbad)
{
    BENCHFUN
    if (mode == 2 && radius < 0.25) { // for gauss sigma less than 0.25 gaussianblur() just calls memcpy => nothing to do here
        return;
    }

    const int width = ncie->W, height = ncie->H;

    constexpr float eps = 1.f;

    JaggedArray<float> tmL(width, height);

    const std::unique_ptr<float[]> badpix(new float[width * height]);

    if (radius >= 0.5) { // for gauss sigma less than 0.25 gaussianblur() just calls memcpy => nothing to do here
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            //luma sh_p
            gaussianBlur(ncie->sh_p, tmL, width, height, radius / 2.0); // low value to avoid artifacts
        }

        //luma badpixels
        constexpr float sh_thr = 4.5f; // low value for luma sh_p to avoid artifacts
        constexpr float shthr = sh_thr / 24.0f; // divide by 24 because we are using a 5x5 grid and centre point is excluded from summation

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            const vfloat shthrv = F2V(shthr);
            const vfloat onev = F2V(1.f);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < height; i++) {
                int j = 0;
                for (; j < 2; j++) {
                    const float shfabs = std::fabs(ncie->sh_p[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = 0; j1 <= j + 2; j1++) {
                            shmed += std::fabs(ncie->sh_p[i1][j1] - tmL[i1][j1]);
                        }
                    }

                    badpix[i * width + j] = shfabs > ((shmed - shfabs) * shthr);
                }

#ifdef __SSE2__

                for (; j < width - 5; j += 4) {
                    const vfloat shfabsv = vabsf(LVFU(ncie->sh_p[i][j]) - LVFU(tmL[i][j]));
                    vfloat shmedv = ZEROV;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            shmedv += vabsf(LVFU(ncie->sh_p[i1][j1]) - LVFU(tmL[i1][j1]));
                        }
                    }

                    STVFU(badpix[i * width + j], vselfzero(vmaskf_gt(shfabsv, (shmedv - shfabsv) * shthrv), onev));
                }
#endif
                for (; j < width - 2; j++) {
                    const float shfabs = std::fabs(ncie->sh_p[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            shmed += std::fabs(ncie->sh_p[i1][j1] - tmL[i1][j1]);
                        }
                    }

                    badpix[i * width + j] = shfabs > ((shmed - shfabs) * shthr);
                }

                for (; j < width; j++) {
                    const float shfabs = std::fabs(ncie->sh_p[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 < width; j1++) {
                            shmed += std::fabs(ncie->sh_p[i1][j1] - tmL[i1][j1]);
                        }
                    }

                    badpix[i * width + j] = shfabs > ((shmed - shfabs) * shthr);
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

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = 0; j1 <= j + 2; j1++) {
                            if (!badpix[i1 * width + j1]) {
                                sum += ncie->sh_p[i1][j1];
                                tot += 1.f;
                                const float dirsh = 1.f / (SQR(ncie->sh_p[i1][j1] - ncie->sh_p[i][j]) + eps);
                                shsum += dirsh * ncie->sh_p[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        ncie->sh_p[i][j] = shsum / norm;
                    } else if (tot > 0.f) {
                        ncie->sh_p[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width - 2; j++) {
                if (badpix[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            if (!badpix[i1 * width + j1]) {
                                sum += ncie->sh_p[i1][j1];
                                tot += 1.f;
                                const float dirsh = 1.f / (SQR(ncie->sh_p[i1][j1] - ncie->sh_p[i][j]) + eps);
                                shsum += dirsh * ncie->sh_p[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        ncie->sh_p[i][j] = shsum / norm;
                    } else if (tot > 0.f) {
                        ncie->sh_p[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width; j++) {
                if (badpix[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 < width; j1++) {
                            if (!badpix[i1 * width + j1]) {
                                sum += ncie->sh_p[i1][j1];
                                tot += 1.f;
                                const float dirsh = 1.f / (SQR(ncie->sh_p[i1][j1] - ncie->sh_p[i][j]) + eps);
                                shsum += dirsh * ncie->sh_p[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        ncie->sh_p[i][j] = shsum / norm;
                    } else if (tot > 0.f) {
                        ncie->sh_p[i][j] = sum / tot;
                    }
                }
            }
        }
    } // end luma badpixels

    if (hotbad) {
        JaggedArray<float> sraa(width, height);
        JaggedArray<float> srbb(width, height);

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {

#ifdef __SSE2__
            const vfloat piDiv180v = F2V(RT_PI_F_180);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < height; i++) {
                int j = 0;
#ifdef __SSE2__

                for (; j < width - 3; j += 4) {
                    const vfloat2 sincosvalv = xsincosf(piDiv180v * LVFU(ncie->h_p[i][j]));
                    STVFU(sraa[i][j], LVFU(ncie->C_p[i][j])*sincosvalv.y);
                    STVFU(srbb[i][j], LVFU(ncie->C_p[i][j])*sincosvalv.x);
                }
#endif
                for (; j < width; j++) {
                    const float2 sincosval = xsincosf(RT_PI_F_180 * ncie->h_p[i][j]);
                    sraa[i][j] = ncie->C_p[i][j] * sincosval.y;
                    srbb[i][j] = ncie->C_p[i][j] * sincosval.x;
                }
            }
        }

        float** const tmaa = tmL; // reuse tmL buffer
        JaggedArray<float> tmbb(width, height);

        if (mode == 2) { // choice of gaussian blur
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                //chroma a and b
                gaussianBlur(sraa, tmaa, width, height, radius);
                gaussianBlur(srbb, tmbb, width, height, radius);
            }

        } else if (mode == 1) { // choice of median
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
#ifdef _OPENMP
                #pragma omp for nowait // nowait because next loop inside this parallel region is independent on this one
#endif

                for (int i = 0; i < height; i++) {
                    const int ip = i < 2 ? i + 2 : i - 2;
                    const int in = i > height - 3 ? i - 2 : i + 2;

                    for (int j = 0; j < width; j++) {
                        const int jp = j < 2 ? j + 2 : j -2;
                        const int jn = j > width - 3 ? j - 2 : j + 2;

                        tmaa[i][j] = median(sraa[ip][jp], sraa[ip][j], sraa[ip][jn], sraa[i][jp], sraa[i][j], sraa[i][jn], sraa[in][jp], sraa[in][j], sraa[in][jn]);
                    }
                }

#ifdef _OPENMP
                #pragma omp for
#endif
                for (int i = 0; i < height; i++) {
                    const int ip = i < 2 ? i + 2 : i - 2;
                    const int in = i > height - 3 ? i - 2 : i + 2;

                    for (int j = 0; j < width; j++) {
                        const int jp = j < 2 ? j + 2 : j -2;
                        const int jn = j > width - 3 ? j - 2 : j + 2;

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

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                const float chroma = SQR(sraa[i][j] - tmaa[i][j]) + SQR(srbb[i][j] - tmbb[i][j]);
                chrommed += chroma;
                badpix[i * width + j] = chroma;
            }
        }

        chrommed /= height * width;

        if (chrommed > 0.0) {
        // now as chrommed is calculated, we postprocess badpix to reduce the number of divisions in future
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
#ifdef __SSE2__
                const vfloat chrommedv = F2V(chrommed);
                const  vfloat onev = F2V(1.f);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

                for (int i = 0; i < height; i++) {
                    int j = 0;
#ifdef __SSE2__
                    for (; j < width - 3; j += 4) {
                        STVFU(badpix[i * width + j], onev / (LVFU(badpix[i * width + j]) + chrommedv));
                    }
#endif
                    for (; j < width; j++) {
                        badpix[i * width + j] = 1.f / (badpix[i * width + j] + chrommed);
                    }
                }
            }

            const float threshfactor = 1.f / ((thresh * chrommed) / 33.f + chrommed);
            const int halfwin = std::ceil(2 * radius) + 1;

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int i = 0; i < height; i++) {
                int j = 0;
                for (; j < halfwin; j++) {

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f, btot = 0.f, norm = 0.f;

                        for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++)
                            for (int j1 = 0; j1 < j + halfwin; j1++) {
                                const float wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if (norm > 0.f) {
                            const float intera = atot / norm;
                            const float interb = btot / norm;
                            const float CC = sqrt(SQR(interb) + SQR(intera));

                            if (CC < chrom) {
                                ncie->h_p[i][j] = xatan2f(interb, intera) / RT_PI_F_180;
                                ncie->C_p[i][j] = CC;
                            }
                        }
                    }
                }

#ifdef __SSE2__
                const vfloat threshfactorv = F2V(threshfactor);
                const vfloat chromv = F2V(chrom);
                const vfloat piDiv180v = F2V(RT_PI_F_180);
                for (; j < width - halfwin - 3; j+=4) {

                    vmask selMask = vmaskf_lt(LVFU(badpix[i * width + j]), threshfactorv);
                    if (_mm_movemask_ps((vfloat)selMask)) {
                        vfloat atotv = ZEROV, btotv = ZEROV, normv = ZEROV;

                        for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++)
                            for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                                const vfloat wtv = LVFU(badpix[i1 * width + j1]);
                                atotv += wtv * LVFU(sraa[i1][j1]);
                                btotv += wtv * LVFU(srbb[i1][j1]);
                                normv += wtv;
                            }

                        selMask = vandm(selMask, vmaskf_gt(normv, ZEROV));
                        if (_mm_movemask_ps((vfloat)selMask)) {
                            const vfloat interav = atotv / normv;
                            const vfloat interbv = btotv / normv;
                            const vfloat CCv = vsqrtf(SQRV(interbv) + SQRV(interav));

                            selMask = vandm(selMask, vmaskf_lt(CCv, chromv));
                            if (_mm_movemask_ps((vfloat)selMask)) {
                                STVFU(ncie->h_p[i][j], vself(selMask, xatan2f(interbv, interav) / piDiv180v, LVFU(ncie->h_p[i][j])));
                                STVFU(ncie->C_p[i][j], vself(selMask, CCv, LVFU(ncie->C_p[i][j])));
                            }
                        }
                    }
                }
#endif
                for (; j < width - halfwin; j++) {

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f, btot = 0.f, norm = 0.f;

                        for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++)
                            for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                                const float wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if (norm > 0.f) {
                            const float intera = atot / norm;
                            const float interb = btot / norm;
                            const float CC = sqrt(SQR(interb) + SQR(intera));

                            if (CC < chrom) {
                                ncie->h_p[i][j] = xatan2f(interb, intera) / RT_PI_F_180;
                                ncie->C_p[i][j] = CC;
                            }
                        }
                    }
                }

                for (; j < width; j++) {

                    if (badpix[i * width + j] < threshfactor) {
                        float atot = 0.f, btot = 0.f, norm = 0.f;

                        for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++)
                            for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                                const float wt = badpix[i1 * width + j1];
                                atot += wt * sraa[i1][j1];
                                btot += wt * srbb[i1][j1];
                                norm += wt;
                            }

                        if (norm > 0.f) {
                            const float intera = atot / norm;
                            const float interb = btot / norm;
                            const float CC = sqrt(SQR(interb) + SQR(intera));

                            if (CC < chrom) {
                                ncie->h_p[i][j] = xatan2f(interb, intera) / RT_PI_F_180;
                                ncie->C_p[i][j] = CC;
                            }
                        }
                    }
                }
            }
        }
    }
}

// CbDL reduce artifacts
void ImProcFunctions::BadpixelsLab(LabImage * lab, double radius, int thresh, float chrom)
{
    BENCHFUN

    if (radius < 0.25) { // for gauss sigma less than 0.25 gaussianblur() just calls memcpy => nothing to do here
        return;
    }

    const int halfwin = std::ceil(2 * radius) + 1;

    const int width = lab->W, height = lab->H;

    constexpr float eps = 1.f;

    JaggedArray<float> tmL(width, height);

    const std::unique_ptr<float[]> badpix(new float[width * height]);

    if (radius >= 0.5) { // for gauss sigma less than 0.25 gaussianblur() just calls memcpy => nothing to do here
        //luma badpixels
        // for bad pixels in L channel we need 0 / != 0 information. Use 1 byte per pixel instead of 4 to reduce memory pressure
        uint8_t *badpixb = reinterpret_cast<uint8_t*>(badpix.get());
        constexpr float sh_thr = 4.5f; // low value for luma sh_p to avoid artifacts
        constexpr float shthr = sh_thr / 24.0f; // divide by 24 because we are using a 5x5 grid and centre point is excluded from summation

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            // blur L channel
            gaussianBlur(lab->L, tmL, width, height, radius / 2.0); // low value to avoid artifacts

#ifdef __SSE2__
            const vfloat shthrv = F2V(shthr);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < height; i++) {
                int j = 0;
                for (; j < 2; j++) {
                    const float shfabs = std::fabs(lab->L[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = 0; j1 <= j + 2; j1++) {
                            shmed += std::fabs(lab->L[i1][j1] - tmL[i1][j1]);
                        }
                    }
                    badpixb[i * width + j] = shfabs > ((shmed - shfabs) * shthr);
                }

#ifdef __SSE2__

                for (; j < width - 5; j += 4) {
                    const vfloat shfabsv = vabsf(LVFU(lab->L[i][j]) - LVFU(tmL[i][j]));
                    vfloat shmedv = ZEROV;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            shmedv += vabsf(LVFU(lab->L[i1][j1]) - LVFU(tmL[i1][j1]));
                        }
                    }
                    uint8_t mask = _mm_movemask_ps((vfloat)vmaskf_gt(shfabsv, (shmedv - shfabsv) * shthrv));
                    badpixb[i * width + j] = mask & 1;
                    badpixb[i * width + j + 1] = mask & 2;
                    badpixb[i * width + j + 2] = mask & 4;
                    badpixb[i * width + j + 3] = mask & 8;
                }
#endif
                for (; j < width - 2; j++) {
                    const float shfabs = std::fabs(lab->L[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            shmed += std::fabs(lab->L[i1][j1] - tmL[i1][j1]);
                        }
                    }
                    badpixb[i * width + j] = shfabs > ((shmed - shfabs) * shthr);
                }

                for (; j < width; j++) {
                    const float shfabs = std::fabs(lab->L[i][j] - tmL[i][j]);
                    float shmed = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 < width; j1++) {
                            shmed += std::fabs(lab->L[i1][j1] - tmL[i1][j1]);
                        }
                    }
                    badpixb[i * width + j] = shfabs > ((shmed - shfabs) * shthr);
                }
            }
        }

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            int j = 0;
            for (; j < 2; j++) {
                if (badpixb[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = 0; j1 <= j + 2; j1++) {
                            if (!badpixb[i1 * width + j1]) {
                                sum += lab->L[i1][j1];
                                tot += 1.f;
                                const float dirsh = 1.f / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps);
                                shsum += dirsh * lab->L[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        lab->L[i][j] = shsum / norm;
                    } else if (tot > 0.f) {
                        lab->L[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width - 2; j++) {
                if (badpixb[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 <= j + 2; j1++) {
                            if (!badpixb[i1 * width + j1]) {
                                sum += lab->L[i1][j1];
                                tot += 1.f;
                                const float dirsh = 1.f / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps);
                                shsum += dirsh * lab->L[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        lab->L[i][j] = shsum / norm;
                    } else if (tot > 0.f) {
                        lab->L[i][j] = sum / tot;
                    }
                }
            }

            for (; j < width; j++) {
                if (badpixb[i * width + j]) {
                    float norm = 0.f, shsum = 0.f, sum = 0.f, tot = 0.f;

                    for (int i1 = std::max(0, i - 2); i1 <= std::min(i + 2, height - 1); i1++) {
                        for (int j1 = j - 2; j1 < width; j1++) {
                            if (!badpixb[i1 * width + j1]) {
                                sum += lab->L[i1][j1];
                                tot += 1.f;
                                const float dirsh = 1.f / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps);
                                shsum += dirsh * lab->L[i1][j1];
                                norm += dirsh;
                            }
                        }
                    }
                    if (norm > 0.f) {
                        lab->L[i][j] = shsum / norm;
                    } else if (tot > 0.f) {
                        lab->L[i][j] = sum / tot;
                    }
                }
            }
        }
    } // end luma badpixels

    float** const tmaa = tmL; // reuse tmL buffer
    JaggedArray<float> tmbb(width, height);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        // blur chroma a and b
        gaussianBlur(lab->a, tmaa, width, height, radius);
        gaussianBlur(lab->b, tmbb, width, height, radius);
    }

    // begin chroma badpixels
    double chrommed = 0.0; // use double precision for large summations

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:chrommed)
#endif

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            const float chroma = SQR(lab->a[i][j] - tmaa[i][j]) + SQR(lab->b[i][j] - tmbb[i][j]);
            chrommed += chroma;
            badpix[i * width + j] = chroma;
        }
    }

    chrommed /= height * width;

    if (chrommed > 0.0) {
    // now as chrommed is calculated, we postprocess badpix to reduce the number of divisions in future

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            const vfloat chrommedv = F2V(chrommed);
            const vfloat onev = F2V(1.f);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < height; i++) {
                int j = 0;
#ifdef __SSE2__
                for (; j < width - 3; j += 4) {
                    STVFU(badpix[i * width + j], onev / (LVFU(badpix[i * width + j]) + chrommedv));
                }
#endif
                for (; j < width; j++) {
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

        for (int i = 0; i < height; i++) {
            int j = 0;
            for (; j < halfwin; j++) {
                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++) {
                        for (int j1 = 0; j1 < j + halfwin; j1++) {
                            const float wt = badpix[i1 * width + j1];
                            atot += wt * lab->a[i1][j1];
                            btot += wt * lab->b[i1][j1];
                            norm += wt;
                        }
                    }
                    if (SQR(atot) + SQR(btot) < chrom * SQR(norm)) {
                        lab->a[i][j] = atot / norm;
                        lab->b[i][j] = btot / norm;
                    }
                }
            }

#ifdef __SSE2__
            const vfloat chromv = F2V(chrom);
            const vfloat threshfactorv = F2V(threshfactor);
            for (; j < width - halfwin - 3; j += 4) {
                vmask selMask = vmaskf_lt(LVFU(badpix[i * width + j]), threshfactorv);
                if (_mm_movemask_ps(reinterpret_cast<vfloat>(selMask))) {
                    vfloat atotv = ZEROV, btotv = ZEROV, normv = ZEROV;

                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++) {
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            const vfloat wtv = LVFU(badpix[i1 * width + j1]);
                            atotv += wtv * LVFU(lab->a[i1][j1]);
                            btotv += wtv * LVFU(lab->b[i1][j1]);
                            normv += wtv;
                        }
                    }
                    selMask = vandm(selMask, vmaskf_lt(SQRV(atotv) + SQR(btotv), chromv * SQRV(normv)));
                    if (_mm_movemask_ps(reinterpret_cast<vfloat>(selMask))) {
                        const vfloat aOrig = LVFU(lab->a[i][j]);
                        const vfloat bOrig = LVFU(lab->b[i][j]);
                        STVFU(lab->a[i][j], vself(selMask, atotv / normv, aOrig));
                        STVFU(lab->b[i][j], vself(selMask, btotv / normv, bOrig));
                    }
                }
            }
#endif
            for (; j < width - halfwin; j++) {

                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++) {
                        for (int j1 = j - halfwin + 1; j1 < j + halfwin; j1++) {
                            const float wt = badpix[i1 * width + j1];
                            atot += wt * lab->a[i1][j1];
                            btot += wt * lab->b[i1][j1];
                            norm += wt;
                        }
                    }
                    if (SQR(atot) + SQR(btot) < chrom * SQR(norm)) {
                        lab->a[i][j] = atot / norm;
                        lab->b[i][j] = btot / norm;
                    }
                }
            }

            for (; j < width; j++) {

                if (badpix[i * width + j] < threshfactor) {
                    float atot = 0.f, btot = 0.f, norm = 0.f;

                    for (int i1 = std::max(0, i - halfwin + 1); i1 < std::min(height, i + halfwin); i1++) {
                        for (int j1 = j - halfwin + 1; j1 < width; j1++) {
                            const float wt = badpix[i1 * width + j1];
                            atot += wt * lab->a[i1][j1];
                            btot += wt * lab->b[i1][j1];
                            norm += wt;
                        }
                    }
                    if (SQR(atot) + SQR(btot) < chrom * SQR(norm)) {
                        lab->a[i][j] = atot / norm;
                        lab->b[i][j] = btot / norm;
                    }
                }
            }
        }
    }
}

}
