////////////////////////////////////////////////////////////////
//
//          Green Equilibration via directional average
//
//  copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//  optimized for speed 2017 Ingo Weyrich <heckflosse67@gmx.de>
//
//
// code dated: August 25, 2017
//
//  green_equil_RT.cc is free software: you can redistribute it and/or modify
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

#include <cmath>
#include <cstdlib>
#include <ctime>

#include "rt_math.h"
#include "rawimagesource.h"
#include "opthelper.h"

namespace rtengine
{

void RawImageSource::green_equilibrate_global(array2D<float> &rawData)
{
    // global correction
    int ng1 = 0, ng2 = 0;
    double avgg1 = 0., avgg2 = 0.;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+: ng1, ng2, avgg1, avgg2) schedule(dynamic,16)
#endif

    for (int i = border; i < H - border; i++) {
        double avgg = 0.;

        for (int j = border + ((FC(i, border) & 1) ^ 1); j < W - border; j += 2) {
            avgg += rawData[i][j];
        }

        int ng = (W - 2 * border + (FC(i, border) & 1)) / 2;

        if (i & 1) {
            avgg2 += avgg;
            ng2 += ng;
        } else {
            avgg1 += avgg;
            ng1 += ng;
        }
    }

    // Avoid division by zero
    if(ng1 == 0 || avgg1 == 0.0) {
        ng1 = 1;
        avgg1 = 1.0;
    }
    if(ng2 == 0 || avgg2 == 0.0) {
        ng2 = 1;
        avgg2 = 1.0;
    }

    double corrg1 = (avgg1 / ng1 + avgg2 / ng2) / 2.0 / (avgg1 / ng1);
    double corrg2 = (avgg1 / ng1 + avgg2 / ng2) / 2.0 / (avgg2 / ng2);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int i = border; i < H - border; i++) {
        double corrg = (i & 1) ? corrg2 : corrg1;

        for (int j = border + ((FC(i, border) & 1) ^ 1); j < W - border; j += 2) {
            rawData[i][j] *= corrg;
        }
    }
}

//void green_equilibrate()//for dcraw implementation
void RawImageSource::green_equilibrate(const GreenEqulibrateThreshold &thresh, array2D<float> &rawData)
{
    // thresh = threshold for performing green equilibration; max percentage difference of G1 vs G2
    // G1-G2 differences larger than this will be assumed to be Nyquist texture, and left untouched

    int height = H, width = W;

    // local variables
    array2D<float> cfa(width / 2 + (width & 1), height);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int i = 0; i < height; ++i) {
        int j = (FC(i, 0) & 1) ^ 1;
#ifdef __SSE2__

        for (; j < width - 7; j += 8) {
            STVFU(cfa[i][j >> 1], LC2VFU(rawData[i][j]));
        }

#endif

        for (; j < width; j += 2) {
            cfa[i][j >> 1] = rawData[i][j];
        }
    }

    constexpr float eps = 1.f; //tolerance to avoid dividing by zero
    // const float thresh6 = 6 * thresh;
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // Fill G interpolated values with border interpolation and input values

    //int vote1, vote2;
    //int counter, vtest;

    //The green equilibration algorithm starts here
    //now smooth the cfa data
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef __SSE2__
        vfloat zd5v = F2V(0.5f);
        vfloat onev = F2V(1.f);
        // vfloat threshv = F2V(thresh);
        // vfloat thresh6v = F2V(thresh6);
        vfloat epsv = F2V(eps);
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int rr = 4; rr < height - 4; rr++) {
            int cc = 5 - (FC(rr, 2) & 1);
#ifdef __SSE2__

            for (; cc < width - 12; cc += 8) {
                //neighbour checking code from Manuel Llorens Garcia
                vfloat o1_1 = LVFU(cfa[rr - 1][(cc - 1) >> 1]);
                vfloat o1_2 = LVFU(cfa[rr - 1][(cc + 1) >> 1]);
                vfloat o1_3 = LVFU(cfa[rr + 1][(cc - 1) >> 1]);
                vfloat o1_4 = LVFU(cfa[rr + 1][(cc + 1) >> 1]);
                vfloat o2_1 = LVFU(cfa[rr - 2][cc >> 1]);
                vfloat o2_2 = LVFU(cfa[rr + 2][cc >> 1]);
                vfloat o2_3 = LVFU(cfa[rr][(cc >> 1) - 1]);
                vfloat o2_4 = LVFU(cfa[rr][(cc >> 1) + 1]);

                vfloat d1 = (o1_1 + o1_2 + o1_3 + o1_4);
                vfloat d2 = (o2_1 + o2_2 + o2_3 + o2_4);

                vfloat c1 = (vabsf(o1_1 - o1_2) + vabsf(o1_1 - o1_3) + vabsf(o1_1 - o1_4) + vabsf(o1_2 - o1_3) + vabsf(o1_3 - o1_4) + vabsf(o1_2 - o1_4));
                vfloat c2 = (vabsf(o2_1 - o2_2) + vabsf(o2_1 - o2_3) + vabsf(o2_1 - o2_4) + vabsf(o2_2 - o2_3) + vabsf(o2_3 - o2_4) + vabsf(o2_2 - o2_4));

                vfloat tfv;
                for (int k = 0; k < 4; ++k) {
                    tfv[k] = thresh(rr, cc + 2 * k);
                }
                vfloat tf6v = F2V(6.f) * tfv;

                vmask mask1 = vmaskf_lt(c1 + c2, tf6v * vabsf(d1 - d2));

                if (_mm_movemask_ps((vfloat)mask1)) {  // if for any of the 4 pixels the condition is true, do the maths for all 4 pixels and mask the unused out at the end
                    //pixel interpolation
                    vfloat gin = LVFU(cfa[rr][cc >> 1]);

                    vfloat gmp2p2 = gin - LVFU(cfa[rr + 2][(cc >> 1) + 1]);
                    vfloat gmm2m2 = gin - LVFU(cfa[rr - 2][(cc >> 1) - 1]);
                    vfloat gmm2p2 = gin - LVFU(cfa[rr - 2][(cc >> 1) + 1]);
                    vfloat gmp2m2 = gin - LVFU(cfa[rr + 2][(cc >> 1) - 1]);

                    vfloat gse = o1_4 + zd5v * gmp2p2;
                    vfloat gnw = o1_1 + zd5v * gmm2m2;
                    vfloat gne = o1_2 + zd5v * gmm2p2;
                    vfloat gsw = o1_3 + zd5v * gmp2m2;

                    vfloat wtse = onev / (epsv + SQRV(gmp2p2) + SQRV(LVFU(cfa[rr + 3][(cc + 3) >> 1]) - o1_4));
                    vfloat wtnw = onev / (epsv + SQRV(gmm2m2) + SQRV(LVFU(cfa[rr - 3][(cc - 3) >> 1]) - o1_1));
                    vfloat wtne = onev / (epsv + SQRV(gmm2p2) + SQRV(LVFU(cfa[rr - 3][(cc + 3) >> 1]) - o1_2));
                    vfloat wtsw = onev / (epsv + SQRV(gmp2m2) + SQRV(LVFU(cfa[rr + 3][(cc - 3) >> 1]) - o1_3));

                    vfloat ginterp = (gse * wtse + gnw * wtnw + gne * wtne + gsw * wtsw) / (wtse + wtnw + wtne + wtsw);

                    vfloat val = vself(vmaskf_lt(ginterp - gin, tfv * (ginterp + gin)), zd5v * (ginterp + gin), gin);
                    val = vself(mask1, val, gin);
                    STC2VFU(rawData[rr][cc], val);
                }
            }

#endif

            for (; cc < width - 6; cc += 2) {
                //neighbour checking code from Manuel Llorens Garcia
                float o1_1 = cfa[rr - 1][(cc - 1) >> 1];
                float o1_2 = cfa[rr - 1][(cc + 1) >> 1];
                float o1_3 = cfa[rr + 1][(cc - 1) >> 1];
                float o1_4 = cfa[rr + 1][(cc + 1) >> 1];
                float o2_1 = cfa[rr - 2][cc >> 1];
                float o2_2 = cfa[rr + 2][cc >> 1];
                float o2_3 = cfa[rr][(cc - 2) >> 1];
                float o2_4 = cfa[rr][(cc + 2) >> 1];

                float d1 = (o1_1 + o1_2) + (o1_3 + o1_4);
                float d2 = (o2_1 + o2_2) + (o2_3 + o2_4);

                float c1 = (fabs(o1_1 - o1_2) + fabs(o1_1 - o1_3) + fabs(o1_1 - o1_4) + fabs(o1_2 - o1_3) + fabs(o1_3 - o1_4) + fabs(o1_2 - o1_4));
                float c2 = (fabs(o2_1 - o2_2) + fabs(o2_1 - o2_3) + fabs(o2_1 - o2_4) + fabs(o2_2 - o2_3) + fabs(o2_3 - o2_4) + fabs(o2_2 - o2_4));

                float tf = thresh(rr, cc);

                if (c1 + c2 < 6 * tf * fabs(d1 - d2)) {
                    //pixel interpolation
                    float gin = cfa[rr][cc >> 1];

                    float gmp2p2 = gin - cfa[rr + 2][(cc + 2) >> 1];
                    float gmm2m2 = gin - cfa[rr - 2][(cc - 2) >> 1];
                    float gmm2p2 = gin - cfa[rr - 2][(cc + 2) >> 1];
                    float gmp2m2 = gin - cfa[rr + 2][(cc - 2) >> 1];

                    float gse = o1_4 + 0.5f * gmp2p2;
                    float gnw = o1_1 + 0.5f * gmm2m2;
                    float gne = o1_2 + 0.5f * gmm2p2;
                    float gsw = o1_3 + 0.5f * gmp2m2;

                    float wtse = 1.f / (eps + SQR(gmp2p2) + SQR(cfa[rr + 3][(cc + 3) >> 1] - o1_4));
                    float wtnw = 1.f / (eps + SQR(gmm2m2) + SQR(cfa[rr - 3][(cc - 3) >> 1] - o1_1));
                    float wtne = 1.f / (eps + SQR(gmm2p2) + SQR(cfa[rr - 3][(cc + 3) >> 1] - o1_2));
                    float wtsw = 1.f / (eps + SQR(gmp2m2) + SQR(cfa[rr + 3][(cc - 3) >> 1] - o1_3));

                    float ginterp = (gse * wtse + gnw * wtnw + gne * wtne + gsw * wtsw) / (wtse + wtnw + wtne + wtsw);

                    if (ginterp - gin < tf * (ginterp + gin)) {
                        rawData[rr][cc] = 0.5f * (ginterp + gin);
                    }
                }
            }
        }
    }
}
}
