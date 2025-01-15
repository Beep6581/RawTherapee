/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Ingo Weyrich (heckflosse67@gmx.de)
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

//
//   Adaptive Homogeneity-Directed interpolation is based on
//   the work of Keigo Hirakawa, Thomas Parks, and Paul Lee.
//   Optimized for speed and reduced memory usage 2018 Ingo Weyrich
//

#include <cstdint>
#include "rtengine.h"
#include "rawimagesource.h"
#include "rt_math.h"
#include "rtgui/multilangmgr.h"
#include "median.h"
//#define BENCHMARK
#include "StopWatch.h"

namespace
{
unsigned fc(const unsigned int cfa[2][2], int r, int c) {
    return cfa[r & 1][c & 1];
}
}
namespace rtengine
{
#define TS 144
void RawImageSource::ahd_demosaic()
{
    BENCHFUN

    const unsigned int cfa[2][2] = {{FC(0,0), FC(0,1)}, {FC(1,0), FC(1,1)}};
    constexpr int dirs[4] = { -1, 1, -TS, TS };
    float xyz_cam[3][3];
    LUTf cbrt(65536);

    int width = W, height = H;

    constexpr float xyz_rgb[3][3] = {        /* XYZ from RGB */
        { 0.412453f, 0.357580f, 0.180423f },
        { 0.212671f, 0.715160f, 0.072169f },
        { 0.019334f, 0.119193f, 0.950227f }
    };

    constexpr float d65_white[3] = { 0.950456f, 1.f, 1.088754f };

    double progress = 0.0;

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), M("TP_RAW_AHD")));
        plistener->setProgress (progress);
    }

    for (int i = 0; i < 65536; i++) {
        const double r = i / 65535.0;
        cbrt[i] = r > 0.008856 ? std::cbrt(r) : 7.787 * r + 16 / 116.0;
    }

    for (int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            xyz_cam[i][j] = 0;
            for (int k = 0; k < 3; k++) {
                xyz_cam[i][j] += xyz_rgb[i][k] * static_cast<float>(imatrices.rgb_cam[k][j]) / d65_white[i];
            }
        }
    }
    border_interpolate(W, H, 5, rawData, red, green, blue);


#ifdef _OPENMP
#pragma omp parallel
#endif
{
    int progresscounter = 0;
    float *buffer = new float[13 * TS * TS]; /* 1053 kB per core */
    auto rgb  = (float(*)[TS][TS][3]) buffer;
    auto lab  = (float(*)[TS][TS][3])(buffer + 6 * TS * TS);
    auto homo = (uint16_t(*)[TS][TS])(buffer + 12 * TS * TS);

#ifdef _OPENMP
    #pragma omp for collapse(2) schedule(dynamic) nowait
#endif
    for (int top = 2; top < height - 5; top += TS - 6) {
        for (int left = 2; left < width - 5; left += TS - 6) {
            //  Interpolate green horizontally and vertically:
            for (int row = top; row < top + TS && row < height - 2; row++) {
                for (int col = left + (fc(cfa, row, left) & 1); col < std::min(left + TS, width - 2); col += 2) {
                    auto pix = &rawData[row][col];
                    float val0 = 0.25f * ((pix[-1] + pix[0] + pix[1]) * 2
                                  - pix[-2] - pix[2]) ;
                    rgb[0][row - top][col - left][1] = median(val0, pix[-1], pix[1]);
                    float val1 = 0.25f * ((pix[-width] + pix[0] + pix[width]) * 2
                                  - pix[-2 * width] - pix[2 * width]) ;
                    rgb[1][row - top][col - left][1] = median(val1, pix[-width], pix[width]);
                }
            }

            //  Interpolate red and blue, and convert to CIELab:
            for (int d = 0; d < 2; d++)
                for (int row = top + 1; row < top + TS - 1 && row < height - 3; row++) {
                    int cng = fc(cfa, row + 1, fc(cfa, row + 1, 0) & 1);
                    for (int col = left + 1; col < std::min(left + TS - 1, width - 3); col++) {
                        auto pix = &rawData[row][col];
                        auto rix = &rgb[d][row - top][col - left];
                        auto lix = lab[d][row - top][col - left];
                        if (fc(cfa, row, col) == 1) {
                            rix[0][2 - cng] = CLIP(pix[0] + (0.5f * (pix[-1] + pix[1]
                                                       - rix[-1][1] - rix[1][1] ) ));
                            rix[0][cng] = CLIP(pix[0] + (0.5f * (pix[-width] + pix[width]
                                                       - rix[-TS][1] - rix[TS][1])));
                            rix[0][1] = pix[0];
                        } else {
                            rix[0][cng] = CLIP(rix[0][1] + (0.25f * (pix[-width - 1] + pix[-width + 1]
                                                        + pix[+width - 1] + pix[+width + 1]
                                                        - rix[-TS - 1][1] - rix[-TS + 1][1]
                                                        - rix[+TS - 1][1] - rix[+TS + 1][1])));
                            rix[0][2 - cng] = pix[0];
                        }
                        float xyz[3] = {};
                        
                        for(unsigned int c = 0; c < 3; ++c) {
                            xyz[0] += xyz_cam[0][c] * rix[0][c];
                            xyz[1] += xyz_cam[1][c] * rix[0][c];
                            xyz[2] += xyz_cam[2][c] * rix[0][c];
                        }

                        xyz[0] = cbrt[xyz[0]];
                        xyz[1] = cbrt[xyz[1]];
                        xyz[2] = cbrt[xyz[2]];

                        lix[0] = 116.f * xyz[1] - 16.f;
                        lix[1] = 500.f * (xyz[0] - xyz[1]);
                        lix[2] = 200.f * (xyz[1] - xyz[2]);
                    }
                }

            //  Build homogeneity maps from the CIELab images:

            for (int row = top + 2; row < top + TS - 2 && row < height - 4; row++) {
                int tr = row - top;
                float ldiff[2][4], abdiff[2][4];

                for (int col = left + 2, tc = 2; col < left + TS - 2 && col < width - 4; col++, tc++) {
                    for (int d = 0; d < 2; d++) {
                        auto lix = &lab[d][tr][tc];

                        for (int i = 0; i < 4; i++) {
                            ldiff[d][i] = std::fabs(lix[0][0] - lix[dirs[i]][0]);
                            abdiff[d][i] = SQR(lix[0][1] - lix[dirs[i]][1])
                                           + SQR(lix[0][2] - lix[dirs[i]][2]);
                        }
                    }

                    float leps = std::min(std::max(ldiff[0][0], ldiff[0][1]),
                                     std::max(ldiff[1][2], ldiff[1][3]));
                    float abeps = std::min(std::max(abdiff[0][0], abdiff[0][1]),
                                      std::max(abdiff[1][2], abdiff[1][3]));

                    for (int d = 0; d < 2; d++) {
                        homo[d][tr][tc] = 0;
                        for (int i = 0; i < 4; i++) {
                            homo[d][tr][tc] += (ldiff[d][i] <= leps) * (abdiff[d][i] <= abeps);
                        }
                    }
                }
            }

            //  Combine the most homogeneous pixels for the final result:
            for (int row = top + 3; row < top + TS - 3 && row < height - 5; row++) {
                int tr = row - top;

                for (int col = left + 3, tc = 3; col < std::min(left + TS - 3, width - 5); col++, tc++) {
                    uint16_t hm0 = 0, hm1 = 0;
                    for (int i = tr - 1; i <= tr + 1; i++)
                        for (int j = tc - 1; j <= tc + 1; j++) {
                            hm0 += homo[0][i][j];
                            hm1 += homo[1][i][j];
                        }

                    if (hm0 != hm1) {
                        int dir = hm1 > hm0;
                        red[row][col] = rgb[dir][tr][tc][0];
                        green[row][col] = rgb[dir][tr][tc][1];
                        blue[row][col] = rgb[dir][tr][tc][2];
                    } else {
                        red[row][col] = 0.5f * (rgb[0][tr][tc][0] + rgb[1][tr][tc][0]);
                        green[row][col] = 0.5f * (rgb[0][tr][tc][1] + rgb[1][tr][tc][1]);
                        blue[row][col] = 0.5f * (rgb[0][tr][tc][2] + rgb[1][tr][tc][2]);
                    }
                }
            }

            if(plistener) {
                progresscounter++;

                if(progresscounter % 32 == 0) {
#ifdef _OPENMP
                    #pragma omp critical (ahdprogress)
#endif
                    {
                        progress += 32.0 * SQR(TS - 6) / (height * width);
                        progress = std::min(progress, 1.0);
                        plistener->setProgress(progress);
                    }
                }
            }
 
        }
    }
    delete [] buffer;
}
    if(plistener) {
        plistener->setProgress (1.0);
    }

}
#undef TS

}
