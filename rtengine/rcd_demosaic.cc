/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017-2020 Luis Sanz Rodriguez (luis.sanz.rodriguez(at)gmail(dot)com) and Ingo Weyrich (heckflosse67@gmx.de)
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
#include <cmath>

#include "rawimagesource.h"
#include "rt_math.h"
#include "rtgui/multilangmgr.h"
#include "StopWatch.h"

using namespace std;

namespace
{
unsigned fc(const unsigned int cfa[2][2], int r, int c) {
    return cfa[r & 1][c & 1];
}
}

namespace rtengine
{

/*
* RATIO CORRECTED DEMOSAICING
* Luis Sanz Rodriguez (luis.sanz.rodriguez(at)gmail(dot)com)
*
* Release 2.3 @ 171125
*
* Original code from https://github.com/LuisSR/RCD-Demosaicing
* Licensed under the GNU GPL version 3
*/

// Tiled version by Ingo Weyrich (heckflosse67@gmx.de)
// Luis Sanz Rodriguez significantly optimised the v 2.3 code and simplified the directional
// coefficients in an exact, shorter and more performant formula.
// In cooperation with Hanno Schwalm (hanno@schwalm-bremen.de) and Luis Sanz Rodriguez this has been tuned for performance.

void RawImageSource::rcd_demosaic(size_t chunkSize, bool measure)
{
    // Test for RGB cfa
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            if (FC(i, j) == 3) {
                // avoid crash
                std::cout << "rcd_demosaic supports only RGB Colour filter arrays. Falling back to igv_interpolate" << std::endl;
                igv_interpolate(W, H);
                return;
            }
        }
    }

    std::unique_ptr<StopWatch> stop;
    if (measure) {
        std::cout << "Demosaicing " << W << "x" << H << " image using rcd with " << chunkSize << " tiles per thread" << std::endl;
        stop.reset(new StopWatch("rcd demosaic"));
    }

    double progress = 0.0;

    if (plistener) {
        plistener->setProgressStr(Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), M("TP_RAW_RCD")));
        plistener->setProgress(progress);
    }
    
    const unsigned int cfarray[2][2] = {{FC(0,0), FC(0,1)}, {FC(1,0), FC(1,1)}};
    constexpr int tileBorder = 9; // avoid tile-overlap errors
    constexpr int rcdBorder = 9;
    constexpr int tileSize = 194;
    constexpr int tileSizeN = tileSize - 2 * tileBorder;
    const int numTh = H / (tileSizeN) + ((H % (tileSizeN)) ? 1 : 0);
    const int numTw = W / (tileSizeN) + ((W % (tileSizeN)) ? 1 : 0);
    constexpr int w1 = tileSize, w2 = 2 * tileSize, w3 = 3 * tileSize, w4 = 4 * tileSize;
    //Tolerance to avoid dividing by zero
    constexpr float eps = 1e-5f;
    constexpr float epssq = 1e-10f;
    constexpr float scale = 65536.f;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    int progresscounter = 0;
    float *const cfa = (float*) calloc(tileSize * tileSize, sizeof *cfa);
    float (*const rgb)[tileSize * tileSize] = (float (*)[tileSize * tileSize])malloc(3 * sizeof *rgb);
    float *const VH_Dir = (float*) calloc(tileSize * tileSize, sizeof *VH_Dir);
    float *const PQ_Dir = (float*) calloc(tileSize * tileSize / 2, sizeof *PQ_Dir);
    float *const lpf = PQ_Dir; // reuse buffer, they don't overlap in usage
    float *const P_CDiff_Hpf = (float*) calloc(tileSize * tileSize / 2, sizeof *P_CDiff_Hpf);
    float *const Q_CDiff_Hpf = (float*) calloc(tileSize * tileSize / 2, sizeof *Q_CDiff_Hpf);

#ifdef _OPENMP
    #pragma omp for schedule(dynamic, chunkSize) collapse(2) nowait
#endif
    for (int tr = 0; tr < numTh; ++tr) {
        for (int tc = 0; tc < numTw; ++tc) {
            const int rowStart = tr * tileSizeN;
            const int rowEnd = std::min(rowStart + tileSize, H);
            if (rowStart + tileBorder == rowEnd - tileBorder) {
                continue;
            }
            const int colStart = tc * tileSizeN;
            const int colEnd = std::min(colStart + tileSize, W);
            if (colStart + tileBorder == colEnd - tileBorder) {
                continue;
            }

            const int tileRows = std::min(rowEnd - rowStart, tileSize);
            const int tilecols = std::min(colEnd - colStart, tileSize);

            for (int row = rowStart; row < rowEnd; row++) {
                const int c0 = fc(cfarray, row, colStart);
                const int c1 = fc(cfarray, row, colStart + 1);
                for (int col = colStart, indx = (row - rowStart) * tileSize; col < colEnd; ++col, ++indx) {
                    cfa[indx] = rgb[c0][indx] = rgb[c1][indx] = LIM01(rawData[row][col] / scale);
                }
            }

            // Step 1: Find cardinal and diagonal interpolation directions
            float bufferV[3][tileSize - 8];

            // Step 1.1: Calculate the square of the vertical and horizontal color difference high pass filter
            for (int row = 3; row < std::min(tileRows - 3, 5); ++row) {
                for (int col = 4, indx = row * tileSize + col; col < tilecols - 4; ++col, ++indx) {
                    bufferV[row - 3][col - 4] = SQR((cfa[indx - w3] - cfa[indx - w1] - cfa[indx + w1] + cfa[indx + w3]) - 3.f * (cfa[indx - w2] + cfa[indx + w2])  + 6.f * cfa[indx]);
                }
            }

            // Step 1.2: Obtain the vertical and horizontal directional discrimination strength
            float bufferH[tileSize - 6] ALIGNED16;
            float* V0 = bufferV[0];
            float* V1 = bufferV[1];
            float* V2 = bufferV[2];
            for (int row = 4; row < tileRows - 4; ++row) {
                for (int col = 3, indx = row * tileSize + col; col < tilecols - 3; ++col, ++indx) {
                    bufferH[col - 3] = SQR((cfa[indx -  3] - cfa[indx -  1] - cfa[indx +  1] + cfa[indx +  3]) - 3.f * (cfa[indx -  2] + cfa[indx +  2]) + 6.f * cfa[indx]);
                }
                for (int col = 4, indx = (row + 1) * tileSize + col; col < tilecols - 4; ++col, ++indx) {
                    V2[col - 4] = SQR((cfa[indx - w3] - cfa[indx - w1] - cfa[indx + w1] + cfa[indx + w3]) - 3.f * (cfa[indx - w2] + cfa[indx + w2])  + 6.f * cfa[indx]);
                }
                for (int col = 4, indx = row * tileSize + col; col < tilecols - 4; ++col, ++indx) {

                    float V_Stat = std::max(epssq, V0[col - 4] + V1[col - 4] + V2[col - 4]);
                    float H_Stat = std::max(epssq, bufferH[col -  4] + bufferH[col - 3] + bufferH[col -  2]);

                    VH_Dir[indx] = V_Stat / (V_Stat + H_Stat);
                }
                // rotate pointers from row0, row1, row2 to row1, row2, row0
                std::swap(V0, V2);
                std::swap(V0, V1);
            }

            // Step 2: Low pass filter incorporating green, red and blue local samples from the raw data
            for (int row = 2; row < tileRows - 2; ++row) {
                for (int col = 2 + (fc(cfarray, row, 0) & 1), indx = row * tileSize + col, lpindx = indx / 2; col < tilecols - 2; col += 2, indx += 2, ++lpindx) {
                    lpf[lpindx] = cfa[indx] +
                                  0.5f * (cfa[indx - w1] + cfa[indx + w1] + cfa[indx - 1] + cfa[indx + 1]) +
                                  0.25f * (cfa[indx - w1 - 1] + cfa[indx - w1 + 1] + cfa[indx + w1 - 1] + cfa[indx + w1 + 1]);
                }
            }

            // Step 3: Populate the green channel at blue and red CFA positions
            for (int row = 4; row < tileRows - 4; ++row) {
                for (int col = 4 + (fc(cfarray, row, 0) & 1), indx = row * tileSize + col, lpindx = indx / 2; col < tilecols - 4; col += 2, indx += 2, ++lpindx) {
                    // Cardinal gradients
                    const float cfai = cfa[indx];
                    const float N_Grad = eps + (std::fabs(cfa[indx - w1] - cfa[indx + w1]) + std::fabs(cfai - cfa[indx - w2])) + (std::fabs(cfa[indx - w1] - cfa[indx - w3]) + std::fabs(cfa[indx - w2] - cfa[indx - w4]));
                    const float S_Grad = eps + (std::fabs(cfa[indx - w1] - cfa[indx + w1]) + std::fabs(cfai - cfa[indx + w2])) + (std::fabs(cfa[indx + w1] - cfa[indx + w3]) + std::fabs(cfa[indx + w2] - cfa[indx + w4]));
                    const float W_Grad = eps + (std::fabs(cfa[indx -  1] - cfa[indx +  1]) + std::fabs(cfai - cfa[indx -  2])) + (std::fabs(cfa[indx -  1] - cfa[indx -  3]) + std::fabs(cfa[indx -  2] - cfa[indx -  4]));
                    const float E_Grad = eps + (std::fabs(cfa[indx -  1] - cfa[indx +  1]) + std::fabs(cfai - cfa[indx +  2])) + (std::fabs(cfa[indx +  1] - cfa[indx +  3]) + std::fabs(cfa[indx +  2] - cfa[indx +  4]));

                    // Cardinal pixel estimations
                    const float lpfi = lpf[lpindx];
                    const float N_Est = cfa[indx - w1] * (lpfi + lpfi) / (eps + lpfi + lpf[lpindx - w1]);
                    const float S_Est = cfa[indx + w1] * (lpfi + lpfi) / (eps + lpfi + lpf[lpindx + w1]);
                    const float W_Est = cfa[indx -  1] * (lpfi + lpfi) / (eps + lpfi + lpf[lpindx -  1]);
                    const float E_Est = cfa[indx +  1] * (lpfi + lpfi) / (eps + lpfi + lpf[lpindx +  1]);

                    // Vertical and horizontal estimations
                    const float V_Est = (S_Grad * N_Est + N_Grad * S_Est) / (N_Grad + S_Grad);
                    const float H_Est = (W_Grad * E_Est + E_Grad * W_Est) / (E_Grad + W_Grad);

                    // G@B and G@R interpolation
                    // Refined vertical and horizontal local discrimination
                    const float VH_Central_Value = VH_Dir[indx];
                    const float VH_Neighbourhood_Value = 0.25f * ((VH_Dir[indx - w1 - 1] + VH_Dir[indx - w1 + 1]) + (VH_Dir[indx + w1 - 1] + VH_Dir[indx + w1 + 1]));

                    const float VH_Disc = std::fabs(0.5f - VH_Central_Value) < std::fabs(0.5f - VH_Neighbourhood_Value) ? VH_Neighbourhood_Value : VH_Central_Value;
                    rgb[1][indx] = intp(VH_Disc, H_Est, V_Est);
                }
            }

            /**
            * STEP 4: Populate the red and blue channels
            */

            // Step 4.0: Calculate the square of the P/Q diagonals color difference high pass filter
            for (int row = 3; row < tileRows - 3; ++row) {
                for (int col = 3, indx = row * tileSize + col, indx2 = indx / 2; col < tilecols - 3; col+=2, indx+=2, indx2++ ) {
                    P_CDiff_Hpf[indx2] = SQR((cfa[indx - w3 - 3] - cfa[indx - w1 - 1] - cfa[indx + w1 + 1] + cfa[indx + w3 + 3]) - 3.f * (cfa[indx - w2 - 2] + cfa[indx + w2 + 2]) + 6.f * cfa[indx]);
                    Q_CDiff_Hpf[indx2] = SQR((cfa[indx - w3 + 3] - cfa[indx - w1 + 1] - cfa[indx + w1 - 1] + cfa[indx + w3 - 3]) - 3.f * (cfa[indx - w2 + 2] + cfa[indx + w2 - 2]) + 6.f * cfa[indx]);
                }
            }

            // Step 4.1: Obtain the P/Q diagonals directional discrimination strength
            for (int row = 4; row < tileRows - 4; ++row) {
                for (int col = 4 + (fc(cfarray, row, 0) & 1), indx = row * tileSize + col, indx2 = indx / 2, indx3 = (indx - w1 - 1) / 2, indx4 = (indx + w1 - 1) / 2; col < tilecols - 4; col += 2, indx += 2, indx2++, indx3++, indx4++ ) {
                    float P_Stat = std::max(epssq, P_CDiff_Hpf[indx3] + P_CDiff_Hpf[indx2] + P_CDiff_Hpf[indx4 + 1]);
                    float Q_Stat = std::max(epssq, Q_CDiff_Hpf[indx3 + 1] + Q_CDiff_Hpf[indx2] + Q_CDiff_Hpf[indx4]);
                    PQ_Dir[indx2] = P_Stat / (P_Stat + Q_Stat);
                }
            }

            // Step 4.2: Populate the red and blue channels at blue and red CFA positions
            for (int row = 4; row < tileRows - 4; ++row) {
                for (int col = 4 + (fc(cfarray, row, 0) & 1), indx = row * tileSize + col, c = 2 - fc(cfarray, row, col), pqindx = indx / 2, pqindx2 = (indx - w1 - 1) / 2, pqindx3 = (indx + w1 - 1) / 2; col < tilecols - 4; col += 2, indx += 2, ++pqindx, ++pqindx2, ++pqindx3) {

                    // Refined P/Q diagonal local discrimination
                    float PQ_Central_Value   = PQ_Dir[pqindx];
                    float PQ_Neighbourhood_Value = 0.25f * (PQ_Dir[pqindx2] + PQ_Dir[pqindx2 + 1] + PQ_Dir[pqindx3] + PQ_Dir[pqindx3 + 1]);

                    float PQ_Disc = (std::fabs(0.5f - PQ_Central_Value) < std::fabs(0.5f - PQ_Neighbourhood_Value)) ? PQ_Neighbourhood_Value : PQ_Central_Value;

                    // Diagonal gradients
                    float NW_Grad = eps + std::fabs(rgb[c][indx - w1 - 1] - rgb[c][indx + w1 + 1]) + std::fabs(rgb[c][indx - w1 - 1] - rgb[c][indx - w3 - 3]) + std::fabs(rgb[1][indx] - rgb[1][indx - w2 - 2]);
                    float NE_Grad = eps + std::fabs(rgb[c][indx - w1 + 1] - rgb[c][indx + w1 - 1]) + std::fabs(rgb[c][indx - w1 + 1] - rgb[c][indx - w3 + 3]) + std::fabs(rgb[1][indx] - rgb[1][indx - w2 + 2]);
                    float SW_Grad = eps + std::fabs(rgb[c][indx - w1 + 1] - rgb[c][indx + w1 - 1]) + std::fabs(rgb[c][indx + w1 - 1] - rgb[c][indx + w3 - 3]) + std::fabs(rgb[1][indx] - rgb[1][indx + w2 - 2]);
                    float SE_Grad = eps + std::fabs(rgb[c][indx - w1 - 1] - rgb[c][indx + w1 + 1]) + std::fabs(rgb[c][indx + w1 + 1] - rgb[c][indx + w3 + 3]) + std::fabs(rgb[1][indx] - rgb[1][indx + w2 + 2]);

                    // Diagonal colour differences
                    float NW_Est = rgb[c][indx - w1 - 1] - rgb[1][indx - w1 - 1];
                    float NE_Est = rgb[c][indx - w1 + 1] - rgb[1][indx - w1 + 1];
                    float SW_Est = rgb[c][indx + w1 - 1] - rgb[1][indx + w1 - 1];
                    float SE_Est = rgb[c][indx + w1 + 1] - rgb[1][indx + w1 + 1];

                    // P/Q estimations
                    float P_Est = (NW_Grad * SE_Est + SE_Grad * NW_Est) / (NW_Grad + SE_Grad);
                    float Q_Est = (NE_Grad * SW_Est + SW_Grad * NE_Est) / (NE_Grad + SW_Grad);

                    // R@B and B@R interpolation
                    rgb[c][indx] = rgb[1][indx] + intp(PQ_Disc, Q_Est, P_Est);
                }
            }

            // Step 4.3: Populate the red and blue channels at green CFA positions
            for (int row = 4; row < tileRows - 4; ++row) {
                for (int col = 4 + (fc(cfarray, row, 1) & 1), indx = row * tileSize + col; col < tilecols - 4; col += 2, indx += 2) {

                    // Refined vertical and horizontal local discrimination
                    float VH_Central_Value = VH_Dir[indx];
                    float VH_Neighbourhood_Value = 0.25f * ((VH_Dir[indx - w1 - 1] + VH_Dir[indx - w1 + 1]) + (VH_Dir[indx + w1 - 1] + VH_Dir[indx + w1 + 1]));

                    float VH_Disc = (std::fabs(0.5f - VH_Central_Value) < std::fabs(0.5f - VH_Neighbourhood_Value)) ? VH_Neighbourhood_Value : VH_Central_Value;
                    float rgb1 = rgb[1][indx];
                    float N1 = eps + std::fabs(rgb1 - rgb[1][indx - w2]);
                    float S1 = eps + std::fabs(rgb1 - rgb[1][indx + w2]);
                    float W1 = eps + std::fabs(rgb1 - rgb[1][indx -  2]);
                    float E1 = eps + std::fabs(rgb1 - rgb[1][indx +  2]);

                    float rgb1mw1 = rgb[1][indx - w1];
                    float rgb1pw1 = rgb[1][indx + w1];
                    float rgb1m1 = rgb[1][indx - 1];
                    float rgb1p1 = rgb[1][indx + 1];
                    for (int c = 0; c <= 2; c += 2) {
                        // Cardinal gradients
                        float SNabs = std::fabs(rgb[c][indx - w1] - rgb[c][indx + w1]);
                        float EWabs = std::fabs(rgb[c][indx -  1] - rgb[c][indx +  1]);
                        float N_Grad = N1 + SNabs + std::fabs(rgb[c][indx - w1] - rgb[c][indx - w3]);
                        float S_Grad = S1 + SNabs + std::fabs(rgb[c][indx + w1] - rgb[c][indx + w3]);
                        float W_Grad = W1 + EWabs + std::fabs(rgb[c][indx -  1] - rgb[c][indx -  3]);
                        float E_Grad = E1 + EWabs + std::fabs(rgb[c][indx +  1] - rgb[c][indx +  3]);

                        // Cardinal colour differences
                        float N_Est = rgb[c][indx - w1] - rgb1mw1;
                        float S_Est = rgb[c][indx + w1] - rgb1pw1;
                        float W_Est = rgb[c][indx -  1] - rgb1m1;
                        float E_Est = rgb[c][indx +  1] - rgb1p1;

                        // Vertical and horizontal estimations
                        float V_Est = (N_Grad * S_Est + S_Grad * N_Est) / (N_Grad + S_Grad);
                        float H_Est = (E_Grad * W_Est + W_Grad * E_Est) / (E_Grad + W_Grad);

                        // R@G and B@G interpolation
                        rgb[c][indx] = rgb1 + intp(VH_Disc, H_Est, V_Est);
                    }
                }
            }

            // For the outermost tiles in all directions we can use a smaller border margin
            const int firstVertical = rowStart + ((tr == 0) ? rcdBorder : tileBorder);
            const int lastVertical = rowEnd - ((tr == numTh - 1) ? rcdBorder : tileBorder);
            const int firstHorizontal = colStart + ((tc == 0) ? rcdBorder : tileBorder);
            const int lastHorizontal =  colEnd - ((tc == numTw - 1) ? rcdBorder : tileBorder);
            for (int row = firstVertical; row < lastVertical; ++row) {
                for (int col = firstHorizontal; col < lastHorizontal; ++col) {
                    int idx = (row - rowStart) * tileSize + col - colStart ;
                    red[row][col] = std::max(0.f, rgb[0][idx] * scale);
                    green[row][col] = std::max(0.f, rgb[1][idx] * scale);
                    blue[row][col] = std::max(0.f, rgb[2][idx] * scale);
                }
            }

            if (plistener) {
                progresscounter++;
                if (progresscounter % 32 == 0) {
#ifdef _OPENMP
                    #pragma omp critical (rcdprogress)
#endif
                    {
                        progress += (double)32 * ((tileSizeN) * (tileSizeN)) / (H * W);
                        progress = progress > 1.0 ? 1.0 : progress;
                        plistener->setProgress(progress);
                    }
                }
            }
        }
    }

    free(cfa);
    free(rgb);
    free(VH_Dir);
    free(PQ_Dir);
    free(P_CDiff_Hpf);
    free(Q_CDiff_Hpf);
}

    border_interpolate(W, H, rcdBorder, rawData, red, green, blue);

    if (plistener) {
        plistener->setProgress(1);
    }
}

} /* namespace */
