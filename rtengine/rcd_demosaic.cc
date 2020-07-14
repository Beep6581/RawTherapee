/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017-2018 Luis Sanz Rodriguez (luis.sanz.rodriguez(at)gmail(dot)com) and Ingo Weyrich (heckflosse67@gmx.de)
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
#include "../rtgui/multilangmgr.h"
#include "opthelper.h"
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

/**
* RATIO CORRECTED DEMOSAICING
* Luis Sanz Rodriguez (luis.sanz.rodriguez(at)gmail(dot)com)
*
* Release 2.3 @ 171125
*
* Original code from https://github.com/LuisSR/RCD-Demosaicing
* Licensed under the GNU GPL version 3
*/
// Tiled version by Ingo Weyrich (heckflosse67@gmx.de)
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
    constexpr int rcdBorder = 9;
    constexpr int tileSize = 214;
    constexpr int tileSizeN = tileSize - 2 * rcdBorder;
    const int numTh = H / (tileSizeN) + ((H % (tileSizeN)) ? 1 : 0);
    const int numTw = W / (tileSizeN) + ((W % (tileSizeN)) ? 1 : 0);
    constexpr int w1 = tileSize, w2 = 2 * tileSize, w3 = 3 * tileSize, w4 = 4 * tileSize;
    //Tolerance to avoid dividing by zero
    static constexpr float eps = 1e-5f;
    static constexpr float epssq = 1e-10f;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    int progresscounter = 0;
    float *cfa = (float*) calloc(tileSize * tileSize, sizeof *cfa);
    float (*rgb)[tileSize * tileSize] = (float (*)[tileSize * tileSize])malloc(3 * sizeof *rgb);
    float *VH_Dir = (float*) calloc(tileSize * tileSize, sizeof *VH_Dir);
    float *PQ_Dir = (float*) calloc(tileSize * tileSize, sizeof *PQ_Dir);
    float *lpf = PQ_Dir; // reuse buffer, they don't overlap in usage

#ifdef _OPENMP
    #pragma omp for schedule(dynamic, chunkSize) collapse(2) nowait
#endif
    for(int tr = 0; tr < numTh; ++tr) {
        for(int tc = 0; tc < numTw; ++tc) {
            const int rowStart = tr * tileSizeN;
            const int rowEnd = std::min(rowStart + tileSize, H);
            if(rowStart + rcdBorder == rowEnd - rcdBorder) {
                continue;
            }
            const int colStart = tc * tileSizeN;
            const int colEnd = std::min(colStart + tileSize, W);
            if(colStart + rcdBorder == colEnd - rcdBorder) {
                continue;
            }

            const int tileRows = std::min(rowEnd - rowStart, tileSize);
            const int tilecols = std::min(colEnd - colStart, tileSize);

            for (int row = rowStart; row < rowEnd; row++) {
                int indx = (row - rowStart) * tileSize;
                int c0 = fc(cfarray, row, colStart);
                int c1 = fc(cfarray, row, colStart + 1);
                int col = colStart;

                for (; col < colEnd - 1; col+=2, indx+=2) {
                    cfa[indx] = rgb[c0][indx] = LIM01(rawData[row][col] / 65535.f);
                    cfa[indx + 1] = rgb[c1][indx + 1] = LIM01(rawData[row][col + 1] / 65535.f);
                }
                if(col < colEnd) {
                    cfa[indx] = rgb[c0][indx] = LIM01(rawData[row][col] / 65535.f);
                }
            }

            /**
            * STEP 1: Find cardinal and diagonal interpolation directions
            */

            for (int row = 4; row < tileRows - 4; row++) {
                for (int col = 4, indx = row * tileSize + col; col < tilecols - 4; col++, indx++) {
                    const float cfai = cfa[indx];
                    //Calculate h/v local discrimination
                    float V_Stat = max(epssq, -18.f * cfai * (cfa[indx - w1] + cfa[indx + w1] + 2.f * (cfa[indx - w2] + cfa[indx + w2]) - cfa[indx - w3] - cfa[indx + w3]) - 2.f * cfai * (cfa[indx - w4] + cfa[indx + w4] - 19.f * cfai) - cfa[indx - w1] * (70.f * cfa[indx + w1] + 12.f * cfa[indx - w2] - 24.f * cfa[indx + w2] + 38.f * cfa[indx - w3] - 16.f * cfa[indx + w3] - 12.f * cfa[indx - w4] + 6.f * cfa[indx + w4] - 46.f * cfa[indx - w1]) + cfa[indx + w1] * (24.f * cfa[indx - w2] - 12.f * cfa[indx + w2] + 16.f * cfa[indx - w3] - 38.f * cfa[indx + w3] - 6.f * cfa[indx - w4] + 12.f * cfa[indx + w4] + 46.f * cfa[indx + w1]) + cfa[indx - w2] * (14.f * cfa[indx + w2] - 12.f * cfa[indx + w3] - 2.f * cfa[indx - w4] + 2.f * cfa[indx + w4] + 11.f * cfa[indx - w2]) + cfa[indx + w2] * (-12.f * cfa[indx - w3] + 2.f * (cfa[indx - w4] - cfa[indx + w4]) + 11.f * cfa[indx + w2]) + cfa[indx - w3] * (2.f * cfa[indx + w3] - 6.f * cfa[indx - w4] + 10.f * cfa[indx - w3]) + cfa[indx + w3] * (-6.f * cfa[indx + w4] + 10.f * cfa[indx + w3]) + cfa[indx - w4] * cfa[indx - w4] + cfa[indx + w4] * cfa[indx + w4]);
                    float H_Stat = max(epssq, -18.f * cfai * (cfa[indx -  1] + cfa[indx +  1] + 2.f * (cfa[indx -  2] + cfa[indx +  2]) - cfa[indx -  3] - cfa[indx +  3]) - 2.f * cfai * (cfa[indx -  4] + cfa[indx +  4] - 19.f * cfai) - cfa[indx -  1] * (70.f * cfa[indx +  1] + 12.f * cfa[indx -  2] - 24.f * cfa[indx +  2] + 38.f * cfa[indx -  3] - 16.f * cfa[indx +  3] - 12.f * cfa[indx -  4] + 6.f * cfa[indx +  4] - 46.f * cfa[indx -  1]) + cfa[indx +  1] * (24.f * cfa[indx -  2] - 12.f * cfa[indx +  2] + 16.f * cfa[indx -  3] - 38.f * cfa[indx +  3] - 6.f * cfa[indx -  4] + 12.f * cfa[indx +  4] + 46.f * cfa[indx +  1]) + cfa[indx -  2] * (14.f * cfa[indx +  2] - 12.f * cfa[indx +  3] - 2.f * cfa[indx -  4] + 2.f * cfa[indx +  4] + 11.f * cfa[indx -  2]) + cfa[indx +  2] * (-12.f * cfa[indx -  3] + 2.f * (cfa[indx -  4] - cfa[indx +  4]) + 11.f * cfa[indx +  2]) + cfa[indx -  3] * (2.f * cfa[indx +  3] - 6.f * cfa[indx -  4] + 10.f * cfa[indx -  3]) + cfa[indx +  3] * (-6.f * cfa[indx +  4] + 10.f * cfa[indx +  3]) + cfa[indx -  4] * cfa[indx -  4] + cfa[indx +  4] * cfa[indx +  4]);

                    VH_Dir[indx] = V_Stat / (V_Stat + H_Stat);
                }
            }

            /**
            * STEP 2: Calculate the low pass filter
            */
            // Step 2.1: Low pass filter incorporating green, red and blue local samples from the raw data

            for (int row = 2; row < tileRows - 2; row++) {
                for (int col = 2 + (fc(cfarray, row, 0) & 1), indx = row * tileSize + col; col < tilecols - 2; col += 2, indx += 2) {
                    lpf[indx>>1] = 0.25f * cfa[indx] + 0.125f * (cfa[indx - w1] + cfa[indx + w1] + cfa[indx - 1] + cfa[indx + 1]) + 0.0625f * (cfa[indx - w1 - 1] + cfa[indx - w1 + 1] + cfa[indx + w1 - 1] + cfa[indx + w1 + 1]);
                }
            }

            /**
            * STEP 3: Populate the green channel
            */
            // Step 3.1: Populate the green channel at blue and red CFA positions
           for (int row = 4; row < tileRows - 4; row++) {
               for (int col = 4 + (fc(cfarray, row, 0) & 1), indx = row * tileSize + col; col < tilecols - 4; col += 2, indx += 2) {

                    // Refined vertical and horizontal local discrimination
                    float VH_Central_Value = VH_Dir[indx];
                    float VH_Neighbourhood_Value = 0.25f * ((VH_Dir[indx - w1 - 1] + VH_Dir[indx - w1 + 1]) + (VH_Dir[indx + w1 - 1] + VH_Dir[indx + w1 + 1]));

                    float VH_Disc = std::fabs(0.5f - VH_Central_Value) < std::fabs(0.5f - VH_Neighbourhood_Value) ? VH_Neighbourhood_Value : VH_Central_Value;

                    // Cardinal gradients
                    float N_Grad = eps + std::fabs(cfa[indx - w1] - cfa[indx + w1]) + std::fabs(cfa[indx] - cfa[indx - w2]) + std::fabs(cfa[indx - w1] - cfa[indx - w3]) + std::fabs(cfa[indx - w2] - cfa[indx - w4]);
                    float S_Grad = eps + std::fabs(cfa[indx - w1] - cfa[indx + w1]) + std::fabs(cfa[indx] - cfa[indx + w2]) + std::fabs(cfa[indx + w1] - cfa[indx + w3]) + std::fabs(cfa[indx + w2] - cfa[indx + w4]);
                    float W_Grad = eps + std::fabs(cfa[indx -  1] - cfa[indx +  1]) + std::fabs(cfa[indx] - cfa[indx -  2]) + std::fabs(cfa[indx -  1] - cfa[indx -  3]) + std::fabs(cfa[indx -  2] - cfa[indx -  4]);
                    float E_Grad = eps + std::fabs(cfa[indx -  1] - cfa[indx +  1]) + std::fabs(cfa[indx] - cfa[indx +  2]) + std::fabs(cfa[indx +  1] - cfa[indx +  3]) + std::fabs(cfa[indx +  2] - cfa[indx +  4]);

                    // Cardinal pixel estimations
                    float N_Est = cfa[indx - w1] * (1.f + (lpf[indx>>1] - lpf[(indx - w2)>>1]) / (eps + lpf[indx>>1] + lpf[(indx - w2)>>1]));
                    float S_Est = cfa[indx + w1] * (1.f + (lpf[indx>>1] - lpf[(indx + w2)>>1]) / (eps + lpf[indx>>1] + lpf[(indx + w2)>>1]));
                    float W_Est = cfa[indx -  1] * (1.f + (lpf[indx>>1] - lpf[(indx -  2)>>1]) / (eps + lpf[indx>>1] + lpf[(indx -  2)>>1]));
                    float E_Est = cfa[indx +  1] * (1.f + (lpf[indx>>1] - lpf[(indx +  2)>>1]) / (eps + lpf[indx>>1] + lpf[(indx +  2)>>1]));

                    // Vertical and horizontal estimations
                    float V_Est = (S_Grad * N_Est + N_Grad * S_Est) / (N_Grad + S_Grad);
                    float H_Est = (W_Grad * E_Est + E_Grad * W_Est) / (E_Grad + W_Grad);

                    // G@B and G@R interpolation
                    rgb[1][indx] = VH_Disc * H_Est + (1.f - VH_Disc) * V_Est;

                }
            }
            /**
            * STEP 4: Populate the red and blue channels
            */

            // Step 4.1: Calculate P/Q diagonal local discrimination
            for (int row = rcdBorder - 4; row < tileRows - rcdBorder + 4; row++) {
                for (int col = rcdBorder - 4 + (fc(cfarray, row, rcdBorder) & 1), indx = row * tileSize + col; col < tilecols - rcdBorder + 4; col += 2, indx += 2) {
                    const float cfai = cfa[indx];

                    float P_Stat = max(epssq, - 18.f * cfai * (cfa[indx - w1 - 1] + cfa[indx + w1 + 1] + 2.f * (cfa[indx - w2 - 2] + cfa[indx + w2 + 2]) - cfa[indx - w3 - 3] - cfa[indx + w3 + 3]) - 2.f * cfai * (cfa[indx - w4 - 4] + cfa[indx + w4 + 4] - 19.f * cfai) - cfa[indx - w1 - 1] * (70.f * cfa[indx + w1 + 1] - 12.f * cfa[indx - w2 - 2] + 24.f * cfa[indx + w2 + 2] - 38.f * cfa[indx - w3 - 3] + 16.f * cfa[indx + w3 + 3] + 12.f * cfa[indx - w4 - 4] - 6.f * cfa[indx + w4 + 4] + 46.f * cfa[indx - w1 - 1]) + cfa[indx + w1 + 1] * (24.f * cfa[indx - w2 - 2] - 12.f * cfa[indx + w2 + 2] + 16.f * cfa[indx - w3 - 3] - 38.f * cfa[indx + w3 + 3] - 6.f * cfa[indx - w4 - 4] + 12.f * cfa[indx + w4 + 4] + 46.f * cfa[indx + w1 + 1]) + cfa[indx - w2 - 2] * (14.f * cfa[indx + w2 + 2] - 12.f * cfa[indx + w3 + 3] - 2.f * (cfa[indx - w4 - 4] - cfa[indx + w4 + 4]) + 11.f * cfa[indx - w2 - 2]) - cfa[indx + w2 + 2] * (12.f * cfa[indx - w3 - 3] + 2.f * (cfa[indx - w4 - 4] - cfa[indx + w4 + 4]) + 11.f * cfa[indx + w2 + 2]) + cfa[indx - w3 - 3] * (2.f * cfa[indx + w3 + 3] - 6.f * cfa[indx - w4 - 4] + 10.f * cfa[indx - w3 - 3]) - cfa[indx + w3 + 3] * (6.f * cfa[indx + w4 + 4] + 10.f * cfa[indx + w3 + 3]) + cfa[indx - w4 - 4] * cfa[indx - w4 - 4] + cfa[indx + w4 + 4] * cfa[indx + w4 + 4]);
                    float Q_Stat = max(epssq, - 18.f * cfai * (cfa[indx + w1 - 1] + cfa[indx - w1 + 1] + 2.f * (cfa[indx + w2 - 2] + cfa[indx - w2 + 2]) - cfa[indx + w3 - 3] - cfa[indx - w3 + 3]) - 2.f * cfai * (cfa[indx + w4 - 4] + cfa[indx - w4 + 4] - 19.f * cfai) - cfa[indx + w1 - 1] * (70.f * cfa[indx - w1 + 1] - 12.f * cfa[indx + w2 - 2] + 24.f * cfa[indx - w2 + 2] - 38.f * cfa[indx + w3 - 3] + 16.f * cfa[indx - w3 + 3] + 12.f * cfa[indx + w4 - 4] - 6.f * cfa[indx - w4 + 4] + 46.f * cfa[indx + w1 - 1]) + cfa[indx - w1 + 1] * (24.f * cfa[indx + w2 - 2] - 12.f * cfa[indx - w2 + 2] + 16.f * cfa[indx + w3 - 3] - 38.f * cfa[indx - w3 + 3] - 6.f * cfa[indx + w4 - 4] + 12.f * cfa[indx - w4 + 4] + 46.f * cfa[indx - w1 + 1]) + cfa[indx + w2 - 2] * (14.f * cfa[indx - w2 + 2] - 12.f * cfa[indx - w3 + 3] - 2.f * (cfa[indx + w4 - 4] - cfa[indx - w4 + 4]) + 11.f * cfa[indx + w2 - 2]) - cfa[indx - w2 + 2] * (12.f * cfa[indx + w3 - 3] + 2.f * (cfa[indx + w4 - 4] - cfa[indx - w4 + 4]) + 11.f * cfa[indx - w2 + 2]) + cfa[indx + w3 - 3] * (2.f * cfa[indx - w3 + 3] - 6.f * cfa[indx + w4 - 4] + 10.f * cfa[indx + w3 - 3]) - cfa[indx - w3 + 3] * (6.f * cfa[indx - w4 + 4] + 10.f * cfa[indx - w3 + 3]) + cfa[indx + w4 - 4] * cfa[indx + w4 - 4] + cfa[indx - w4 + 4] * cfa[indx - w4 + 4]);

                    PQ_Dir[indx] = P_Stat / (P_Stat + Q_Stat);

                }
            }

            // Step 4.2: Populate the red and blue channels at blue and red CFA positions
            for (int row = rcdBorder - 3; row < tileRows - rcdBorder + 3; row++) {
                for (int col = rcdBorder - 3 + (fc(cfarray, row, rcdBorder - 1) & 1), indx = row * tileSize + col, c = 2 - fc(cfarray, row, col); col < tilecols - rcdBorder + 3; col += 2, indx += 2) {

                    // Refined P/Q diagonal local discrimination
                    float PQ_Central_Value   = PQ_Dir[indx];
                    float PQ_Neighbourhood_Value = 0.25f * (PQ_Dir[indx - w1 - 1] + PQ_Dir[indx - w1 + 1] + PQ_Dir[indx + w1 - 1] + PQ_Dir[indx + w1 + 1]);

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
                    rgb[c][indx] = rgb[1][indx] + (1.f - PQ_Disc) * P_Est + PQ_Disc * Q_Est;

                }
            }

            // Step 4.3: Populate the red and blue channels at green CFA positions
            for (int row = rcdBorder; row < tileRows - rcdBorder; row++) {
                for (int col = rcdBorder + (fc(cfarray, row, rcdBorder - 1) & 1), indx = row * tileSize + col; col < tilecols - rcdBorder; col += 2, indx += 2) {

                    // Refined vertical and horizontal local discrimination
                    float VH_Central_Value   = VH_Dir[indx];
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
                        rgb[c][indx] = rgb1 + (1.f - VH_Disc) * V_Est + VH_Disc * H_Est;

                    }
                }
            }

            for (int row = rowStart + rcdBorder; row < rowEnd - rcdBorder; ++row) {
                for (int col = colStart + rcdBorder; col < colEnd - rcdBorder; ++col) {
                    int idx = (row - rowStart) * tileSize + col - colStart ;
                    red[row][col] = std::max(0.f, rgb[0][idx] * 65535.f);
                    green[row][col] = std::max(0.f, rgb[1][idx] * 65535.f);
                    blue[row][col] = std::max(0.f, rgb[2][idx] * 65535.f);
                }
            }

            if(plistener) {
                progresscounter++;
                if(progresscounter % 32 == 0) {
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
}

    border_interpolate(W, H, rcdBorder, rawData, red, green, blue);

    if (plistener) {
        plistener->setProgress(1);
    }
    // -------------------------------------------------------------------------
}

} /* namespace */
