////////////////////////////////////////////////////////////////
//
//          CFA denoise by wavelet transform, FT filtering
//
//  copyright (c) 2008-2012  Emil Martinec <ejmartin@uchicago.edu>
//
//
//  code dated: March 9, 2012
//
//  FTblockDN.cc is free software: you can redistribute it and/or modify
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
#include <fftw3.h>
#include "../rtgui/threadutils.h"
#include "rtengine.h"
#include "improcfun.h"
#include "LUT.h"
#include "array2D.h"
#include "iccmatrices.h"
#include "boxblur.h"
#include "rt_math.h"
#include "mytime.h"
#include "sleef.c"
#include "opthelper.h"
#include "cplx_wavelet_dec.h"
#include "median.h"
#include "iccstore.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define TS 64       // Tile size
#define offset 25   // shift between tiles
#define fTS ((TS/2+1))  // second dimension of Fourier tiles
#define blkrad 1    // radius of block averaging

#define epsilon 0.001f/(TS*TS) //tolerance

namespace rtengine
{

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
 Structure of the algorithm:

 1. Compute an initial denoise of the image via undecimated wavelet transform
 and universal thresholding modulated by user input.
 2. Decompose the residual image into TSxTS size tiles, shifting by 'offset' each step
 (so roughly each pixel is in (TS/offset)^2 tiles); Discrete Cosine transform the tiles.
 3. Filter the DCT data to pick out patterns missed by the wavelet denoise
 4. Inverse DCT the denoised tile data and combine the tiles into a denoised output image.

 */

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


extern const Settings* settings;

void ImProcFunctions::Median_Denoise(float **src, float **dst, const int width, const int height, const Median medianType, const int iterations, const int numThreads, float **buffer)
{
    int border = 1;

    switch (medianType) {
        case Median::TYPE_3X3_SOFT:
        case Median::TYPE_3X3_STRONG: {
            border = 1;
            break;
        }

        case Median::TYPE_5X5_SOFT: {
            border = 2;
            break;
        }

        case Median::TYPE_5X5_STRONG: {
            border = 2;
            break;
        }

        case Median::TYPE_7X7: {
            border = 3;
            break;
        }

        case Median::TYPE_9X9: {
            border = 4;
            break;
        }
    }

    float **allocBuffer = nullptr;
    float **medBuffer[2];
    medBuffer[0] = src;

    // we need a buffer if src == dst or if (src != dst && iterations > 1)
    if (src == dst || iterations > 1) {
        if (buffer == nullptr) { // we didn't get a buufer => create one
            allocBuffer = new float*[height];

            for (int i = 0; i < height; ++i) {
                allocBuffer[i] = new float[width];
            }

            medBuffer[1] = allocBuffer;
        } else { // we got a buffer => use it
            medBuffer[1] = buffer;
        }
    } else { // we can write directly into destination
        medBuffer[1] = dst;
    }

    float ** medianIn, ** medianOut;
    int BufferIndex = 0;

    for (int iteration = 1; iteration <= iterations; ++iteration) {
        medianIn = medBuffer[BufferIndex];
        medianOut = medBuffer[BufferIndex ^ 1];

        if (iteration == 1) { // upper border
            for (int i = 0; i < border; ++i) {
                for (int j = 0; j < width; ++j) {
                    medianOut[i][j] = medianIn[i][j];
                }
            }
        }

#ifdef _OPENMP
        #pragma omp parallel for num_threads(numThreads) if (numThreads>1) schedule(dynamic,16)
#endif

        for (int i = border; i < height - border; ++i) {
            int j = 0;

            for (; j < border; ++j) {
                medianOut[i][j] = medianIn[i][j];
            }

            switch (medianType) {
                case Median::TYPE_3X3_SOFT: {
                    for (; j < width - border; ++j) {
                        medianOut[i][j] = median(
                                              medianIn[i - 1][j],
                                              medianIn[i][j - 1],
                                              medianIn[i][j],
                                              medianIn[i][j + 1],
                                              medianIn[i + 1][j]
                                          );
                    }

                    break;
                }

                case Median::TYPE_3X3_STRONG: {
                    for (; j < width - border; ++j) {
                        medianOut[i][j] = median(
                                              medianIn[i - 1][j - 1],
                                              medianIn[i - 1][j],
                                              medianIn[i - 1][j + 1],
                                              medianIn[i][j - 1],
                                              medianIn[i][j],
                                              medianIn[i][j + 1],
                                              medianIn[i + 1][j - 1],
                                              medianIn[i + 1][j],
                                              medianIn[i + 1][j + 1]
                                          );
                    }

                    break;
                }

                case Median::TYPE_5X5_SOFT: {
                    for (; j < width - border; ++j) {
                        medianOut[i][j] = median(
                                              medianIn[i - 2][j],
                                              medianIn[i - 1][j - 1],
                                              medianIn[i - 1][j],
                                              medianIn[i - 1][j + 1],
                                              medianIn[i][j - 2],
                                              medianIn[i][j - 1],
                                              medianIn[i][j],
                                              medianIn[i][j + 1],
                                              medianIn[i][j + 2],
                                              medianIn[i + 1][j - 1],
                                              medianIn[i + 1][j],
                                              medianIn[i + 1][j + 1],
                                              medianIn[i + 2][j]
                                          );
                    }

                    break;
                }

                case Median::TYPE_5X5_STRONG: {
#ifdef __SSE2__

                    for (; j < width - border - 3; j += 4) {
                        STVFU(
                            medianOut[i][j],
                            median(
                                LVFU(medianIn[i - 2][j - 2]),
                                LVFU(medianIn[i - 2][j - 1]),
                                LVFU(medianIn[i - 2][j]),
                                LVFU(medianIn[i - 2][j + 1]),
                                LVFU(medianIn[i - 2][j + 2]),
                                LVFU(medianIn[i - 1][j - 2]),
                                LVFU(medianIn[i - 1][j - 1]),
                                LVFU(medianIn[i - 1][j]),
                                LVFU(medianIn[i - 1][j + 1]),
                                LVFU(medianIn[i - 1][j + 2]),
                                LVFU(medianIn[i][j - 2]),
                                LVFU(medianIn[i][j - 1]),
                                LVFU(medianIn[i][j]),
                                LVFU(medianIn[i][j + 1]),
                                LVFU(medianIn[i][j + 2]),
                                LVFU(medianIn[i + 1][j - 2]),
                                LVFU(medianIn[i + 1][j - 1]),
                                LVFU(medianIn[i + 1][j]),
                                LVFU(medianIn[i + 1][j + 1]),
                                LVFU(medianIn[i + 1][j + 2]),
                                LVFU(medianIn[i + 2][j - 2]),
                                LVFU(medianIn[i + 2][j - 1]),
                                LVFU(medianIn[i + 2][j]),
                                LVFU(medianIn[i + 2][j + 1]),
                                LVFU(medianIn[i + 2][j + 2])
                            )
                        );
                    }

#endif

                    for (; j < width - border; ++j) {
                        medianOut[i][j] = median(
                                              medianIn[i - 2][j - 2],
                                              medianIn[i - 2][j - 1],
                                              medianIn[i - 2][j],
                                              medianIn[i - 2][j + 1],
                                              medianIn[i - 2][j + 2],
                                              medianIn[i - 1][j - 2],
                                              medianIn[i - 1][j - 1],
                                              medianIn[i - 1][j],
                                              medianIn[i - 1][j + 1],
                                              medianIn[i - 1][j + 2],
                                              medianIn[i][j - 2],
                                              medianIn[i][j - 1],
                                              medianIn[i][j],
                                              medianIn[i][j + 1],
                                              medianIn[i][j + 2],
                                              medianIn[i + 1][j - 2],
                                              medianIn[i + 1][j - 1],
                                              medianIn[i + 1][j],
                                              medianIn[i + 1][j + 1],
                                              medianIn[i + 1][j + 2],
                                              medianIn[i + 2][j - 2],
                                              medianIn[i + 2][j - 1],
                                              medianIn[i + 2][j],
                                              medianIn[i + 2][j + 1],
                                              medianIn[i + 2][j + 2]
                                          );
                    }

                    break;
                }

                case Median::TYPE_7X7: {
#ifdef __SSE2__
                    std::array<vfloat, 49> vpp ALIGNED16;

                    for (; j < width - border - 3; j += 4) {
                        for (int kk = 0, ii = -border; ii <= border; ++ii) {
                            for (int jj = -border; jj <= border; ++jj, ++kk) {
                                vpp[kk] = LVFU(medianIn[i + ii][j + jj]);
                            }
                        }

                        STVFU(medianOut[i][j], median(vpp));
                    }

#endif

                    std::array<float, 49> pp;

                    for (; j < width - border; ++j) {
                        for (int kk = 0, ii = -border; ii <= border; ++ii) {
                            for (int jj = -border; jj <= border; ++jj, ++kk) {
                                pp[kk] = medianIn[i + ii][j + jj];
                            }
                        }

                        medianOut[i][j] = median(pp);
                    }

                    break;
                }

                case Median::TYPE_9X9: {
#ifdef __SSE2__
                    std::array<vfloat, 81> vpp ALIGNED16;

                    for (; j < width - border - 3; j += 4) {
                        for (int kk = 0, ii = -border; ii <= border; ++ii) {
                            for (int jj = -border; jj <= border; ++jj, ++kk) {
                                vpp[kk] = LVFU(medianIn[i + ii][j + jj]);
                            }
                        }

                        STVFU(medianOut[i][j], median(vpp));
                    }

#endif

                    std::array<float, 81> pp;

                    for (; j < width - border; ++j) {
                        for (int kk = 0, ii = -border; ii <= border; ++ii) {
                            for (int jj = -border; jj <= border; ++jj, ++kk) {
                                pp[kk] = medianIn[i + ii][j + jj];
                            }
                        }

                        medianOut[i][j] = median(pp);
                    }

                    for (; j < width; ++j) {
                        medianOut[i][j] = medianIn[i][j];
                    }

                    break;
                }
            }

            for (; j < width; ++j) {
                medianOut[i][j] = medianIn[i][j];
            }
        }

        if (iteration == 1) { // lower border
            for (int i = height - border; i < height; ++i) {
                for (int j = 0; j < width; ++j) {
                    medianOut[i][j] = medianIn[i][j];
                }
            }
        }

        BufferIndex ^= 1; // swap buffers
    }

    if (medianOut != dst) {
#ifdef _OPENMP
        #pragma omp parallel for num_threads(numThreads) if (numThreads>1)
#endif

        for (int i = border; i < height - border; ++i) {
            for (int j = border; j < width - border; ++j) {
                dst[i][j] = medianOut[i][j];
            }
        }
    }

    if (allocBuffer != nullptr) { // we allocated memory, so let's free it now
        for (int i = 0; i < height; ++i) {
            delete[] allocBuffer[i];
        }

        delete[] allocBuffer;
    }
}

void ImProcFunctions::Tile_calc(int tilesize, int overlap, int kall, int imwidth, int imheight, int &numtiles_W, int &numtiles_H, int &tilewidth, int &tileheight, int &tileWskip, int &tileHskip)

{
    if (kall == 2) {

        if (imwidth < tilesize) {
            numtiles_W = 1;
            tileWskip = imwidth;
            tilewidth = imwidth;
        } else {
            numtiles_W = ceil((static_cast<float>(imwidth)) / (tilesize - overlap));
            tilewidth  = ceil((static_cast<float>(imwidth)) / (numtiles_W)) + overlap;
            tilewidth += (tilewidth & 1);
            tileWskip = tilewidth - overlap;
        }

        if (imheight < tilesize) {
            numtiles_H = 1;
            tileHskip = imheight;
            tileheight = imheight;
        } else {
            numtiles_H = ceil((static_cast<float>(imheight)) / (tilesize - overlap));
            tileheight = ceil((static_cast<float>(imheight)) / (numtiles_H)) + overlap;
            tileheight += (tileheight & 1);
            tileHskip = tileheight - overlap;
        }
    }

    if (kall == 0) {
        numtiles_W = 1;
        tileWskip = imwidth;
        tilewidth = imwidth;
        numtiles_H = 1;
        tileHskip = imheight;
        tileheight = imheight;
    }

    //  printf("Nw=%d NH=%d tileW=%d tileH=%d\n",numtiles_W,numtiles_H,tileWskip,tileHskip);
}

int denoiseNestedLevels = 1;
enum nrquality {QUALITY_STANDARD, QUALITY_HIGH};

SSEFUNCTION void ImProcFunctions::RGB_denoise(int kall, Imagefloat * src, Imagefloat * dst, Imagefloat * calclum, float * ch_M, float *max_r, float *max_b, bool isRAW, const procparams::DirPyrDenoiseParams & dnparams, const double expcomp, const NoiseCurve & noiseLCurve, const NoiseCurve & noiseCCurve, float &chaut, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &nresi, float &highresi)
{
//#ifdef _DEBUG
    MyTime t1e, t2e;
    t1e.set();

//#endif
    if (dnparams.luma == 0 && dnparams.chroma == 0  && !dnparams.median && !noiseLCurve && !noiseCCurve) {
        //nothing to do; copy src to dst or do nothing in case src == dst
        if (src != dst) {
            src->copyData(dst);
        }

        if (calclum) {
            delete calclum;
            calclum = nullptr;
        }

        return;
    }

    static MyMutex FftwMutex;
    MyMutex::MyLock lock(FftwMutex);

    const nrquality nrQuality = (dnparams.smethod == "shal") ? QUALITY_STANDARD : QUALITY_HIGH;//shrink method
    const float qhighFactor = (nrQuality == QUALITY_HIGH) ? 1.f / static_cast<float>( settings->nrhigh) : 1.0f;
    const bool useNoiseCCurve = (noiseCCurve && noiseCCurve.getSum() > 5.f);
    const bool useNoiseLCurve = (noiseLCurve && noiseLCurve.getSum() >= 7.f);
    const bool autoch = (settings->leveldnautsimpl == 1 && (dnparams.Cmethod == "AUT" || dnparams.Cmethod == "PRE")) || (settings->leveldnautsimpl == 0 && (dnparams.C2method == "AUTO" || dnparams.C2method == "PREV"));

    float** lumcalc;
    float* lumcalcBuffer;
    float** ccalc;
    float* ccalcBuffer;

    bool ponder = false;
    float ponderCC = 1.f;

    if (settings->leveldnautsimpl == 1 && params->dirpyrDenoise.Cmethod == "PON") {
        ponder = true;
        ponderCC = 0.5f;
    }

    if (settings->leveldnautsimpl == 1 && params->dirpyrDenoise.Cmethod == "PRE") {
        ponderCC = 0.5f;
    }

    if (settings->leveldnautsimpl == 0 && params->dirpyrDenoise.Cmethod == "PREV") {
        ponderCC = 0.5f;
    }

    int metchoice = 0;

    if (dnparams.methodmed == "Lonly") {
        metchoice = 1;
    } else if (dnparams.methodmed == "Lab") {
        metchoice = 2;
    } else if (dnparams.methodmed == "ab") {
        metchoice = 3;
    } else if (dnparams.methodmed == "Lpab") {
        metchoice = 4;
    }

    const bool denoiseMethodRgb = (dnparams.dmethod == "RGB");
    // init luma noisevarL
    const float noiseluma = static_cast<float>( dnparams.luma);
    const float noisevarL = (useNoiseLCurve && (denoiseMethodRgb || !isRAW)) ? static_cast<float> (SQR(((noiseluma + 1.0) / 125.0) * (10. + (noiseluma + 1.0) / 25.0))) : static_cast<float> (SQR((noiseluma / 125.0) * (1.0 + noiseluma / 25.0)));
    const bool denoiseLuminance = (noisevarL > 0.00001f);

    //printf("NL=%f \n",noisevarL);
    if (useNoiseLCurve || useNoiseCCurve) {
        int hei = calclum->getHeight();
        int wid = calclum->getWidth();
        TMatrix wprofi = iccStore->workingSpaceMatrix (params->icm.working);

        const float wpi[3][3] = {
            {static_cast<float>(wprofi[0][0]), static_cast<float>(wprofi[0][1]), static_cast<float>(wprofi[0][2])},
            {static_cast<float>(wprofi[1][0]), static_cast<float>(wprofi[1][1]), static_cast<float>(wprofi[1][2])},
            {static_cast<float>(wprofi[2][0]), static_cast<float>(wprofi[2][1]), static_cast<float>(wprofi[2][2])}
        };
        lumcalcBuffer = new float[hei * wid];
        lumcalc = new float*[(hei)];

        for (int i = 0; i < hei; ++i) {
            lumcalc[i] = lumcalcBuffer + (i * wid);
        }

        ccalcBuffer = new float[hei * wid];
        ccalc = new float*[(hei)];

        for (int i = 0; i < hei; ++i) {
            ccalc[i] = ccalcBuffer + (i * wid);
        }

        float cn100Precalc;

        if (useNoiseCCurve) {
            cn100Precalc = SQR(1.f + ponderCC * (4.f * noiseCCurve[100.f / 60.f]));
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int ii = 0; ii < hei; ++ii) {
            for (int jj = 0; jj < wid; ++jj) {
                float LLum, AAum, BBum;

                float RL = calclum->r(ii, jj);
                float GL = calclum->g(ii, jj);
                float BL = calclum->b(ii, jj);
                // determine luminance and chrominance for noisecurves
                float XL, YL, ZL;
                Color::rgbxyz(RL, GL, BL, XL, YL, ZL, wpi);
                Color::XYZ2Lab(XL, YL, ZL, LLum, AAum, BBum);

                if (useNoiseLCurve) {
                    float epsi = 0.01f;

                    if (LLum < 2.f) {
                        LLum = 2.f;    //avoid divided by zero
                    }

                    if (LLum > 32768.f) {
                        LLum = 32768.f;    // not strictly necessary
                    }

                    float kinterm = epsi + noiseLCurve[xdivf(LLum, 15) * 500.f];
                    kinterm *= 100.f;
                    kinterm += noiseluma;
                    lumcalc[ii][jj] = SQR((kinterm / 125.f) * (1.f + kinterm / 25.f));
                }

                if (useNoiseCCurve) {
                    float cN = sqrtf(SQR(AAum) + SQR(BBum));

                    if (cN > 100) {
                        ccalc[ii][jj] = SQR(1.f + ponderCC * (4.f * noiseCCurve[cN / 60.f]));
                    } else {
                        ccalc[ii][jj] = cn100Precalc;
                    }
                }
            }
        }

        delete calclum;
        calclum = nullptr;
    }

    const short int imheight = src->getHeight(), imwidth = src->getWidth();

    if (dnparams.luma != 0 || dnparams.chroma != 0 || dnparams.methodmed == "Lab" || dnparams.methodmed == "Lonly") {
        // gamma transform for input data
        float gam = dnparams.gamma;
        float gamthresh = 0.001f;

        if (!isRAW) {//reduce gamma under 1 for Lab mode ==> TIF and JPG
            if (gam < 1.9f) {
                gam = 1.f - (1.9f - gam) / 3.f;    //minimum gamma 0.7
            } else if (gam >= 1.9f && gam <= 3.f) {
                gam = (1.4f / 1.1f) * gam - 1.41818f;
            }
        }


        LUTf gamcurve(65536, LUT_CLIP_BELOW);
        float gamslope = exp(log(static_cast<double>(gamthresh)) / gam) / gamthresh;

        if (denoiseMethodRgb) {
            Color::gammaf2lut(gamcurve, gam, gamthresh, gamslope, 65535.f, 32768.f);
        } else {
            Color::gammanf2lut(gamcurve, gam, 65535.f, 32768.f);
        }

        // inverse gamma transform for output data
        float igam = 1.f / gam;
        float igamthresh = gamthresh * gamslope;
        float igamslope = 1.f / gamslope;

        LUTf igamcurve(65536, LUT_CLIP_BELOW);

        if (denoiseMethodRgb) {
            Color::gammaf2lut(igamcurve, igam, igamthresh, igamslope, 32768.f, 65535.f);
        } else {
            Color::gammanf2lut(igamcurve, igam, 32768.f, 65535.f);
        }

        const float gain = pow (2.0f, float(expcomp));
        float noisevar_Ldetail = SQR(static_cast<float>(SQR(100. - dnparams.Ldetail) + 50.*(100. - dnparams.Ldetail)) * TS * 0.5f);

        array2D<float> tilemask_in(TS, TS);
        array2D<float> tilemask_out(TS, TS);

        if (denoiseLuminance) {
            const int border = MAX(2, TS / 16);

            for (int i = 0; i < TS; ++i) {
                float i1 = abs((i > TS / 2 ? i - TS + 1 : i));
                float vmask = (i1 < border ? SQR(sin((rtengine::RT_PI * i1) / (2 * border))) : 1.0f);
                float vmask2 = (i1 < 2 * border ? SQR(sin((rtengine::RT_PI * i1) / (2 * border))) : 1.0f);

                for (int j = 0; j < TS; ++j) {
                    float j1 = abs((j > TS / 2 ? j - TS + 1 : j));
                    tilemask_in[i][j] = (vmask * (j1 < border ? SQR(sin((rtengine::RT_PI * j1) / (2 * border))) : 1.0f)) + epsilon;
                    tilemask_out[i][j] = (vmask2 * (j1 < 2 * border ? SQR(sin((rtengine::RT_PI * j1) / (2 * border))) : 1.0f)) + epsilon;

                }
            }
        }

        int tilesize;
        int overlap;

        if (settings->leveldnti == 0) {
            tilesize = 1024;
            overlap = 128;
        }

        if (settings->leveldnti == 1) {
            tilesize = 768;
            overlap = 96;
        }

        int numTries = 0;

        if (ponder) {
            printf("Tiled denoise processing caused by Automatic Multizone mode\n");
        }

        bool memoryAllocationFailed = false;

        do {
            ++numTries;

            if (numTries == 2) {
                printf("1st denoise pass failed due to insufficient memory, starting 2nd (tiled) pass now...\n");
            }

            int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;

            Tile_calc (tilesize, overlap, (options.rgbDenoiseThreadLimit == 0 && !ponder) ? (numTries == 1 ? 0 : 2) : 2, imwidth, imheight, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);
            memoryAllocationFailed = false;
            const int numtiles = numtiles_W * numtiles_H;

            //output buffer
            Imagefloat * dsttmp;

            if (numtiles == 1) {
                dsttmp = dst;
            } else {
                dsttmp = new Imagefloat(imwidth, imheight);
#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int i = 0; i < imheight; ++i) {
                    for (int j = 0; j < imwidth; ++j) {
                        dsttmp->r(i, j) = 0.f;
                        dsttmp->g(i, j) = 0.f;
                        dsttmp->b(i, j) = 0.f;
                    }
                }
            }

            //now we have tile dimensions, overlaps
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            // According to FFTW-Doc 'it is safe to execute the same plan in parallel by multiple threads', so we now create 4 plans
            // outside the parallel region and use them inside the parallel region.

            // calculate max size of numblox_W.
            int max_numblox_W = ceil((static_cast<float>(MIN(imwidth, tilewidth))) / (offset)) + 2 * blkrad;
            // calculate min size of numblox_W.
            int min_numblox_W = ceil((static_cast<float>((MIN(imwidth, ((numtiles_W - 1) * tileWskip) + tilewidth)) - ((numtiles_W - 1) * tileWskip))) / (offset)) + 2 * blkrad;

            // these are needed only for creation of the plans and will be freed before entering the parallel loop
            fftwf_plan plan_forward_blox[2];
            fftwf_plan plan_backward_blox[2];

            if (denoiseLuminance) {
                float *Lbloxtmp  = reinterpret_cast<float*>( fftwf_malloc(max_numblox_W * TS * TS * sizeof (float)));
                float *fLbloxtmp = reinterpret_cast<float*>( fftwf_malloc(max_numblox_W * TS * TS * sizeof (float)));

                int nfwd[2] = {TS, TS};

                //for DCT:
                fftw_r2r_kind fwdkind[2] = {FFTW_REDFT10, FFTW_REDFT10};
                fftw_r2r_kind bwdkind[2] = {FFTW_REDFT01, FFTW_REDFT01};

                // Creating the plans with FFTW_MEASURE instead of FFTW_ESTIMATE speeds up the execute a bit
                plan_forward_blox[0]  = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, Lbloxtmp, nullptr, 1, TS * TS, fLbloxtmp, nullptr, 1, TS * TS, fwdkind, FFTW_MEASURE || FFTW_DESTROY_INPUT  );
                plan_backward_blox[0] = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, fLbloxtmp, nullptr, 1, TS * TS, Lbloxtmp, nullptr, 1, TS * TS, bwdkind, FFTW_MEASURE || FFTW_DESTROY_INPUT  );
                plan_forward_blox[1]  = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, Lbloxtmp, nullptr, 1, TS * TS, fLbloxtmp, nullptr, 1, TS * TS, fwdkind, FFTW_MEASURE || FFTW_DESTROY_INPUT  );
                plan_backward_blox[1] = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, fLbloxtmp, nullptr, 1, TS * TS, Lbloxtmp, nullptr, 1, TS * TS, bwdkind, FFTW_MEASURE || FFTW_DESTROY_INPUT  );
                fftwf_free (Lbloxtmp);
                fftwf_free (fLbloxtmp);
            }

#ifndef _OPENMP
            int numthreads = 1;
#else
            // Calculate number of tiles. If less than omp_get_max_threads(), then limit num_threads to number of tiles
            int numthreads = MIN(numtiles, omp_get_max_threads());

            if (options.rgbDenoiseThreadLimit > 0) {
                numthreads = MIN(numthreads, options.rgbDenoiseThreadLimit);
            }

#ifdef _RT_NESTED_OPENMP
            denoiseNestedLevels = omp_get_max_threads() / numthreads;
            bool oldNested = omp_get_nested();

            if (denoiseNestedLevels < 2) {
                denoiseNestedLevels = 1;
            } else {
                omp_set_nested(true);
            }

            if (options.rgbDenoiseThreadLimit > 0)
                while(denoiseNestedLevels * numthreads > options.rgbDenoiseThreadLimit) {
                    denoiseNestedLevels--;
                }

#endif

            if (settings->verbose) {
                printf("RGB_denoise uses %d main thread(s) and up to %d nested thread(s) for each main thread\n", numthreads, denoiseNestedLevels);
            }

#endif
            float *LbloxArray[denoiseNestedLevels * numthreads];
            float *fLbloxArray[denoiseNestedLevels * numthreads];

            if (numtiles > 1 && denoiseLuminance) {
                for (int i = 0; i < denoiseNestedLevels * numthreads; ++i) {
                    LbloxArray[i]  = reinterpret_cast<float*>( fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
                    fLbloxArray[i] = reinterpret_cast<float*>( fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
                }
            }

            TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
            //inverse matrix user select
            const float wip[3][3] = {
                {static_cast<float>(wiprof[0][0]), static_cast<float>(wiprof[0][1]), static_cast<float>(wiprof[0][2])},
                {static_cast<float>(wiprof[1][0]), static_cast<float>(wiprof[1][1]), static_cast<float>(wiprof[1][2])},
                {static_cast<float>(wiprof[2][0]), static_cast<float>(wiprof[2][1]), static_cast<float>(wiprof[2][2])}
            };

            TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);

            const float wp[3][3] = {
                {static_cast<float>(wprof[0][0]), static_cast<float>(wprof[0][1]), static_cast<float>(wprof[0][2])},
                {static_cast<float>(wprof[1][0]), static_cast<float>(wprof[1][1]), static_cast<float>(wprof[1][2])},
                {static_cast<float>(wprof[2][0]), static_cast<float>(wprof[2][1]), static_cast<float>(wprof[2][2])}
            };

            // begin tile processing of image
#ifdef _OPENMP
            #pragma omp parallel num_threads(numthreads) if (numthreads>1)
#endif
            {
                int pos;
                float* noisevarlum;
                float* noisevarchrom;

                if (numtiles == 1 && isRAW && (useNoiseCCurve || useNoiseLCurve)) {
                    noisevarlum = lumcalcBuffer;
                    noisevarchrom = ccalcBuffer;
                } else {
                    noisevarlum = new float[((tileheight + 1) / 2) * ((tilewidth + 1) / 2)];
                    noisevarchrom = new float[((tileheight + 1) / 2) * ((tilewidth + 1) / 2)];
                }

#ifdef _OPENMP
                #pragma omp for schedule(dynamic) collapse(2)
#endif

                for (int tiletop = 0; tiletop < imheight; tiletop += tileHskip) {
                    for (int tileleft = 0; tileleft < imwidth ; tileleft += tileWskip) {
                        //printf("titop=%d tileft=%d\n",tiletop/tileHskip, tileleft/tileWskip);
                        pos = (tiletop / tileHskip) * numtiles_W + tileleft / tileWskip ;
                        int tileright = MIN(imwidth, tileleft + tilewidth);
                        int tilebottom = MIN(imheight, tiletop + tileheight);
                        int width  = tileright - tileleft;
                        int height = tilebottom - tiletop;
                        int width2 = (width + 1) / 2;
                        float realred, realblue;
                        float interm_med = static_cast<float>( dnparams.chroma) / 10.0;
                        float intermred, intermblue;

                        if (dnparams.redchro > 0.) {
                            intermred = (dnparams.redchro / 10.);
                        } else {
                            intermred = static_cast<float>( dnparams.redchro) / 7.0;    //increase slower than linear for more sensit
                        }

                        if (dnparams.bluechro > 0.) {
                            intermblue = (dnparams.bluechro / 10.);
                        } else {
                            intermblue = static_cast<float>( dnparams.bluechro) / 7.0;    //increase slower than linear for more sensit
                        }

                        if (ponder && kall == 2) {
                            interm_med = ch_M[pos] / 10.f;
                            intermred = max_r[pos] / 10.f;
                            intermblue = max_b[pos] / 10.f;
                        }

                        if (ponder && kall == 0) {
                            interm_med = 0.01f;
                            intermred = 0.f;
                            intermblue = 0.f;
                        }

                        realred = interm_med + intermred;

                        if (realred <= 0.f) {
                            realred = 0.001f;
                        }

                        realblue = interm_med + intermblue;

                        if (realblue <= 0.f) {
                            realblue = 0.001f;
                        }

                        const float noisevarab_r = SQR(realred);
                        const float noisevarab_b = SQR(realblue);

                        //input L channel
                        array2D<float> *Lin = nullptr;
                        //wavelet denoised image
                        LabImage * labdn = new LabImage(width, height);

                        //fill tile from image; convert RGB to "luma/chroma"
                        const float maxNoiseVarab = max(noisevarab_b, noisevarab_r);

                        if (isRAW) {//image is raw; use channel differences for chroma channels

                            if (!denoiseMethodRgb) { //lab mode
                                //modification Jacques feb 2013 and july 2014
#ifdef _RT_NESTED_OPENMP
                                #pragma omp parallel for num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif

                                for (int i = tiletop; i < tilebottom; ++i) {
                                    int i1 = i - tiletop;

                                    for (int j = tileleft; j < tileright; ++j) {
                                        int j1 = j - tileleft;
                                        float R_ = gain * src->r(i, j);
                                        float G_ = gain * src->g(i, j);
                                        float B_ = gain * src->b(i, j);

                                        R_ = Color::denoiseIGammaTab[R_];
                                        G_ = Color::denoiseIGammaTab[G_];
                                        B_ = Color::denoiseIGammaTab[B_];

                                        //apply gamma noise standard (slider)
                                        R_ = R_ < 65535.f ? gamcurve[R_] : (Color::gammanf(R_ / 65535.f, gam) * 32768.f);
                                        G_ = G_ < 65535.f ? gamcurve[G_] : (Color::gammanf(G_ / 65535.f, gam) * 32768.f);
                                        B_ = B_ < 65535.f ? gamcurve[B_] : (Color::gammanf(B_ / 65535.f, gam) * 32768.f);

                                        //true conversion xyz=>Lab
                                        float X, Y, Z;
                                        Color::rgbxyz(R_, G_, B_, X, Y, Z, wp);

                                        //convert to Lab
                                        float L, a, b;
                                        Color::XYZ2Lab(X, Y, Z, L, a, b);

                                        labdn->L[i1][j1] = L;
                                        labdn->a[i1][j1] = a;
                                        labdn->b[i1][j1] = b;

                                        if (((i1 | j1) & 1) == 0) {
                                            if (numTries == 1) {
                                                noisevarlum[(i1 >> 1)*width2 + (j1 >> 1)] = useNoiseLCurve ? lumcalc[i >> 1][j >> 1] : noisevarL;
                                                noisevarchrom[(i1 >> 1)*width2 + (j1 >> 1)] = useNoiseCCurve ? maxNoiseVarab * ccalc[i >> 1][j >> 1] : 1.f;
                                            } else {
                                                noisevarlum[(i1 >> 1)*width2 + (j1 >> 1)] = lumcalc[i >> 1][j >> 1];
                                                noisevarchrom[(i1 >> 1)*width2 + (j1 >> 1)] = ccalc[i >> 1][j >> 1];
                                            }
                                        }

                                        //end chroma
                                    }
                                }
                            } else {//RGB mode
#ifdef _RT_NESTED_OPENMP
                                #pragma omp parallel for num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif

                                for (int i = tiletop; i < tilebottom; ++i) {
                                    int i1 = i - tiletop;

                                    for (int j = tileleft; j < tileright; ++j) {
                                        int j1 = j - tileleft;

                                        float X = gain * src->r(i, j);
                                        float Y = gain * src->g(i, j);
                                        float Z = gain * src->b(i, j);
                                        //conversion colorspace to determine luminance with no gamma
                                        X = X < 65535.f ? gamcurve[X] : (Color::gammaf(X / 65535.f, gam, gamthresh, gamslope) * 32768.f);
                                        Y = Y < 65535.f ? gamcurve[Y] : (Color::gammaf(Y / 65535.f, gam, gamthresh, gamslope) * 32768.f);
                                        Z = Z < 65535.f ? gamcurve[Z] : (Color::gammaf(Z / 65535.f, gam, gamthresh, gamslope) * 32768.f);
                                        //end chroma
                                        labdn->L[i1][j1] = Y;
                                        labdn->a[i1][j1] = (X - Y);
                                        labdn->b[i1][j1] = (Y - Z);

                                        if (((i1 | j1) & 1) == 0) {
                                            if (numTries == 1) {
                                                noisevarlum[(i1 >> 1)*width2 + (j1 >> 1)] = useNoiseLCurve ? lumcalc[i >> 1][j >> 1] : noisevarL;
                                                noisevarchrom[(i1 >> 1)*width2 + (j1 >> 1)] = useNoiseCCurve ? maxNoiseVarab * ccalc[i >> 1][j >> 1] : 1.f;
                                            } else {
                                                noisevarlum[(i1 >> 1)*width2 + (j1 >> 1)] = lumcalc[i >> 1][j >> 1];
                                                noisevarchrom[(i1 >> 1)*width2 + (j1 >> 1)] = ccalc[i >> 1][j >> 1];
                                            }
                                        }
                                    }
                                }
                            }
                        } else {//image is not raw; use Lab parametrization
#ifdef _RT_NESTED_OPENMP
                            #pragma omp parallel for num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif

                            for (int i = tiletop; i < tilebottom; ++i) {
                                int i1 = i - tiletop;

                                for (int j = tileleft; j < tileright; ++j) {
                                    int j1 = j - tileleft;
                                    float L, a, b;
                                    float rLum = src->r(i, j) ; //for denoise curves
                                    float gLum = src->g(i, j) ;
                                    float bLum = src->b(i, j) ;

                                    //use gamma sRGB, not good if TIF (JPG) Output profil not with gamma sRGB  (eg : gamma =1.0, or 1.8...)
                                    //very difficult to solve !
                                    // solution ==> save TIF with gamma sRGB and re open
                                    float rtmp = Color::igammatab_srgb[ src->r(i, j) ];
                                    float gtmp = Color::igammatab_srgb[ src->g(i, j) ];
                                    float btmp = Color::igammatab_srgb[ src->b(i, j) ];
                                    //modification Jacques feb 2013
                                    // gamma slider different from raw
                                    rtmp = rtmp < 65535.f ? gamcurve[rtmp] : (Color::gammanf(rtmp / 65535.f, gam) * 32768.f);
                                    gtmp = gtmp < 65535.f ? gamcurve[gtmp] : (Color::gammanf(gtmp / 65535.f, gam) * 32768.f);
                                    btmp = btmp < 65535.f ? gamcurve[btmp] : (Color::gammanf(btmp / 65535.f, gam) * 32768.f);

                                    float X, Y, Z;
                                    Color::rgbxyz(rtmp, gtmp, btmp, X, Y, Z, wp);

                                    //convert Lab
                                    Color::XYZ2Lab(X, Y, Z, L, a, b);
                                    labdn->L[i1][j1] = L;
                                    labdn->a[i1][j1] = a;
                                    labdn->b[i1][j1] = b;

                                    if (((i1 | j1) & 1) == 0) {
                                        float Llum, alum, blum;

                                        if (useNoiseLCurve || useNoiseCCurve) {
                                            float XL, YL, ZL;
                                            Color::rgbxyz(rLum, gLum, bLum, XL, YL, ZL, wp);
                                            Color::XYZ2Lab(XL, YL, ZL, Llum, alum, blum);
                                        }

                                        if (useNoiseLCurve) {
                                            float kN = Llum;
                                            float epsi = 0.01f;

                                            if (kN < 2.f) {
                                                kN = 2.f;
                                            }

                                            if (kN > 32768.f) {
                                                kN = 32768.f;
                                            }

                                            float kinterm = epsi + noiseLCurve[xdivf(kN, 15) * 500.f];
                                            float ki = kinterm * 100.f;
                                            ki += noiseluma;
                                            noisevarlum[(i1 >> 1)*width2 + (j1 >> 1)] = SQR((ki / 125.f) * (1.f + ki / 25.f));
                                        } else {
                                            noisevarlum[(i1 >> 1)*width2 + (j1 >> 1)] = noisevarL;
                                        }

                                        if (useNoiseCCurve) {
                                            float aN = alum;
                                            float bN = blum;
                                            float cN = sqrtf(SQR(aN) + SQR(bN));

                                            if (cN < 100.f) {
                                                cN = 100.f;    //avoid divided by zero ???
                                            }

                                            float Cinterm = 1.f + ponderCC * 4.f * noiseCCurve[cN / 60.f];
                                            noisevarchrom[(i1 >> 1)*width2 + (j1 >> 1)] = maxNoiseVarab * SQR(Cinterm);
                                        } else {
                                            noisevarchrom[(i1 >> 1)*width2 + (j1 >> 1)] = 1.f;
                                        }
                                    }
                                }
                            }
                        }

                        //now perform basic wavelet denoise
                        //arguments 4 and 5 of wavelet decomposition are max number of wavelet decomposition levels;
                        //and whether to subsample the image after wavelet filtering.  Subsampling is coded as
                        //binary 1 or 0 for each level, eg subsampling = 0 means no subsampling, 1 means subsample
                        //the first level only, 7 means subsample the first three levels, etc.
                        //actual implementation only works with subsampling set to 1
                        float interm_medT = static_cast<float>( dnparams.chroma) / 10.0;
                        bool execwavelet = true;

                        if (!denoiseLuminance && interm_medT < 0.05f && dnparams.median && (dnparams.methodmed == "Lab" || dnparams.methodmed == "Lonly")) {
                            execwavelet = false;    //do not exec wavelet if sliders luminance and chroma are very small and median need
                        }

                        //we considered user don't want wavelet
                        if (settings->leveldnautsimpl == 1 && dnparams.Cmethod != "MAN") {
                            execwavelet = true;
                        }

                        if (settings->leveldnautsimpl == 0 && dnparams.C2method != "MANU") {
                            execwavelet = true;
                        }

                        if (execwavelet) {//gain time if user choose only median  sliders L <=1  slider chrom master < 1
                            wavelet_decomposition* Ldecomp;
                            wavelet_decomposition* adecomp;

                            int levwav = 5;
                            float maxreal = max(realred, realblue);

                            //increase the level of wavelet if user increase much or very much sliders
                            if (maxreal < 8.f) {
                                levwav = 5;
                            } else if (maxreal < 10.f) {
                                levwav = 6;
                            } else if (maxreal < 15.f) {
                                levwav = 7;
                            } else {
                                levwav = 8;    //maximum ==> I have increase Maxlevel in cplx_wavelet_dec.h from 8 to 9
                            }

                            if (nrQuality == QUALITY_HIGH) {
                                levwav += settings->nrwavlevel;    //increase level for enhanced mode
                            }

                            if (levwav > 8) {
                                levwav = 8;
                            }

                            int minsizetile = min(tilewidth, tileheight);
                            int maxlev2 = 8;

                            if (minsizetile < 256) {
                                maxlev2 = 7;
                            }

                            if (minsizetile < 128) {
                                maxlev2 = 6;
                            }

                            if (minsizetile < 64) {
                                maxlev2 = 5;
                            }

                            levwav = min(maxlev2, levwav);

                            //  if (settings->verbose) printf("levwavelet=%i  noisevarA=%f noisevarB=%f \n",levwav, noisevarab_r, noisevarab_b);
                            Ldecomp = new wavelet_decomposition (labdn->L[0], labdn->W, labdn->H, levwav, 1, 1, max(1, denoiseNestedLevels));

                            if (Ldecomp->memoryAllocationFailed) {
                                memoryAllocationFailed = true;
                            }

                            float madL[8][3];

                            if (!memoryAllocationFailed) {
                                // precalculate madL, because it's used in adecomp and bdecomp
                                int maxlvl = Ldecomp->maxlevel();
#ifdef _RT_NESTED_OPENMP
                                #pragma omp parallel for schedule(dynamic) collapse(2) num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif

                                for (int lvl = 0; lvl < maxlvl; ++lvl) {
                                    for (int dir = 1; dir < 4; ++dir) {
                                        // compute median absolute deviation (MAD) of detail coefficients as robust noise estimator
                                        int Wlvl_L = Ldecomp->level_W(lvl);
                                        int Hlvl_L = Ldecomp->level_H(lvl);

                                        float ** WavCoeffs_L = Ldecomp->level_coeffs(lvl);

                                        if (!denoiseMethodRgb) {
                                            madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                                        } else {
                                            madL[lvl][dir - 1] = SQR(MadRgb(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                                        }
                                    }
                                }
                            }

                            float chresid = 0.f;
                            float chresidtemp = 0.f;
                            float chmaxresid = 0.f;
                            float chmaxresidtemp = 0.f;

                            adecomp = new wavelet_decomposition (labdn->a[0], labdn->W, labdn->H, levwav, 1, 1, max(1, denoiseNestedLevels));

                            if (adecomp->memoryAllocationFailed) {
                                memoryAllocationFailed = true;
                            }

                            if (!memoryAllocationFailed) {
                                if (nrQuality == QUALITY_STANDARD) {
                                    if (!WaveletDenoiseAllAB(*Ldecomp, *adecomp, noisevarchrom, madL, noisevarab_r, useNoiseCCurve, autoch, denoiseMethodRgb)) { //enhance mode
                                        memoryAllocationFailed = true;
                                    }
                                } else { /*if (nrQuality==QUALITY_HIGH)*/
                                    if (!WaveletDenoiseAll_BiShrinkAB(*Ldecomp, *adecomp, noisevarchrom, madL, noisevarab_r, useNoiseCCurve, autoch, denoiseMethodRgb)) { //enhance mode
                                        memoryAllocationFailed = true;
                                    }

                                    if (!memoryAllocationFailed) {
                                        if (!WaveletDenoiseAllAB(*Ldecomp, *adecomp, noisevarchrom, madL, noisevarab_r, useNoiseCCurve, autoch, denoiseMethodRgb)) {
                                            memoryAllocationFailed = true;
                                        }
                                    }
                                }
                            }

                            if (!memoryAllocationFailed) {
                                if (kall == 0) {
                                    Noise_residualAB(*adecomp, chresid, chmaxresid, denoiseMethodRgb);
                                    chresidtemp = chresid;
                                    chmaxresidtemp = chmaxresid;
                                }

                                adecomp->reconstruct(labdn->a[0]);
                            }

                            delete adecomp;

                            if (!memoryAllocationFailed) {
                                wavelet_decomposition* bdecomp = new wavelet_decomposition (labdn->b[0], labdn->W, labdn->H, levwav, 1, 1, max(1, denoiseNestedLevels));

                                if (bdecomp->memoryAllocationFailed) {
                                    memoryAllocationFailed = true;
                                }

                                if (!memoryAllocationFailed) {
                                    if (nrQuality == QUALITY_STANDARD) {
                                        if (!WaveletDenoiseAllAB(*Ldecomp, *bdecomp, noisevarchrom, madL, noisevarab_b, useNoiseCCurve, autoch, denoiseMethodRgb)) { //enhance mode
                                            memoryAllocationFailed = true;
                                        }
                                    } else { /*if (nrQuality==QUALITY_HIGH)*/
                                        if (!WaveletDenoiseAll_BiShrinkAB(*Ldecomp, *bdecomp, noisevarchrom, madL, noisevarab_b, useNoiseCCurve, autoch, denoiseMethodRgb)) { //enhance mode
                                            memoryAllocationFailed = true;
                                        }

                                        if (!memoryAllocationFailed) {
                                            if (!WaveletDenoiseAllAB(*Ldecomp, *bdecomp, noisevarchrom, madL, noisevarab_b, useNoiseCCurve, autoch, denoiseMethodRgb)) {
                                                memoryAllocationFailed = true;
                                            }
                                        }
                                    }
                                }

                                if (!memoryAllocationFailed) {
                                    if (kall == 0) {
                                        Noise_residualAB(*bdecomp, chresid, chmaxresid, denoiseMethodRgb);
                                        chresid += chresidtemp;
                                        chmaxresid += chmaxresidtemp;
                                        chresid = sqrt(chresid / (6 * (levwav)));
                                        highresi = chresid + 0.66f * (sqrt(chmaxresid) - chresid); //evaluate sigma
                                        nresi = chresid;
                                    }

                                    bdecomp->reconstruct(labdn->b[0]);
                                }

                                delete bdecomp;

                                if (!memoryAllocationFailed) {
                                    if (denoiseLuminance) {
                                        int edge = 0;

                                        if (nrQuality == QUALITY_STANDARD) {
                                            if (!WaveletDenoiseAllL(*Ldecomp, noisevarlum, madL, nullptr, edge)) { //enhance mode
                                                memoryAllocationFailed = true;
                                            }
                                        } else { /*if (nrQuality==QUALITY_HIGH)*/
                                            if (!WaveletDenoiseAll_BiShrinkL(*Ldecomp, noisevarlum, madL)) { //enhance mode
                                                memoryAllocationFailed = true;
                                            }

                                            if (!memoryAllocationFailed) {
                                                if (!WaveletDenoiseAllL(*Ldecomp, noisevarlum, madL, nullptr, edge)) {
                                                    memoryAllocationFailed = true;
                                                }
                                            }
                                        }

                                        if (!memoryAllocationFailed) {
                                            // copy labdn->L to Lin before it gets modified by reconstruction
                                            Lin = new array2D<float>(width, height);
#ifdef _RT_NESTED_OPENMP
                                            #pragma omp parallel for num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif

                                            for (int i = 0; i < height; ++i) {
                                                for (int j = 0; j < width; ++j) {
                                                    (*Lin)[i][j] = labdn->L[i][j];
                                                }
                                            }

                                            Ldecomp->reconstruct(labdn->L[0]);
                                        }
                                    }
                                }
                            }

                            delete Ldecomp;
                        }

                        if (!memoryAllocationFailed) {
                            if ((metchoice == 1 || metchoice == 2 || metchoice == 3 || metchoice == 4) && dnparams.median) {
                                float** tmL;
                                int wid = labdn->W;
                                int hei = labdn->H;
                                tmL = new float*[hei];

                                for (int i = 0; i < hei; ++i) {
                                    tmL[i] = new float[wid];
                                }

                                Median medianTypeL = Median::TYPE_3X3_SOFT;
                                Median medianTypeAB = Median::TYPE_3X3_SOFT;

                                if (dnparams.medmethod == "soft") {
                                    if (metchoice != 4) {
                                        medianTypeL = medianTypeAB = Median::TYPE_3X3_SOFT;
                                    } else {
                                        medianTypeL = Median::TYPE_3X3_SOFT;
                                        medianTypeAB = Median::TYPE_3X3_SOFT;
                                    }
                                } else if (dnparams.medmethod == "33") {
                                    if (metchoice != 4) {
                                        medianTypeL = medianTypeAB = Median::TYPE_3X3_STRONG;
                                    } else {
                                        medianTypeL = Median::TYPE_3X3_SOFT;
                                        medianTypeAB = Median::TYPE_3X3_STRONG;
                                    }
                                } else if (dnparams.medmethod == "55soft") {
                                    if (metchoice != 4) {
                                        medianTypeL = medianTypeAB = Median::TYPE_5X5_SOFT;
                                    } else {
                                        medianTypeL = Median::TYPE_3X3_SOFT;
                                        medianTypeAB = Median::TYPE_5X5_SOFT;
                                    }
                                } else if (dnparams.medmethod == "55") {
                                    if (metchoice != 4) {
                                        medianTypeL = medianTypeAB = Median::TYPE_5X5_STRONG;
                                    } else {
                                        medianTypeL = Median::TYPE_3X3_STRONG;
                                        medianTypeAB = Median::TYPE_5X5_STRONG;
                                    }
                                } else if (dnparams.medmethod == "77") {
                                    if (metchoice != 4) {
                                        medianTypeL = medianTypeAB = Median::TYPE_7X7;
                                    } else {
                                        medianTypeL = Median::TYPE_3X3_STRONG;
                                        medianTypeAB = Median::TYPE_7X7;
                                    }
                                } else if (dnparams.medmethod == "99") {
                                    if (metchoice != 4) {
                                        medianTypeL = medianTypeAB = Median::TYPE_9X9;
                                    } else {
                                        medianTypeL = Median::TYPE_5X5_SOFT;
                                        medianTypeAB = Median::TYPE_9X9;
                                    }
                                }

                                if (metchoice == 1 || metchoice == 2 || metchoice == 4) {
                                    Median_Denoise(labdn->L, labdn->L, wid, hei, medianTypeL, dnparams.passes, denoiseNestedLevels, tmL);
                                }

                                if (metchoice == 2 || metchoice == 3 || metchoice == 4) {
                                    Median_Denoise(labdn->a, labdn->a, wid, hei, medianTypeAB, dnparams.passes, denoiseNestedLevels, tmL);
                                    Median_Denoise(labdn->b, labdn->b, wid, hei, medianTypeAB, dnparams.passes, denoiseNestedLevels, tmL);
                                }

                                for (int i = 0; i < hei; ++i) {
                                    delete[] tmL[i];
                                }

                                delete[] tmL;
                            }

                            //wavelet denoised L channel
                            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            // now do detail recovery using block DCT to detect
                            // patterns missed by wavelet denoise
                            // blocks are not the same thing as tiles!

                            // calculation for detail recovery blocks
                            const int numblox_W = ceil((static_cast<float>(width)) / (offset)) + 2 * blkrad;
                            const int numblox_H = ceil((static_cast<float>(height)) / (offset)) + 2 * blkrad;



                            // end of tiling calc

                            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            // Main detail recovery algorithm: Block loop
                            //DCT block data storage

                            if (denoiseLuminance /*&& execwavelet*/) {
                                //residual between input and denoised L channel
                                array2D<float> Ldetail(width, height, ARRAY2D_CLEAR_DATA);
                                //pixel weight
                                array2D<float> totwt(width, height, ARRAY2D_CLEAR_DATA); //weight for combining DCT blocks

                                if (numtiles == 1) {
                                    for (int i = 0; i < denoiseNestedLevels * numthreads; ++i) {
                                        LbloxArray[i]  = reinterpret_cast<float*>( fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
                                        fLbloxArray[i] = reinterpret_cast<float*>( fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
                                    }
                                }

#ifdef _RT_NESTED_OPENMP
                                int masterThread = omp_get_thread_num();
#endif
#ifdef _RT_NESTED_OPENMP
                                #pragma omp parallel num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif
                                {
#ifdef _RT_NESTED_OPENMP
                                    int subThread = masterThread * denoiseNestedLevels + omp_get_thread_num();
#else
                                    int subThread = 0;
#endif
                                    float blurbuffer[TS * TS] ALIGNED64;
                                    float *Lblox = LbloxArray[subThread];
                                    float *fLblox = fLbloxArray[subThread];
                                    float pBuf[width + TS + 2 * blkrad * offset] ALIGNED16;
                                    float nbrwt[TS * TS] ALIGNED64;
#ifdef _RT_NESTED_OPENMP
                                    #pragma omp for
#endif

                                    for (int vblk = 0; vblk < numblox_H; ++vblk) {

                                        int top = (vblk - blkrad) * offset;
                                        float * datarow = pBuf + blkrad * offset;

                                        for (int i = 0; i < TS; ++i) {
                                            int row = top + i;
                                            int rr = row;

                                            if (row < 0) {
                                                rr = MIN(-row, height - 1);
                                            } else if (row >= height) {
                                                rr = MAX(0, 2 * height - 2 - row);
                                            }

                                            for (int j = 0; j < labdn->W; ++j) {
                                                datarow[j] = ((*Lin)[rr][j] - labdn->L[rr][j]);
                                            }

                                            for (int j = -blkrad * offset; j < 0; ++j) {
                                                datarow[j] = datarow[MIN(-j, width - 1)];
                                            }

                                            for (int j = width; j < width + TS + blkrad * offset; ++j) {
                                                datarow[j] = datarow[MAX(0, 2 * width - 2 - j)];
                                            }//now we have a padded data row

                                            //now fill this row of the blocks with Lab high pass data
                                            for (int hblk = 0; hblk < numblox_W; ++hblk) {
                                                int left = (hblk - blkrad) * offset;
                                                int indx = (hblk) * TS; //index of block in malloc

                                                if (top + i >= 0 && top + i < height) {
                                                    int j;

                                                    for (j = 0; j < min((-left), TS); ++j) {
                                                        Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                                                    }

                                                    for (; j < min(TS, width - left); ++j) {
                                                        Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                                                        totwt[top + i][left + j] += tilemask_in[i][j] * tilemask_out[i][j];
                                                    }

                                                    for (; j < TS; ++j) {
                                                        Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                                                    }
                                                } else {
                                                    for (int j = 0; j < TS; ++j) {
                                                        Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                                                    }
                                                }

                                            }

                                        }//end of filling block row

                                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        //fftwf_print_plan (plan_forward_blox);
                                        if (numblox_W == max_numblox_W) {
                                            fftwf_execute_r2r(plan_forward_blox[0], Lblox, fLblox);    // DCT an entire row of tiles
                                        } else {
                                            fftwf_execute_r2r(plan_forward_blox[1], Lblox, fLblox);    // DCT an entire row of tiles
                                        }

                                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        // now process the vblk row of blocks for noise reduction


                                        for (int hblk = 0; hblk < numblox_W; ++hblk) {
                                            RGBtile_denoise (fLblox, hblk, noisevar_Ldetail, nbrwt, blurbuffer);
                                        }//end of horizontal block loop

                                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                        //now perform inverse FT of an entire row of blocks
                                        if (numblox_W == max_numblox_W) {
                                            fftwf_execute_r2r(plan_backward_blox[0], fLblox, Lblox);    //for DCT
                                        } else {
                                            fftwf_execute_r2r(plan_backward_blox[1], fLblox, Lblox);    //for DCT
                                        }

                                        int topproc = (vblk - blkrad) * offset;

                                        //add row of blocks to output image tile
                                        RGBoutput_tile_row (Lblox, Ldetail, tilemask_out, height, width, topproc);

                                        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                    }//end of vertical block loop

                                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                }
                                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef _RT_NESTED_OPENMP
                                #pragma omp parallel for num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif

                                for (int i = 0; i < height; ++i) {
                                    for (int j = 0; j < width; ++j) {
                                        //may want to include masking threshold for large hipass data to preserve edges/detail
                                        labdn->L[i][j] += Ldetail[i][j] / totwt[i][j]; //note that labdn initially stores the denoised hipass data
                                    }
                                }
                            }

                            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            // transform denoised "Lab" to output RGB

                            //calculate mask for feathering output tile overlaps
                            float Vmask[height + 1] ALIGNED16;
                            float Hmask[width + 1] ALIGNED16;
                            float newGain;

                            if (numtiles > 1) {
                                for (int i = 0; i < height; ++i) {
                                    Vmask[i] = 1;
                                }

                                newGain = 1.f;

                                if (isRAW) {
                                    newGain = gain;
                                }

                                for (int j = 0; j < width; ++j) {
                                    Hmask[j] = 1.f / newGain;
                                }

                                for (int i = 0; i < overlap; ++i) {
                                    float mask = SQR(xsinf((rtengine::RT_PI * i) / (2 * overlap)));

                                    if (tiletop > 0) {
                                        Vmask[i] = mask;
                                    }

                                    if (tilebottom < imheight) {
                                        Vmask[height - i] = mask;
                                    }

                                    if (tileleft > 0) {
                                        Hmask[i] = mask / newGain;
                                    }

                                    if (tileright < imwidth) {
                                        Hmask[width - i] = mask / newGain;
                                    }
                                }
                            } else {
                                newGain = isRAW ? 1.f / gain : 1.f;;
                            }

                            //convert back to RGB and write to destination array
                            if (isRAW) {
                                if (!denoiseMethodRgb) {//Lab mode
                                    realred /= 100.f;
                                    realblue /= 100.f;

#ifdef _RT_NESTED_OPENMP
                                    #pragma omp parallel for schedule(dynamic,16) num_threads(denoiseNestedLevels)
#endif

                                    for (int i = tiletop; i < tilebottom; ++i) {
                                        int i1 = i - tiletop;

                                        for (int j = tileleft; j < tileright; ++j) {
                                            int j1 = j - tileleft;
                                            //modification Jacques feb 2013
                                            //true conversion Lab==>xyz
                                            float L = labdn->L[i1][j1];
                                            float a = labdn->a[i1][j1];
                                            float b = labdn->b[i1][j1];
                                            float c_h = SQR(a) + SQR(b);

                                            if (c_h > 9000000.f) {
                                                a *= 1.f + qhighFactor * realred;
                                                b *= 1.f + qhighFactor * realblue;
                                            }

                                            //convert XYZ
                                            float X, Y, Z;
                                            Color::Lab2XYZ(L, a, b, X, Y, Z);
                                            //apply inverse gamma noise
                                            float r_, g_, b_;
                                            Color::xyz2rgb(X, Y, Z, r_, g_, b_, wip);
                                            //inverse gamma standard (slider)
                                            r_ = r_ < 32768.f ? igamcurve[r_] : (Color::gammanf(r_ / 32768.f, igam) * 65535.f);
                                            g_ = g_ < 32768.f ? igamcurve[g_] : (Color::gammanf(g_ / 32768.f, igam) * 65535.f);
                                            b_ = b_ < 32768.f ? igamcurve[b_] : (Color::gammanf(b_ / 32768.f, igam) * 65535.f);

                                            //readapt arbitrary gamma (inverse from beginning)
                                            r_ = Color::denoiseGammaTab[r_];
                                            g_ = Color::denoiseGammaTab[g_];
                                            b_ = Color::denoiseGammaTab[b_];

                                            if (numtiles == 1) {
                                                dsttmp->r(i, j) = newGain * r_;
                                                dsttmp->g(i, j) = newGain * g_;
                                                dsttmp->b(i, j) = newGain * b_;
                                            } else {
                                                float factor = Vmask[i1] * Hmask[j1];
                                                dsttmp->r(i, j) += factor * r_;
                                                dsttmp->g(i, j) += factor * g_;
                                                dsttmp->b(i, j) += factor * b_;
                                            }
                                        }
                                    }
                                } else {//RGB mode
#ifdef _RT_NESTED_OPENMP
                                    #pragma omp parallel for num_threads(denoiseNestedLevels)
#endif

                                    for (int i = tiletop; i < tilebottom; ++i) {
                                        int i1 = i - tiletop;

                                        for (int j = tileleft; j < tileright; ++j) {
                                            int j1 = j - tileleft;
                                            float c_h = sqrt(SQR(labdn->a[i1][j1]) + SQR(labdn->b[i1][j1]));

                                            if (c_h > 3000.f) {
                                                labdn->a[i1][j1] *= 1.f + qhighFactor * realred / 100.f;
                                                labdn->b[i1][j1] *= 1.f + qhighFactor * realblue / 100.f;
                                            }

                                            float Y = labdn->L[i1][j1];
                                            float X = (labdn->a[i1][j1]) + Y;
                                            float Z = Y - (labdn->b[i1][j1]);


                                            X = X < 32768.f ? igamcurve[X] : (Color::gammaf(X / 32768.f, igam, igamthresh, igamslope) * 65535.f);
                                            Y = Y < 32768.f ? igamcurve[Y] : (Color::gammaf(Y / 32768.f, igam, igamthresh, igamslope) * 65535.f);
                                            Z = Z < 32768.f ? igamcurve[Z] : (Color::gammaf(Z / 32768.f, igam, igamthresh, igamslope) * 65535.f);

                                            if (numtiles == 1) {
                                                dsttmp->r(i, j) = newGain * X;
                                                dsttmp->g(i, j) = newGain * Y;
                                                dsttmp->b(i, j) = newGain * Z;
                                            } else {
                                                float factor = Vmask[i1] * Hmask[j1];
                                                dsttmp->r(i, j) += factor * X;
                                                dsttmp->g(i, j) += factor * Y;
                                                dsttmp->b(i, j) += factor * Z;
                                            }
                                        }
                                    }

                                }
                            } else {
#ifdef _RT_NESTED_OPENMP
                                #pragma omp parallel for num_threads(denoiseNestedLevels)
#endif

                                for (int i = tiletop; i < tilebottom; ++i) {
                                    int i1 = i - tiletop;

                                    for (int j = tileleft; j < tileright; ++j) {
                                        int j1 = j - tileleft;
                                        //modification Jacques feb 2013
                                        float L = labdn->L[i1][j1];
                                        float a = labdn->a[i1][j1];
                                        float b = labdn->b[i1][j1];
                                        float c_h = sqrt(SQR(a) + SQR(b));

                                        if (c_h > 3000.f) {
                                            a *= 1.f + qhighFactor * realred / 100.f;
                                            b *= 1.f + qhighFactor * realblue / 100.f;
                                        }

                                        float X, Y, Z;
                                        Color::Lab2XYZ(L, a, b, X, Y, Z);

                                        float r_, g_, b_;
                                        Color::xyz2rgb(X, Y, Z, r_, g_, b_, wip);
                                        //gamma slider is different from Raw
                                        r_ = r_ < 32768.f ? igamcurve[r_] : (Color::gammanf(r_ / 32768.f, igam) * 65535.f);
                                        g_ = g_ < 32768.f ? igamcurve[g_] : (Color::gammanf(g_ / 32768.f, igam) * 65535.f);
                                        b_ = b_ < 32768.f ? igamcurve[b_] : (Color::gammanf(b_ / 32768.f, igam) * 65535.f);

                                        if (numtiles == 1) {
                                            dsttmp->r(i, j) = newGain * r_;
                                            dsttmp->g(i, j) = newGain * g_;
                                            dsttmp->b(i, j) = newGain * b_;
                                        } else {
                                            float factor = Vmask[i1] * Hmask[j1];
                                            dsttmp->r(i, j) += factor * r_;
                                            dsttmp->g(i, j) += factor * g_;
                                            dsttmp->b(i, j) += factor * b_;
                                        }
                                    }
                                }
                            }

                            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        }

                        delete labdn;
                        delete Lin;

                    }//end of tile row
                }//end of tile loop

                if (numtiles > 1 || !isRAW || (!useNoiseCCurve && !useNoiseLCurve)) {
                    delete[] noisevarlum;
                    delete[] noisevarchrom;
                }

            }

            if (denoiseLuminance) {
                for (int i = 0; i < denoiseNestedLevels * numthreads; ++i) {
                    fftwf_free(LbloxArray[i]);
                    fftwf_free(fLbloxArray[i]);
                }
            }

#ifdef _RT_NESTED_OPENMP
            omp_set_nested(oldNested);
#endif

            //copy denoised image to output
            if (numtiles > 1) {
                if (!memoryAllocationFailed) {
                    dsttmp->copyData(dst);
                } else if (dst != src) {
                    src->copyData(dst);
                }

                delete dsttmp;
            }

            if (!isRAW && !memoryAllocationFailed) {//restore original image gamma
#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int i = 0; i < dst->getHeight(); ++i) {
                    for (int j = 0; j < dst->getWidth(); ++j) {
                        dst->r(i, j) = Color::gammatab_srgb[ dst->r(i, j) ];
                        dst->g(i, j) = Color::gammatab_srgb[ dst->g(i, j) ];
                        dst->b(i, j) = Color::gammatab_srgb[ dst->b(i, j) ];
                    }
                }
            }

            if (denoiseLuminance) {
                // destroy the plans
                fftwf_destroy_plan(plan_forward_blox[0]);
                fftwf_destroy_plan(plan_backward_blox[0]);
                fftwf_destroy_plan(plan_forward_blox[1]);
                fftwf_destroy_plan(plan_backward_blox[1]);
                fftwf_cleanup();
            }
        } while(memoryAllocationFailed && numTries < 2 && (options.rgbDenoiseThreadLimit == 0) && !ponder);

        if (memoryAllocationFailed) {
            printf("tiled denoise failed due to isufficient memory. Output is not denoised!\n");
        }

    }


//median 3x3 in complement on RGB
    if (dnparams.methodmed == "RGB" && dnparams.median) {
//printf("RGB den\n");
        int wid = dst->getWidth(), hei = dst->getHeight();
        float** tm;
        tm = new float*[hei];

        for (int i = 0; i < hei; ++i) {
            tm[i] = new float[wid];
        }

        Imagefloat *source;

        if (dnparams.luma == 0 && dnparams.chroma == 0) {
            source = dst;
        } else {
            source = src;
        }

        int methmed = 0;
        int border = 1;

        if (dnparams.rgbmethod == "soft") {
            methmed = 0;
        } else if (dnparams.rgbmethod == "33") {
            methmed = 1;
        } else if (dnparams.rgbmethod == "55") {
            methmed = 3;
            border = 2;
        } else if (dnparams.rgbmethod == "55soft") {
            methmed = 2;
            border = 2;
        }

        for (int iteration = 1; iteration <= dnparams.passes; ++iteration) {

            #pragma omp parallel
            {
                if (methmed < 2)
                {
                    #pragma omp for

                    for (int i = 1; i < hei - 1; ++i) {
                        if (methmed == 0) {
                            for (int j = 1; j < wid - 1; ++j) {
                                tm[i][j] = median(source->r(i, j), source->r(i - 1, j), source->r(i + 1, j), source->r(i, j + 1), source->r(i, j - 1)); //3x3 soft
                            }
                        } else {
                            for (int j = 1; j < wid - 1; ++j) {
                                tm[i][j] = median(source->r(i, j), source->r(i - 1, j), source->r(i + 1, j), source->r(i, j + 1), source->r(i, j - 1), source->r(i - 1, j - 1), source->r(i - 1, j + 1), source->r(i + 1, j - 1), source->r(i + 1, j + 1)); //3x3
                            }
                        }
                    }
                } else {
                    #pragma omp for

                    for (int i = 2; i < hei - 2; ++i)
                    {
                        if (methmed == 3) {
                            for (int j = 2; j < wid - 2; ++j) {
                                tm[i][j] = median(source->r(i, j), source->r(i - 1, j), source->r(i + 1, j), source->r(i, j + 1), source->r(i, j - 1), source->r(i - 1, j - 1), source->r(i - 1, j + 1), source->r(i + 1, j - 1), source->r(i + 1, j + 1),
                                source->r(i - 2, j), source->r(i + 2, j), source->r(i, j + 2), source->r(i, j - 2), source->r(i - 2, j - 2), source->r(i - 2, j + 2), source->r(i + 2, j - 2), source->r(i + 2, j + 2),
                                source->r(i - 2, j + 1), source->r(i + 2, j + 1), source->r(i - 1, j + 2), source->r(i - 1, j - 2), source->r(i - 2, j - 1), source->r(i + 2, j - 1), source->r(i + 1, j + 2), source->r(i + 1, j - 2));//5x5
                            }
                        } else {
                            for (int j = 2; j < wid - 2; ++j) {
                                tm[i][j] = median(
                                               source->r(i, j),
                                               source->r(i - 1, j),
                                               source->r(i + 1, j),
                                               source->r(i, j + 1),
                                               source->r(i, j - 1),
                                               source->r(i - 1, j - 1),
                                               source->r(i - 1, j + 1),
                                               source->r(i + 1, j - 1),
                                               source->r(i + 1, j + 1),
                                               source->r(i + 2, j),
                                               source->r(i - 2, j),
                                               source->r(i, j + 2),
                                               source->r(i, j - 2)
                                           ); // 5x5 soft
                            }
                        }
                    }
                }
#ifdef _OPENMP
                #pragma omp for nowait
#endif

                for (int i = border; i < hei - border; ++i)
                {
                    for (int j = border; j < wid - border; ++j) {
                        dst->r(i, j) = tm[i][j];
                    }
                }

                if (methmed < 2)
                {
                    #pragma omp for

                    for (int i = 1; i < hei - 1; ++i) {
                        if (methmed == 0) {
                            for (int j = 1; j < wid - 1; ++j) {
                                tm[i][j] = median(source->b(i, j), source->b(i - 1, j), source->b(i + 1, j), source->b(i, j + 1), source->b(i, j - 1));
                            }
                        } else {
                            for (int j = 1; j < wid - 1; ++j) {
                                tm[i][j] = median(source->b(i, j), source->b(i - 1, j), source->b(i + 1, j), source->b(i, j + 1), source->b(i, j - 1), source->b(i - 1, j - 1), source->b(i - 1, j + 1), source->b(i + 1, j - 1), source->b(i + 1, j + 1));
                            }
                        }
                    }
                } else {
                    #pragma omp for

                    for (int i = 2; i < hei - 2; ++i)
                    {
                        if (methmed == 3) {
                            for (int j = 2; j < wid - 2; ++j) {
                                tm[i][j] = median(source->b(i, j), source->b(i - 1, j), source->b(i + 1, j), source->b(i, j + 1), source->b(i, j - 1), source->b(i - 1, j - 1), source->b(i - 1, j + 1), source->b(i + 1, j - 1), source->b(i + 1, j + 1),
                                source->b(i - 2, j), source->b(i + 2, j), source->b(i, j + 2), source->b(i, j - 2), source->b(i - 2, j - 2), source->b(i - 2, j + 2), source->b(i + 2, j - 2), source->b(i + 2, j + 2),
                                source->b(i - 2, j + 1), source->b(i + 2, j + 1), source->b(i - 1, j + 2), source->b(i - 1, j - 2), source->b(i - 2, j - 1), source->b(i + 2, j - 1), source->b(i + 1, j + 2), source->b(i + 1, j - 2)); // 5x5
                            }
                        } else {
                            for (int j = 2; j < wid - 2; ++j) {
                                tm[i][j] = median(
                                               source->b(i, j),
                                               source->b(i - 1, j),
                                               source->b(i + 1, j),
                                               source->b(i, j + 1),
                                               source->b(i, j - 1),
                                               source->b(i - 1, j - 1),
                                               source->b(i - 1, j + 1),
                                               source->b(i + 1, j - 1),
                                               source->b(i + 1, j + 1),
                                               source->b(i + 2, j),
                                               source->b(i - 2, j),
                                               source->b(i, j + 2),
                                               source->b(i, j - 2)
                                           ); // 5x5 soft
                            }
                        }
                    }
                }

#ifdef _OPENMP
                #pragma omp for nowait
#endif

                for (int i = border; i < hei - border; ++i)
                {
                    for (int j = border; j < wid - border; ++j) {
                        dst->b(i, j) = tm[i][j];
                    }
                }


                if (methmed < 2)
                {
                    #pragma omp for

                    for (int i = 1; i < hei - 1; ++i) {
                        if (methmed == 0) {
                            for (int j = 1; j < wid - 1; ++j) {
                                tm[i][j] = median(source->g(i, j), source->g(i - 1, j), source->g(i + 1, j), source->g(i, j + 1), source->g(i, j - 1));
                            }
                        } else {
                            for (int j = 1; j < wid - 1; ++j) {
                                tm[i][j] = median(source->g(i, j), source->g(i - 1, j), source->g(i + 1, j), source->g(i, j + 1), source->g(i, j - 1), source->g(i - 1, j - 1), source->g(i - 1, j + 1), source->g(i + 1, j - 1), source->g(i + 1, j + 1));
                            }
                        }
                    }
                } else {
                    #pragma omp for

                    for (int i = 2; i < hei - 2; ++i)
                    {
                        if (methmed == 3) {
                            for (int j = 2; j < wid - 2; ++j) {
                                tm[i][j] = median(source->g(i, j), source->g(i - 1, j), source->g(i + 1, j), source->g(i, j + 1), source->g(i, j - 1), source->g(i - 1, j - 1), source->g(i - 1, j + 1), source->g(i + 1, j - 1), source->g(i + 1, j + 1),
                                source->g(i - 2, j), source->g(i + 2, j), source->g(i, j + 2), source->g(i, j - 2), source->g(i - 2, j - 2), source->g(i - 2, j + 2), source->g(i + 2, j - 2), source->g(i + 2, j + 2),
                                source->g(i - 2, j + 1), source->g(i + 2, j + 1), source->g(i - 1, j + 2), source->g(i - 1, j - 2), source->g(i - 2, j - 1), source->g(i + 2, j - 1), source->g(i + 1, j + 2), source->g(i + 1, j - 2)); // 5x5
                            }
                        } else {
                            for (int j = 2; j < wid - 2; ++j) {
                                tm[i][j] = median(
                                               source->g(i, j),
                                               source->g(i - 1, j),
                                               source->g(i + 1, j),
                                               source->g(i, j + 1),
                                               source->g(i, j - 1),
                                               source->g(i - 1, j - 1),
                                               source->g(i - 1, j + 1),
                                               source->g(i + 1, j - 1),
                                               source->g(i + 1, j + 1),
                                               source->g(i + 2, j),
                                               source->g(i - 2, j),
                                               source->g(i, j + 2),
                                               source->g(i, j - 2)
                                           ); // 5x5 soft
                            }
                        }
                    }
                }

#ifdef _OPENMP
                #pragma omp for
#endif

                for (int i = border; i < hei - border; ++i)
                {
                    for (int j = border; j < wid - border; ++j) {
                        dst->g(i, j) = tm[i][j];
                    }
                }
            }
        }

        for (int i = 0; i < hei; ++i) {
            delete[] tm[i];
        }

        delete[] tm;

    }

    //end median
    if (noiseLCurve || useNoiseCCurve) {
        delete[] lumcalcBuffer;
        delete[] lumcalc;
        delete[] ccalcBuffer;
        delete[] ccalc;
    }

//#ifdef _DEBUG
    if (settings->verbose) {
        t2e.set();
        printf("Denoise performed in %d usec:\n", t2e.etime(t1e));
    }

//#endif

}//end of main RGB_denoise


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SSEFUNCTION void ImProcFunctions::RGBtile_denoise (float * fLblox, int hblproc, float noisevar_Ldetail, float * nbrwt, float * blurbuffer) //for DCT
{
    int blkstart = hblproc * TS * TS;

    boxabsblur(fLblox + blkstart, nbrwt, 3, 3, TS, TS, blurbuffer); //blur neighbor weights for more robust estimation //for DCT

#ifdef __SSE2__
    __m128  tempv;
    __m128  noisevar_Ldetailv = _mm_set1_ps(noisevar_Ldetail);
    __m128  onev = _mm_set1_ps(1.0f);

    for (int n = 0; n < TS * TS; n += 4) { //for DCT
        tempv  = onev - xexpf(-SQRV(LVF(nbrwt[n])) / noisevar_Ldetailv);
        _mm_storeu_ps(&fLblox[blkstart + n], LVFU(fLblox[blkstart + n]) * tempv);
    }//output neighbor averaged result

#else

    for (int n = 0; n < TS * TS; ++n) { //for DCT
        fLblox[blkstart + n] *= (1 - xexpf(-SQR(nbrwt[n]) / noisevar_Ldetail));
    }//output neighbor averaged result

#endif

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //printf("vblk=%d  hlk=%d  wsqave=%f   ||   ",vblproc,hblproc,wsqave);

}//end of function tile_denoise


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void ImProcFunctions::RGBoutput_tile_row (float *bloxrow_L, float ** Ldetail, float ** tilemask_out, int height, int width, int top)
{
    const int numblox_W = ceil((static_cast<float>(width)) / (offset));
    const float DCTnorm = 1.0f / (4 * TS * TS); //for DCT

    int imin = MAX(0, -top);
    int bottom = MIN(top + TS, height);
    int imax = bottom - top;

    //add row of tiles to output image
    for (int i = imin; i < imax; ++i) {
        for (int hblk = 0; hblk < numblox_W; ++hblk) {
            int left = (hblk - blkrad) * offset;
            int right  = MIN(left + TS, width);
            int jmin = MAX(0, -left);
            int jmax = right - left;
            int indx = hblk * TS;

            for (int j = jmin; j < jmax; ++j) { // this loop gets auto vectorized by gcc
                Ldetail[top + i][left + j] += tilemask_out[i][j] * bloxrow_L[(indx + i) * TS + j] * DCTnorm; //for DCT

            }
        }
    }
}
/*
#undef TS
#undef fTS
#undef offset
#undef epsilon
*/

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

float ImProcFunctions::MadMax(float * DataList, int & max, int datalen)
{

    //computes Median Absolute Deviation and Maximum of DataList
    //DataList values should mostly have abs val < 65535

    int * histo = new int[65536];

    //memset(histo, 0, 65536*sizeof(histo));
    for (int i = 0; i < 65536; ++i) {
        histo[i] = 0;
    }

    //calculate histogram of absolute values of HH wavelet coeffs
    for (int i = 0; i < datalen; ++i) {
        histo[MAX(0, MIN(65535, abs((int)DataList[i])))]++;
    }

    //find median of histogram
    int median = 0, count = 0;

    while (count < datalen / 2) {
        count += histo[median];
        ++median;
    }

    //find max of histogram
    max = 65535;

    while (histo[max] == 0) {
        max--;
    }

    int count_ = count - histo[median - 1];

    delete[] histo;

    // interpolate
    return (((median - 1) + (datalen / 2 - count_) / (static_cast<float>(count - count_))) / 0.6745);

}

float ImProcFunctions::Mad(float * DataList, const int datalen)
{
    if(datalen <= 1) { // Avoid possible buffer underrun
        return 0;
    }

    //computes Median Absolute Deviation
    //DataList values should mostly have abs val < 256 because we are in Lab mode
    int histo[256] ALIGNED64 = {0};

    //calculate histogram of absolute values of wavelet coeffs
    for (int i = 0; i < datalen; ++i) {
        histo[min(255, abs(static_cast<int>(DataList[i])))]++;
    }

    //find median of histogram
    int median = 0, count = 0;

    while (count < datalen / 2) {
        count += histo[median];
        ++median;
    }

    int count_ = count - histo[median - 1];

    // interpolate
    return (((median - 1) + (datalen / 2 - count_) / (static_cast<float>(count - count_))) / 0.6745);
}

float ImProcFunctions::MadRgb(float * DataList, const int datalen)
{
    if(datalen <= 1) { // Avoid possible buffer underrun
        return 0;
    }

    //computes Median Absolute Deviation
    //DataList values should mostly have abs val < 65536 because we are in RGB mode
    int * histo = new int[65536];

    for (int i = 0; i < 65536; ++i) {
        histo[i] = 0;
    }

    //calculate histogram of absolute values of wavelet coeffs
    int i;

    for (i = 0; i < datalen; ++i) {
        histo[min(65535, abs(static_cast<int>(DataList[i])))]++;
    }

    //find median of histogram
    int median = 0, count = 0;

    while (count < datalen / 2) {
        count += histo[median];
        ++median;
    }

    int count_ = count - histo[median - 1];

    // interpolate
    delete[] histo;
    return (((median - 1) + (datalen / 2 - count_) / (static_cast<float>(count - count_))) / 0.6745);
}



void ImProcFunctions::Noise_residualAB(wavelet_decomposition &WaveletCoeffs_ab, float &chresid, float &chmaxresid, bool denoiseMethodRgb)
{
    int maxlvl = WaveletCoeffs_ab.maxlevel();
    float resid = 0.f;
    float madC;
    float maxresid = 0.f;

    for (int lvl = 0; lvl < maxlvl; ++lvl) {
        // compute median absolute deviation (MAD) of detail coefficients as robust noise estimator

        int Wlvl_ab = WaveletCoeffs_ab.level_W(lvl);
        int Hlvl_ab = WaveletCoeffs_ab.level_H(lvl);

        float ** WavCoeffs_ab = WaveletCoeffs_ab.level_coeffs(lvl);

        for (int dir = 1; dir < 4; ++dir) {
            if (denoiseMethodRgb) {
                madC = SQR(MadRgb(WavCoeffs_ab[dir], Wlvl_ab * Hlvl_ab));
            } else {
                madC = SQR(Mad(WavCoeffs_ab[dir], Wlvl_ab * Hlvl_ab));
            }

            resid += madC;

            if (madC > maxresid) {
                maxresid = madC;
            }
        }
    }

    chresid = resid;
    chmaxresid = maxresid;
}

SSEFUNCTION bool ImProcFunctions::WaveletDenoiseAll_BiShrinkL(wavelet_decomposition &WaveletCoeffs_L, float *noisevarlum, float madL[8][3])
{
    int maxlvl = min(WaveletCoeffs_L.maxlevel(), 5);
    const float eps = 0.01f;

    int maxWL = 0, maxHL = 0;

    for (int lvl = 0; lvl < maxlvl; ++lvl) {
        if (WaveletCoeffs_L.level_W(lvl) > maxWL) {
            maxWL = WaveletCoeffs_L.level_W(lvl);
        }

        if (WaveletCoeffs_L.level_H(lvl) > maxHL) {
            maxHL = WaveletCoeffs_L.level_H(lvl);
        }
    }

    bool memoryAllocationFailed = false;
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif
    {
        float *buffer[3];
        buffer[0] = new (std::nothrow) float[maxWL * maxHL + 32];
        buffer[1] = new (std::nothrow) float[maxWL * maxHL + 64];
        buffer[2] = new (std::nothrow) float[maxWL * maxHL + 96];

        if (buffer[0] == nullptr || buffer[1] == nullptr || buffer[2] == nullptr) {
            memoryAllocationFailed = true;
        }

        if (!memoryAllocationFailed) {

#ifdef _RT_NESTED_OPENMP
            #pragma omp for schedule(dynamic) collapse(2)
#endif

            for (int lvl = maxlvl - 1; lvl >= 0; lvl--) { //for levels less than max, use level diff to make edge mask
                for (int dir = 1; dir < 4; ++dir) {
                    int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
                    int Hlvl_L = WaveletCoeffs_L.level_H(lvl);

                    float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);

                    if (lvl == maxlvl - 1) {
                        int edge = 0;
                        ShrinkAllL(WaveletCoeffs_L, buffer, lvl, dir, noisevarlum, madL[lvl], nullptr, edge);
                    } else {
                        //simple wavelet shrinkage
                        float * sfave = buffer[0] + 32;
                        float * sfaved = buffer[2] + 96;
                        float * blurBuffer = buffer[1] + 64;

                        float mad_Lr = madL[lvl][dir - 1];

                        float levelFactor = mad_Lr * 5.f / (lvl + 1);
#ifdef __SSE2__
                        __m128 mad_Lv;
                        __m128 ninev = _mm_set1_ps(9.0f);
                        __m128 epsv = _mm_set1_ps(eps);
                        __m128 mag_Lv;
                        __m128 levelFactorv = _mm_set1_ps(levelFactor);
                        int coeffloc_L;

                        for (coeffloc_L = 0; coeffloc_L < Hlvl_L * Wlvl_L - 3; coeffloc_L += 4) {
                            mad_Lv = LVFU(noisevarlum[coeffloc_L]) * levelFactorv;
                            mag_Lv = SQRV(LVFU(WavCoeffs_L[dir][coeffloc_L]));
                            _mm_storeu_ps(&sfave[coeffloc_L], mag_Lv / (mag_Lv + mad_Lv * xexpf(-mag_Lv / (mad_Lv * ninev)) + epsv));
                        }

                        for (; coeffloc_L < Hlvl_L * Wlvl_L; ++coeffloc_L) {
                            float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L]);
                            sfave[coeffloc_L] = mag_L / (mag_L + levelFactor * noisevarlum[coeffloc_L] * xexpf(-mag_L / (9.f * levelFactor * noisevarlum[coeffloc_L])) + eps);
                        }

#else

                        for (int i = 0; i < Hlvl_L; ++i) {
                            for (int j = 0; j < Wlvl_L; ++j) {

                                int coeffloc_L = i * Wlvl_L + j;
                                float mag_L = SQR(WavCoeffs_L[dir][coeffloc_L]);
                                sfave[coeffloc_L] = mag_L / (mag_L + levelFactor * noisevarlum[coeffloc_L] * xexpf(-mag_L / (9.f * levelFactor * noisevarlum[coeffloc_L])) + eps);
                            }
                        }

#endif
                        boxblur(sfave, sfaved, blurBuffer, lvl + 2, lvl + 2, Wlvl_L, Hlvl_L); //increase smoothness by locally averaging shrinkage
#ifdef __SSE2__
                        __m128 sfavev;
                        __m128 sf_Lv;

                        for (coeffloc_L = 0; coeffloc_L < Hlvl_L * Wlvl_L - 3; coeffloc_L += 4) {
                            sfavev = LVFU(sfaved[coeffloc_L]);
                            sf_Lv = LVFU(sfave[coeffloc_L]);
                            _mm_storeu_ps(&WavCoeffs_L[dir][coeffloc_L], LVFU(WavCoeffs_L[dir][coeffloc_L]) * (SQRV(sfavev) + SQRV(sf_Lv)) / (sfavev + sf_Lv + epsv));
                            //use smoothed shrinkage unless local shrinkage is much less
                        }

                        // few remaining pixels
                        for (; coeffloc_L < Hlvl_L * Wlvl_L; ++coeffloc_L) {
                            float sf_L = sfave[coeffloc_L];
                            //use smoothed shrinkage unless local shrinkage is much less
                            WavCoeffs_L[dir][coeffloc_L] *= (SQR(sfaved[coeffloc_L]) + SQR(sf_L)) / (sfaved[coeffloc_L] + sf_L + eps);
                        }//now luminance coeffs are denoised

#else

                        for (int i = 0; i < Hlvl_L; ++i) {
                            for (int j = 0; j < Wlvl_L; ++j) {
                                int coeffloc_L = i * Wlvl_L + j;
                                float sf_L = sfave[coeffloc_L];
                                //use smoothed shrinkage unless local shrinkage is much less
                                WavCoeffs_L[dir][coeffloc_L] *= (SQR(sfaved[coeffloc_L]) + SQR(sf_L)) / (sfaved[coeffloc_L] + sf_L + eps);
                            }//now luminance coeffs are denoised
                        }

#endif
                    }
                }
            }
        }

        for (int i = 2; i >= 0; i--) {
            if (buffer[i] != nullptr) {
                delete[] buffer[i];
            }
        }

    }
    return (!memoryAllocationFailed);
}

SSEFUNCTION bool ImProcFunctions::WaveletDenoiseAll_BiShrinkAB(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_ab,
        float *noisevarchrom, float madL[8][3], float noisevar_ab, const bool useNoiseCCurve, bool autoch, bool denoiseMethodRgb)
{
    int maxlvl = WaveletCoeffs_L.maxlevel();

    if (autoch && noisevar_ab <= 0.001f) {
        noisevar_ab = 0.02f;
    }

    float madab[8][3];

    int maxWL = 0, maxHL = 0;

    for (int lvl = 0; lvl < maxlvl; ++lvl) {
        if (WaveletCoeffs_L.level_W(lvl) > maxWL) {
            maxWL = WaveletCoeffs_L.level_W(lvl);
        }

        if (WaveletCoeffs_L.level_H(lvl) > maxHL) {
            maxHL = WaveletCoeffs_L.level_H(lvl);
        }
    }

    bool memoryAllocationFailed = false;
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif
    {
        float *buffer[3];
        buffer[0] = new (std::nothrow) float[maxWL * maxHL + 32];
        buffer[1] = new (std::nothrow) float[maxWL * maxHL + 64];
        buffer[2] = new (std::nothrow) float[maxWL * maxHL + 96];

        if (buffer[0] == nullptr || buffer[1] == nullptr || buffer[2] == nullptr) {
            memoryAllocationFailed = true;
        }

        if (!memoryAllocationFailed) {


#ifdef _RT_NESTED_OPENMP
            #pragma omp for schedule(dynamic) collapse(2)
#endif

            for (int lvl = 0; lvl < maxlvl; ++lvl) {
                for (int dir = 1; dir < 4; ++dir) {
                    // compute median absolute deviation (MAD) of detail coefficients as robust noise estimator
                    int Wlvl_ab = WaveletCoeffs_ab.level_W(lvl);
                    int Hlvl_ab = WaveletCoeffs_ab.level_H(lvl);
                    float ** WavCoeffs_ab = WaveletCoeffs_ab.level_coeffs(lvl);

                    if (!denoiseMethodRgb) {
                        madab[lvl][dir - 1] = SQR(Mad(WavCoeffs_ab[dir], Wlvl_ab * Hlvl_ab));
                    } else {
                        madab[lvl][dir - 1] = SQR(MadRgb(WavCoeffs_ab[dir], Wlvl_ab * Hlvl_ab));
                    }
                }
            }

#ifdef _RT_NESTED_OPENMP
            #pragma omp for schedule(dynamic) collapse(2)
#endif

            for (int lvl = maxlvl - 1; lvl >= 0; lvl--) { //for levels less than max, use level diff to make edge mask
                for (int dir = 1; dir < 4; ++dir) {
                    int Wlvl_ab = WaveletCoeffs_ab.level_W(lvl);
                    int Hlvl_ab = WaveletCoeffs_ab.level_H(lvl);

                    float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
                    float ** WavCoeffs_ab = WaveletCoeffs_ab.level_coeffs(lvl);

                    if (lvl == maxlvl - 1) {
                        ShrinkAllAB(WaveletCoeffs_L, WaveletCoeffs_ab, buffer, lvl, dir, noisevarchrom, noisevar_ab, useNoiseCCurve, autoch, denoiseMethodRgb, madL[lvl], madab[lvl], true);
                    } else {
                        //simple wavelet shrinkage

                        float mad_Lr = madL[lvl][dir - 1];
                        float mad_abr = useNoiseCCurve ? noisevar_ab * madab[lvl][dir - 1] : SQR(noisevar_ab) * madab[lvl][dir - 1];

                        if (noisevar_ab > 0.001f) {

#ifdef __SSE2__
                            __m128 onev = _mm_set1_ps(1.f);
                            __m128 mad_abrv = _mm_set1_ps(mad_abr);
                            __m128 rmad_Lm9v = onev / _mm_set1_ps(mad_Lr * 9.f);
                            __m128 mad_abv;
                            __m128 mag_Lv, mag_abv;
                            __m128 tempabv;
                            int coeffloc_ab;

                            for (coeffloc_ab = 0; coeffloc_ab < Hlvl_ab * Wlvl_ab - 3; coeffloc_ab += 4) {
                                mad_abv = LVFU(noisevarchrom[coeffloc_ab]) * mad_abrv;

                                tempabv = LVFU(WavCoeffs_ab[dir][coeffloc_ab]);
                                mag_Lv = LVFU(WavCoeffs_L[dir][coeffloc_ab]);
                                mag_abv = SQRV(tempabv);
                                mag_Lv = SQRV(mag_Lv) * rmad_Lm9v;
                                _mm_storeu_ps(&WavCoeffs_ab[dir][coeffloc_ab], tempabv * SQRV((onev - xexpf(-(mag_abv / mad_abv) - (mag_Lv)))));
                            }

                            // few remaining pixels
                            for (; coeffloc_ab < Hlvl_ab * Wlvl_ab; ++coeffloc_ab) {
                                float mag_L = SQR(WavCoeffs_L[dir][coeffloc_ab ]);
                                float mag_ab = SQR(WavCoeffs_ab[dir][coeffloc_ab]);
                                WavCoeffs_ab[dir][coeffloc_ab] *= SQR(1.f - xexpf(-(mag_ab / (noisevarchrom[coeffloc_ab] * mad_abr)) - (mag_L / (9.f * mad_Lr)))/*satfactor_a*/);
                            }//now chrominance coefficients are denoised

#else

                            for (int i = 0; i < Hlvl_ab; ++i) {
                                for (int j = 0; j < Wlvl_ab; ++j) {
                                    int coeffloc_ab = i * Wlvl_ab + j;

                                    float mag_L = SQR(WavCoeffs_L[dir][coeffloc_ab ]);
                                    float mag_ab = SQR(WavCoeffs_ab[dir][coeffloc_ab]);

                                    WavCoeffs_ab[dir][coeffloc_ab] *= SQR(1.f - xexpf(-(mag_ab / (noisevarchrom[coeffloc_ab] * mad_abr)) - (mag_L / (9.f * mad_Lr)))/*satfactor_a*/);

                                }
                            }//now chrominance coefficients are denoised

#endif
                        }

                    }
                }
            }
        }

        for (int i = 2; i >= 0; i--) {
            if (buffer[i] != nullptr) {
                delete[] buffer[i];
            }
        }

    }
    return (!memoryAllocationFailed);
}


bool ImProcFunctions::WaveletDenoiseAllL(wavelet_decomposition &WaveletCoeffs_L, float *noisevarlum, float madL[8][3], float * vari, int edge)//mod JD

{

    int maxlvl = min(WaveletCoeffs_L.maxlevel(), 5);

    if (edge == 1) {
        maxlvl = 4;    //for refine denoise edge wavelet
    }

    int maxWL = 0, maxHL = 0;

    for (int lvl = 0; lvl < maxlvl; ++lvl) {
        if (WaveletCoeffs_L.level_W(lvl) > maxWL) {
            maxWL = WaveletCoeffs_L.level_W(lvl);
        }

        if (WaveletCoeffs_L.level_H(lvl) > maxHL) {
            maxHL = WaveletCoeffs_L.level_H(lvl);
        }
    }

    bool memoryAllocationFailed = false;
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif
    {
        float *buffer[4];
        buffer[0] = new (std::nothrow) float[maxWL * maxHL + 32];
        buffer[1] = new (std::nothrow) float[maxWL * maxHL + 64];
        buffer[2] = new (std::nothrow) float[maxWL * maxHL + 96];
        buffer[3] = new (std::nothrow) float[maxWL * maxHL + 128];

        if (buffer[0] == nullptr || buffer[1] == nullptr || buffer[2] == nullptr || buffer[3] == nullptr) {
            memoryAllocationFailed = true;
        }

        if (!memoryAllocationFailed) {
#ifdef _RT_NESTED_OPENMP
            #pragma omp for schedule(dynamic) collapse(2)
#endif

            for (int lvl = 0; lvl < maxlvl; ++lvl) {
                for (int dir = 1; dir < 4; ++dir) {
                    ShrinkAllL(WaveletCoeffs_L, buffer, lvl, dir, noisevarlum, madL[lvl], vari, edge);
                }
            }
        }

        for (int i = 3; i >= 0; i--) {
            if (buffer[i] != nullptr) {
                delete[] buffer[i];
            }
        }
    }
    return (!memoryAllocationFailed);
}


bool ImProcFunctions::WaveletDenoiseAllAB(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_ab,
        float *noisevarchrom, float madL[8][3], float noisevar_ab, const bool useNoiseCCurve, bool autoch, bool denoiseMethodRgb)//mod JD

{

    int maxlvl = WaveletCoeffs_L.maxlevel();
    int maxWL = 0, maxHL = 0;

    for (int lvl = 0; lvl < maxlvl; ++lvl) {
        if (WaveletCoeffs_L.level_W(lvl) > maxWL) {
            maxWL = WaveletCoeffs_L.level_W(lvl);
        }

        if (WaveletCoeffs_L.level_H(lvl) > maxHL) {
            maxHL = WaveletCoeffs_L.level_H(lvl);
        }
    }

    bool memoryAllocationFailed = false;
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel num_threads(denoiseNestedLevels) if (denoiseNestedLevels>1)
#endif
    {
        float *buffer[3];
        buffer[0] = new (std::nothrow) float[maxWL * maxHL + 32];
        buffer[1] = new (std::nothrow) float[maxWL * maxHL + 64];
        buffer[2] = new (std::nothrow) float[maxWL * maxHL + 96];

        if (buffer[0] == nullptr || buffer[1] == nullptr || buffer[2] == nullptr) {
            memoryAllocationFailed = true;
        }

        if (!memoryAllocationFailed) {
#ifdef _RT_NESTED_OPENMP
            #pragma omp for schedule(dynamic) collapse(2)
#endif

            for (int lvl = 0; lvl < maxlvl; ++lvl) {
                for (int dir = 1; dir < 4; ++dir) {
                    ShrinkAllAB(WaveletCoeffs_L, WaveletCoeffs_ab, buffer, lvl, dir, noisevarchrom, noisevar_ab, useNoiseCCurve, autoch, denoiseMethodRgb, madL[lvl]);
                }
            }
        }

        for (int i = 2; i >= 0; i--) {
            if (buffer[i] != nullptr) {
                delete[] buffer[i];
            }
        }
    }
    return (!memoryAllocationFailed);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SSEFUNCTION void ImProcFunctions::ShrinkAllL(wavelet_decomposition &WaveletCoeffs_L, float **buffer, int level, int dir,
        float *noisevarlum, float * madL, float * vari, int edge)

{
    //simple wavelet shrinkage
    const float eps = 0.01f;

    float * sfave = buffer[0] + 32;
    float * sfaved = buffer[1] + 64;
    float * blurBuffer = buffer[2] + 96;

    int W_L = WaveletCoeffs_L.level_W(level);
    int H_L = WaveletCoeffs_L.level_H(level);

    float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(level);
//      printf("OK lev=%d\n",level);
    float mad_L = madL[dir - 1] ;

    if (edge == 1 && vari) {
        noisevarlum = blurBuffer;       // we need one buffer, but fortunately we don't have to allocate a new one because we can use blurBuffer

        for (int i = 0; i < W_L * H_L; ++i) {
            noisevarlum[i] = vari[level];
        }
    }

    float levelFactor = mad_L * 5.f / static_cast<float>(level + 1);
#ifdef __SSE2__
    __m128  magv;
    __m128 levelFactorv = _mm_set1_ps(levelFactor);
    __m128  mad_Lv;
    __m128  ninev = _mm_set1_ps(9.0f);
    __m128  epsv = _mm_set1_ps(eps);
    int i;

    for (i = 0; i < W_L * H_L - 3; i += 4) {
        mad_Lv = LVFU(noisevarlum[i]) * levelFactorv;
        magv = SQRV(LVFU(WavCoeffs_L[dir][i]));
        _mm_storeu_ps(&sfave[i], magv / (magv + mad_Lv * xexpf(-magv / (ninev * mad_Lv)) + epsv));
    }

    // few remaining pixels
    for (; i < W_L * H_L; ++i) {
        float mag = SQR(WavCoeffs_L[dir][i]);
        sfave[i] = mag / (mag + levelFactor * noisevarlum[i] * xexpf(-mag / (9 * levelFactor * noisevarlum[i])) + eps);
    }

#else

    for (int i = 0; i < W_L * H_L; ++i) {

        float mag = SQR(WavCoeffs_L[dir][i]);
        float shrinkfactor = mag / (mag + levelFactor * noisevarlum[i] * xexpf(-mag / (9 * levelFactor * noisevarlum[i])) + eps);
        sfave[i] = shrinkfactor;
    }

#endif
    boxblur(sfave, sfaved, blurBuffer, level + 2, level + 2, W_L, H_L); //increase smoothness by locally averaging shrinkage

#ifdef __SSE2__
    __m128  sfv;

    for (i = 0; i < W_L * H_L - 3; i += 4) {
        sfv = LVFU(sfave[i]);
        //use smoothed shrinkage unless local shrinkage is much less
        _mm_storeu_ps(&WavCoeffs_L[dir][i], _mm_loadu_ps(&WavCoeffs_L[dir][i]) * (SQRV(LVFU(sfaved[i])) + SQRV(sfv)) / (LVFU(sfaved[i]) + sfv + epsv));
    }

    // few remaining pixels
    for (; i < W_L * H_L; ++i) {
        float sf = sfave[i];

        //use smoothed shrinkage unless local shrinkage is much less
        WavCoeffs_L[dir][i] *= (SQR(sfaved[i]) + SQR(sf)) / (sfaved[i] + sf + eps);
    }//now luminance coefficients are denoised

#else

    for (int i = 0; i < W_L * H_L; ++i) {
        float sf = sfave[i];

        //use smoothed shrinkage unless local shrinkage is much less
        WavCoeffs_L[dir][i] *= (SQR(sfaved[i]) + SQR(sf)) / (sfaved[i] + sf + eps);

    }//now luminance coefficients are denoised

#endif
}


SSEFUNCTION void ImProcFunctions::ShrinkAllAB(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_ab, float **buffer, int level, int dir,
        float *noisevarchrom, float noisevar_ab, const bool useNoiseCCurve, bool autoch,
        bool denoiseMethodRgb, float * madL, float * madaab,  bool madCalculated)

{
    //simple wavelet shrinkage
    const float eps = 0.01f;

    if (autoch && noisevar_ab <= 0.001f) {
        noisevar_ab = 0.02f;
    }

    float * sfaveab = buffer[0] + 32;
    float * sfaveabd = buffer[1] + 64;
    float * blurBuffer = buffer[2] + 96;

    int W_ab = WaveletCoeffs_ab.level_W(level);
    int H_ab = WaveletCoeffs_ab.level_H(level);

    float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(level);
    float ** WavCoeffs_ab = WaveletCoeffs_ab.level_coeffs(level);

    float madab;
    float mad_L = madL[dir - 1];

    if (madCalculated) {
        madab = madaab[dir - 1];
    } else {
        if (!denoiseMethodRgb) {
            madab = SQR(Mad(WavCoeffs_ab[dir], W_ab * H_ab));
        } else {
            madab = SQR(MadRgb(WavCoeffs_ab[dir], W_ab * H_ab));
        }
    }

    if (noisevar_ab > 0.001f) {
        madab = useNoiseCCurve ? madab : madab * noisevar_ab;
#ifdef __SSE2__
        __m128 onev = _mm_set1_ps(1.f);
        __m128 mad_abrv = _mm_set1_ps(madab);

        __m128 rmadLm9v = onev / _mm_set1_ps(mad_L * 9.f);
        __m128 mad_abv ;
        __m128 mag_Lv, mag_abv;
        int coeffloc_ab;

        for (coeffloc_ab = 0; coeffloc_ab < H_ab * W_ab - 3; coeffloc_ab += 4) {
            mad_abv = LVFU(noisevarchrom[coeffloc_ab]) * mad_abrv;

            mag_Lv = LVFU(WavCoeffs_L[dir][coeffloc_ab]);
            mag_abv = SQRV(LVFU(WavCoeffs_ab[dir][coeffloc_ab]));
            mag_Lv = (SQRV(mag_Lv)) * rmadLm9v;
            _mm_storeu_ps(&sfaveab[coeffloc_ab], (onev - xexpf(-(mag_abv / mad_abv) - (mag_Lv))));
        }

        // few remaining pixels
        for (; coeffloc_ab < H_ab * W_ab; ++coeffloc_ab) {
            float mag_L = SQR(WavCoeffs_L[dir][coeffloc_ab]);
            float mag_ab = SQR(WavCoeffs_ab[dir][coeffloc_ab]);
            sfaveab[coeffloc_ab] = (1.f - xexpf(-(mag_ab / (noisevarchrom[coeffloc_ab] * madab)) - (mag_L / (9.f * mad_L))));
        }//now chrominance coefficients are denoised

#else

        for (int i = 0; i < H_ab; ++i) {
            for (int j = 0; j < W_ab; ++j) {
                int coeffloc_ab = i * W_ab + j;
                float mag_L = SQR(WavCoeffs_L[dir][coeffloc_ab]);
                float mag_ab = SQR(WavCoeffs_ab[dir][coeffloc_ab]);
                sfaveab[coeffloc_ab] = (1.f - xexpf(-(mag_ab / (noisevarchrom[coeffloc_ab] * madab)) - (mag_L / (9.f * mad_L))));
            }
        }//now chrominance coefficients are denoised

#endif

        boxblur(sfaveab, sfaveabd, blurBuffer, level + 2, level + 2, W_ab, H_ab); //increase smoothness by locally averaging shrinkage
#ifdef __SSE2__
        __m128 epsv = _mm_set1_ps(eps);
        __m128 sfabv;
        __m128 sfaveabv;

        for (coeffloc_ab = 0; coeffloc_ab < H_ab * W_ab - 3; coeffloc_ab += 4) {
            sfabv = LVFU(sfaveab[coeffloc_ab]);
            sfaveabv = LVFU(sfaveabd[coeffloc_ab]);

            //use smoothed shrinkage unless local shrinkage is much less
            _mm_storeu_ps(&WavCoeffs_ab[dir][coeffloc_ab], LVFU(WavCoeffs_ab[dir][coeffloc_ab]) * (SQRV(sfaveabv) + SQRV(sfabv)) / (sfaveabv + sfabv + epsv));
        }

        // few remaining pixels
        for (; coeffloc_ab < H_ab * W_ab; ++coeffloc_ab) {
            //modification Jacques feb 2013
            float sfab = sfaveab[coeffloc_ab];

            //use smoothed shrinkage unless local shrinkage is much less
            WavCoeffs_ab[dir][coeffloc_ab] *= (SQR(sfaveabd[coeffloc_ab]) + SQR(sfab)) / (sfaveabd[coeffloc_ab] + sfab + eps);
        }//now chrominance coefficients are denoised

#else

        for (int i = 0; i < H_ab; ++i) {
            for (int j = 0; j < W_ab; ++j) {
                int coeffloc_ab = i * W_ab + j;
                float sfab = sfaveab[coeffloc_ab];

                //use smoothed shrinkage unless local shrinkage is much less
                WavCoeffs_ab[dir][coeffloc_ab] *= (SQR(sfaveabd[coeffloc_ab]) + SQR(sfab)) / (sfaveabd[coeffloc_ab] + sfab + eps);
            }//now chrominance coefficients are denoised
        }

#endif
    }

}

SSEFUNCTION void ImProcFunctions::ShrinkAll_info(float ** WavCoeffs_a, float ** WavCoeffs_b, int level,
        int W_ab, int H_ab, int skip_ab, float **noisevarlum, float **noisevarchrom, float **noisevarhue, int width, int height, float noisevar_abr, float noisevar_abb, LabImage * noi, float &chaut, int &Nb, float &redaut, float &blueaut,
        float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut,  bool autoch, int schoice, int lvl, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc,
        float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool denoiseMethodRgb, bool multiThread)
{

    //simple wavelet shrinkage
    if (lvl == 1) { //only one time
        float chro = 0.f;
        float dev = 0.f;
        float devL = 0.f;
        int nc = 0;
        int nL = 0;
        int nry = 0;
        float lume = 0.f;
        float red_yel = 0.f;
        float skin_c = 0.f;
        int nsk = 0;

        for (int i = 0; i < H_ab; ++i) {
            for (int j = 0; j < W_ab; ++j) {
                chro += noisevarchrom[i][j];
                ++nc;
                dev += SQR(noisevarchrom[i][j] - (chro / nc));

                if (noisevarhue[i][j] > -0.8f && noisevarhue[i][j] < 2.0f && noisevarchrom[i][j] > 10000.f) {//saturated red yellow
                    red_yel += noisevarchrom[i][j];
                    ++nry;
                }

                if (noisevarhue[i][j] > 0.f && noisevarhue[i][j] < 1.6f && noisevarchrom[i][j] < 10000.f) {//skin
                    skin_c += noisevarchrom[i][j];
                    ++nsk;
                }

                lume += noisevarlum[i][j];
                ++nL;
                devL += SQR(noisevarlum[i][j] - (lume / nL));
            }
        }

        if (nc > 0) {
            chromina = chro / nc;
            sigma = sqrt(dev / nc);
            nsknc = static_cast<float>(nsk) / static_cast<float>(nc);
        } else {
            nsknc = static_cast<float>(nsk);
        }

        if (nL > 0) {
            lumema = lume / nL;
            sigma_L = sqrt(devL / nL);
        }

        if (nry > 0) {
            redyel = red_yel / nry;
        }

        if (nsk > 0) {
            skinc = skin_c / nsk;
        }
    }

    const float reduc = (schoice == 2) ? static_cast<float>( settings->nrhigh) : 1.f;

    for (int dir = 1; dir < 4; ++dir) {
        float mada, madb;

        if (!denoiseMethodRgb) {
            mada = SQR(Mad(WavCoeffs_a[dir], W_ab * H_ab));
        } else {
            mada = SQR(MadRgb(WavCoeffs_a[dir], W_ab * H_ab));
        }

        chred += mada;

        if (mada > maxchred) {
            maxchred = mada;
        }

        if (mada < minchred) {
            minchred = mada;
        }

        maxredaut = sqrt(reduc * maxchred);
        minredaut = sqrt(reduc * minchred);

        if (!denoiseMethodRgb) {
            madb = SQR(Mad(WavCoeffs_b[dir], W_ab * H_ab));
        } else {
            madb = SQR(MadRgb(WavCoeffs_b[dir], W_ab * H_ab));
        }

        chblue += madb;

        if (madb > maxchblue) {
            maxchblue = madb;
        }

        if (madb < minchblue) {
            minchblue = madb;
        }

        maxblueaut = sqrt(reduc * maxchblue);
        minblueaut = sqrt(reduc * minchblue);

        chau += (mada + madb);
        ++nb;
        //here evaluation of automatic
        chaut = sqrt(reduc * chau / (nb + nb));
        redaut = sqrt(reduc * chred / nb);
        blueaut = sqrt(reduc * chblue / nb);
        Nb = nb;
    }

}


void ImProcFunctions::WaveletDenoiseAll_info(int levwav, wavelet_decomposition &WaveletCoeffs_a,
        wavelet_decomposition &WaveletCoeffs_b, float **noisevarlum, float **noisevarchrom, float **noisevarhue, int width, int height, float noisevar_abr, float noisevar_abb, LabImage * noi, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut, int schoice, bool autoch,
        float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc, float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool denoiseMethodRgb, bool multiThread)
{

    int maxlvl = levwav;

    for (int lvl = 0; lvl < maxlvl; ++lvl) {

        int Wlvl_ab = WaveletCoeffs_a.level_W(lvl);
        int Hlvl_ab = WaveletCoeffs_a.level_H(lvl);

        int skip_ab = WaveletCoeffs_a.level_stride(lvl);

        float ** WavCoeffs_a = WaveletCoeffs_a.level_coeffs(lvl);
        float ** WavCoeffs_b = WaveletCoeffs_b.level_coeffs(lvl);

        ShrinkAll_info(WavCoeffs_a, WavCoeffs_b, lvl, Wlvl_ab, Hlvl_ab,
                       skip_ab, noisevarlum,  noisevarchrom, noisevarhue, width, height, noisevar_abr, noisevar_abb, noi, chaut, Nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut,
                       autoch, schoice, lvl, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc, maxchred, maxchblue, minchred, minchblue, nb, chau, chred, chblue, denoiseMethodRgb, multiThread);

    }
}

SSEFUNCTION void ImProcFunctions::RGB_denoise_infoGamCurve(const procparams::DirPyrDenoiseParams & dnparams, bool isRAW, LUTf &gamcurve, float &gam, float &gamthresh, float &gamslope)
{
    gam = dnparams.gamma;
    gamthresh = 0.001f;

    if (!isRAW) {//reduce gamma under 1 for Lab mode ==> TIF and JPG
        if (gam < 1.9f) {
            gam = 1.f - (1.9f - gam) / 3.f;    //minimum gamma 0.7
        } else if (gam >= 1.9f && gam <= 3.f) {
            gam = (1.4f / 1.1f) * gam - 1.41818f;
        }
    }

    bool denoiseMethodRgb = (dnparams.dmethod == "RGB");

    if (denoiseMethodRgb) {
        gamslope = exp(log(static_cast<double>(gamthresh)) / gam) / gamthresh;
        Color::gammaf2lut(gamcurve, gam, gamthresh, gamslope, 65535.f, 32768.f);
    } else {
        Color::gammanf2lut(gamcurve, gam, 65535.f, 32768.f);
    }
}

void ImProcFunctions::calcautodn_info (float &chaut, float &delta, int Nb, int levaut, float maxmax, float lumema, float chromina, int mode, int lissage, float redyel, float skinc, float nsknc)
{

    float reducdelta = 1.f;

    if (params->dirpyrDenoise.smethod == "shalbi") {
        reducdelta = static_cast<float>(settings->nrhigh);
    }

    chaut = (chaut * Nb - maxmax) / (Nb - 1); //suppress maximum for chaut calcul

    if ((redyel > 5000.f || skinc > 1000.f) && nsknc < 0.4f  && chromina > 3000.f) {
        chaut *= 0.45f;    //reduct action in red zone, except skin for high / med chroma
    } else if ((redyel > 12000.f || skinc > 1200.f) && nsknc < 0.3f && chromina > 3000.f) {
        chaut *= 0.3f;
    }

    if (mode == 0 || mode == 2) { //Preview or Auto multizone
        if (chromina > 10000.f) {
            chaut *= 0.7f;    //decrease action for high chroma  (visible noise)
        } else if (chromina > 6000.f) {
            chaut *= 0.9f;
        } else if (chromina < 3000.f) {
            chaut *= 1.2f;    //increase action in low chroma==> 1.2  /==>2.0 ==> curve CC
        } else if (chromina < 2000.f) {
            chaut *= 1.5f;    //increase action in low chroma==> 1.5 / ==>2.7
        }

        if (lumema < 2500.f) {
            chaut *= 1.3f;    //increase action for low light
        } else if (lumema < 5000.f) {
            chaut *= 1.2f;
        } else if (lumema > 20000.f) {
            chaut *= 0.9f;    //decrease for high light
        }
    } else if (mode == 1) {//auto ==> less coefficient because interaction
        if (chromina > 10000.f) {
            chaut *= 0.8f;    //decrease action for high chroma  (visible noise)
        } else if (chromina > 6000.f) {
            chaut *= 0.9f;
        } else if (chromina < 3000.f) {
            chaut *= 1.5f;    //increase action in low chroma
        } else if (chromina < 2000.f) {
            chaut *= 2.2f;    //increase action in low chroma
        }

        if (lumema < 2500.f) {
            chaut *= 1.2f;    //increase action for low light
        } else if (lumema < 5000.f) {
            chaut *= 1.1f;
        } else if (lumema > 20000.f) {
            chaut *= 0.9f;    //decrease for high light
        }
    }

    if (levaut == 0) { //Low denoise
        if (chaut > 300.f) {
            chaut = 0.714286f * chaut + 85.71428f;
        }
    }

    delta = maxmax - chaut;
    delta *= reducdelta;

    if (lissage == 1 || lissage == 2) {
        if (chaut < 200.f && delta < 200.f) {
            delta *= 0.95f;
        } else if (chaut < 200.f && delta < 400.f) {
            delta *= 0.5f;
        } else if (chaut < 200.f && delta >= 400.f) {
            delta = 200.f;
        } else if (chaut < 400.f && delta < 400.f) {
            delta *= 0.4f;
        } else if (chaut < 400.f && delta >= 400.f) {
            delta = 120.f;
        } else if (chaut < 550.f) {
            delta *= 0.15f;
        } else if (chaut < 650.f) {
            delta *= 0.1f;
        } else /*if (chaut >= 650.f)*/ {
            delta *= 0.07f;
        }

        if (mode == 0 || mode == 2) { //Preview or Auto multizone
            if (chromina < 6000.f) {
                delta *= 1.4f;    //increase maxi
            }

            if (lumema < 5000.f) {
                delta *= 1.4f;
            }
        } else if (mode == 1) { //Auto
            if (chromina < 6000.f) {
                delta *= 1.2f;    //increase maxi
            }

            if (lumema < 5000.f) {
                delta *= 1.2f;
            }
        }
    }

    if (lissage == 0) {
        if (chaut < 200.f && delta < 200.f) {
            delta *= 0.95f;
        } else if (chaut < 200.f && delta < 400.f) {
            delta *= 0.7f;
        } else if (chaut < 200.f && delta >= 400.f) {
            delta = 280.f;
        } else if (chaut < 400.f && delta < 400.f) {
            delta *= 0.6f;
        } else if (chaut < 400.f && delta >= 400.f) {
            delta = 200.f;
        } else if (chaut < 550.f) {
            delta *= 0.3f;
        } else if (chaut < 650.f) {
            delta *= 0.2f;
        } else /*if (chaut >= 650.f)*/ {
            delta *= 0.15f;
        }

        if (mode == 0 || mode == 2) { //Preview or Auto multizone
            if (chromina < 6000.f) {
                delta *= 1.4f;    //increase maxi
            }

            if (lumema < 5000.f) {
                delta *= 1.4f;
            }
        } else if (mode == 1) { //Auto
            if (chromina < 6000.f) {
                delta *= 1.2f;    //increase maxi
            }

            if (lumema < 5000.f) {
                delta *= 1.2f;
            }
        }
    }

}

SSEFUNCTION void ImProcFunctions::RGB_denoise_info(Imagefloat * src, Imagefloat * provicalc, const bool isRAW, LUTf &gamcurve, float gam, float gamthresh, float gamslope, const procparams::DirPyrDenoiseParams & dnparams, const double expcomp, float &chaut, int &Nb,  float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc, bool multiThread)
{
    if ((settings->leveldnautsimpl == 1 && dnparams.Cmethod == "MAN") || (settings->leveldnautsimpl == 0 && dnparams.C2method == "MANU")) {
        //nothing to do
        return;
    }

    int hei, wid;
    float** lumcalc;
    float** acalc;
    float** bcalc;
    hei = provicalc->getHeight();
    wid = provicalc->getWidth();
    TMatrix wprofi = iccStore->workingSpaceMatrix (params->icm.working);

    const float wpi[3][3] = {
        {static_cast<float>(wprofi[0][0]), static_cast<float>(wprofi[0][1]), static_cast<float>(wprofi[0][2])},
        {static_cast<float>(wprofi[1][0]), static_cast<float>(wprofi[1][1]), static_cast<float>(wprofi[1][2])},
        {static_cast<float>(wprofi[2][0]), static_cast<float>(wprofi[2][1]), static_cast<float>(wprofi[2][2])}
    };

    lumcalc = new float*[hei];

    for (int i = 0; i < hei; ++i) {
        lumcalc[i] = new float[wid];
    }

    acalc = new float*[hei];

    for (int i = 0; i < hei; ++i) {
        acalc[i] = new float[wid];
    }

    bcalc = new float*[hei];

    for (int i = 0; i < hei; ++i) {
        bcalc[i] = new float[wid];
    }

#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for if (multiThread)
#endif

    for (int ii = 0; ii < hei; ++ii) {
        for (int jj = 0; jj < wid; ++jj) {
            float LLum, AAum, BBum;
            float RL = provicalc->r(ii, jj);
            float GL = provicalc->g(ii, jj);
            float BL = provicalc->b(ii, jj);
            // determine luminance for noisecurve
            float XL, YL, ZL;
            Color::rgbxyz(RL, GL, BL, XL, YL, ZL, wpi);
            Color::XYZ2Lab(XL, YL, ZL, LLum, AAum, BBum);
            lumcalc[ii][jj] = LLum;
            acalc[ii][jj] = AAum;
            bcalc[ii][jj] = BBum;
        }
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    const int imheight = src->getHeight(), imwidth = src->getWidth();

    bool denoiseMethodRgb = (dnparams.dmethod == "RGB");

    const float gain = pow (2.0f, float(expcomp));

    int tilesize;
    int overlap;

    if (settings->leveldnti == 0) {
        tilesize = 1024;
        overlap = 128;
    }

    if (settings->leveldnti == 1) {
        tilesize = 768;
        overlap = 96;
    }

    int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;

    //always no Tiles
    int kall = 0;
    Tile_calc (tilesize, overlap, kall, imwidth, imheight, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TMatrix wprof = iccStore->workingSpaceMatrix (params->icm.working);
    const float wp[3][3] = {
        {static_cast<float>(wprof[0][0]), static_cast<float>(wprof[0][1]), static_cast<float>(wprof[0][2])},
        {static_cast<float>(wprof[1][0]), static_cast<float>(wprof[1][1]), static_cast<float>(wprof[1][2])},
        {static_cast<float>(wprof[2][0]), static_cast<float>(wprof[2][1]), static_cast<float>(wprof[2][2])}
    };

    float chau = 0.f;
    float chred = 0.f;
    float chblue = 0.f;
    float maxchred = 0.f;
    float maxchblue = 0.f;
    float minchred = 100000000.f;
    float minchblue = 100000000.f;
    int nb = 0;
    int comptlevel = 0;

    for (int tiletop = 0; tiletop < imheight; tiletop += tileHskip) {
        for (int tileleft = 0; tileleft < imwidth; tileleft += tileWskip) {

            int tileright = MIN(imwidth, tileleft + tilewidth);
            int tilebottom = MIN(imheight, tiletop + tileheight);
            int width  = tileright - tileleft;
            int height = tilebottom - tiletop;
            LabImage * labdn = new LabImage(width, height);
            float** noisevarlum = new float*[(height + 1) / 2];

            for (int i = 0; i < (height + 1) / 2; ++i) {
                noisevarlum[i] = new float[(width + 1) / 2];
            }

            float** noisevarchrom = new float*[(height + 1) / 2];

            for (int i = 0; i < (height + 1) / 2; ++i) {
                noisevarchrom[i] = new float[(width + 1) / 2];
            }

            float** noisevarhue = new float*[(height + 1) / 2];

            for (int i = 0; i < (height + 1) / 2; ++i) {
                noisevarhue[i] = new float[(width + 1) / 2];
            }

            // init luma noisevarL
            float noisevarab_b, noisevarab_r;

            float realred, realblue;
            float interm_med = static_cast<float>( dnparams.chroma) / 10.0;
            float intermred, intermblue;

            if (dnparams.redchro > 0.) {
                intermred = (dnparams.redchro / 10.);
            } else {
                intermred = static_cast<float>( dnparams.redchro) / 7.0;    //increase slower than linear for more sensit
            }

            if (dnparams.bluechro > 0.) {
                intermblue = (dnparams.bluechro / 10.);
            } else {
                intermblue = static_cast<float>( dnparams.bluechro) / 7.0;    //increase slower than linear for more sensit
            }

            realred = interm_med + intermred;

            if (realred < 0.f) {
                realred = 0.001f;
            }

            realblue = interm_med + intermblue;

            if (realblue < 0.f) {
                realblue = 0.001f;
            }

            //TODO: implement using AlignedBufferMP
            //fill tile from image; convert RGB to "luma/chroma"

            if (isRAW) {//image is raw; use channel differences for chroma channels
#ifdef _RT_NESTED_OPENMP
                #pragma omp parallel for if (multiThread)
#endif

                for (int i = tiletop; i < tilebottom; i += 2) {
                    int i1 = i - tiletop;
#ifdef __SSE2__
                    __m128 aNv, bNv;
                    __m128 c100v = _mm_set1_ps(100.f);
                    int j;

                    for (j = tileleft; j < tileright - 7; j += 8) {
                        int j1 = j - tileleft;
                        aNv = LVFU(acalc[i >> 1][j >> 1]);
                        bNv = LVFU(bcalc[i >> 1][j >> 1]);
                        _mm_storeu_ps(&noisevarhue[i1 >> 1][j1 >> 1], xatan2f(bNv, aNv));
                        _mm_storeu_ps(&noisevarchrom[i1 >> 1][j1 >> 1], _mm_max_ps(c100v, _mm_sqrt_ps(SQRV(aNv) + SQRV(bNv))));
                    }

                    for (; j < tileright; j += 2) {
                        int j1 = j - tileleft;
                        float aN = acalc[i >> 1][j >> 1];
                        float bN = bcalc[i >> 1][j >> 1];
                        float cN = sqrtf(SQR(aN) + SQR(bN));
                        noisevarhue[i1 >> 1][j1 >> 1] = xatan2f(bN, aN);

                        if (cN < 100.f) {
                            cN = 100.f;    //avoid divided by zero
                        }

                        noisevarchrom[i1 >> 1][j1 >> 1] = cN;
                    }

#else

                    for (int j = tileleft; j < tileright; j += 2) {
                        int j1 = j - tileleft;
                        float aN = acalc[i >> 1][j >> 1];
                        float bN = bcalc[i >> 1][j >> 1];
                        float cN = sqrtf(SQR(aN) + SQR(bN));
                        float hN = xatan2f(bN, aN);

                        if (cN < 100.f) {
                            cN = 100.f;    //avoid divided by zero
                        }

                        noisevarchrom[i1 >> 1][j1 >> 1] = cN;
                        noisevarhue[i1 >> 1][j1 >> 1] = hN;
                    }

#endif
                }

#ifdef _RT_NESTED_OPENMP
                #pragma omp parallel for if (multiThread)
#endif

                for (int i = tiletop; i < tilebottom; i += 2) {
                    int i1 = i - tiletop;

                    for (int j = tileleft; j < tileright; j += 2) {
                        int j1 = j - tileleft;
                        float Llum = lumcalc[i >> 1][j >> 1];
                        Llum = Llum < 2.f ? 2.f : Llum; //avoid divided by zero ?
                        Llum = Llum > 32768.f ? 32768.f : Llum; // not strictly necessary
                        noisevarlum[i1 >> 1][j1 >> 1] = Llum;
                    }
                }

                if (!denoiseMethodRgb) { //lab mode, modification Jacques feb 2013 and july 2014

#ifdef _RT_NESTED_OPENMP
                    #pragma omp parallel for if (multiThread)
#endif

                    for (int i = tiletop; i < tilebottom; ++i) {
                        int i1 = i - tiletop;

                        for (int j = tileleft; j < tileright; ++j) {
                            int j1 = j - tileleft;
                            float R_ = gain * src->r(i, j);
                            float G_ = gain * src->g(i, j);
                            float B_ = gain * src->b(i, j);

                            R_ = Color::denoiseIGammaTab[R_];
                            G_ = Color::denoiseIGammaTab[G_];
                            B_ = Color::denoiseIGammaTab[B_];

                            //apply gamma noise standard (slider)
                            R_ = R_ < 65535.f ? gamcurve[R_] : (Color::gammanf(R_ / 65535.f, gam) * 32768.f);
                            G_ = G_ < 65535.f ? gamcurve[G_] : (Color::gammanf(G_ / 65535.f, gam) * 32768.f);
                            B_ = B_ < 65535.f ? gamcurve[B_] : (Color::gammanf(B_ / 65535.f, gam) * 32768.f);
                            //true conversion xyz=>Lab
                            float X, Y, Z;
                            Color::rgbxyz(R_, G_, B_, X, Y, Z, wp);

                            //convert to Lab
                            float L, a, b;
                            Color::XYZ2Lab(X, Y, Z, L, a, b);

                            labdn->a[i1][j1] = a;
                            labdn->b[i1][j1] = b;
                        }
                    }
                } else { //RGB mode

                    for (int i = tiletop/*, i1=0*/; i < tilebottom; ++i/*, ++i1*/) {
                        int i1 = i - tiletop;

                        for (int j = tileleft/*, j1=0*/; j < tileright; ++j/*, ++j1*/) {
                            int j1 = j - tileleft;

                            float X = gain * src->r(i, j);
                            float Y = gain * src->g(i, j);
                            float Z = gain * src->b(i, j);

                            X = X < 65535.f ? gamcurve[X] : (Color::gammaf(X / 65535.f, gam, gamthresh, gamslope) * 32768.f);
                            Y = Y < 65535.f ? gamcurve[Y] : (Color::gammaf(Y / 65535.f, gam, gamthresh, gamslope) * 32768.f);
                            Z = Z < 65535.f ? gamcurve[Z] : (Color::gammaf(Z / 65535.f, gam, gamthresh, gamslope) * 32768.f);

                            labdn->a[i1][j1] = (X - Y);
                            labdn->b[i1][j1] = (Y - Z);
                        }
                    }
                }

            } else {//image is not raw; use Lab parametrization
                for (int i = tiletop/*, i1=0*/; i < tilebottom; ++i/*, ++i1*/) {
                    int i1 = i - tiletop;

                    for (int j = tileleft/*, j1=0*/; j < tileright; ++j/*, ++j1*/) {
                        int j1 = j - tileleft;
                        float L, a, b;
                        float rLum = src->r(i, j) ; //for luminance denoise curve
                        float gLum = src->g(i, j) ;
                        float bLum = src->b(i, j) ;

                        //use gamma sRGB, not good if TIF (JPG) Output profil not with gamma sRGB  (eg : gamma =1.0, or 1.8...)
                        //very difficult to solve !
                        // solution ==> save TIF with gamma sRGB and re open
                        float rtmp = Color::igammatab_srgb[ src->r(i, j) ];
                        float gtmp = Color::igammatab_srgb[ src->g(i, j) ];
                        float btmp = Color::igammatab_srgb[ src->b(i, j) ];
                        //modification Jacques feb 2013
                        // gamma slider different from raw
                        rtmp = rtmp < 65535.f ? gamcurve[rtmp] : (Color::gammanf(rtmp / 65535.f, gam) * 32768.f);
                        gtmp = gtmp < 65535.f ? gamcurve[gtmp] : (Color::gammanf(gtmp / 65535.f, gam) * 32768.f);
                        btmp = btmp < 65535.f ? gamcurve[btmp] : (Color::gammanf(btmp / 65535.f, gam) * 32768.f);

                        float X, Y, Z;
                        Color::rgbxyz(rtmp, gtmp, btmp, X, Y, Z, wp);

                        //convert Lab
                        Color::XYZ2Lab(X, Y, Z, L, a, b);

                        if (((i1 | j1) & 1) == 0) {
                            float Llum, alum, blum;
                            float XL, YL, ZL;
                            Color::rgbxyz(rLum, gLum, bLum, XL, YL, ZL, wp);
                            Color::XYZ2Lab(XL, YL, ZL, Llum, alum, blum);
                            float kN = Llum;

                            if (kN < 2.f) {
                                kN = 2.f;
                            }

                            if (kN > 32768.f) {
                                kN = 32768.f;
                            }

                            noisevarlum[i1 >> 1][j1 >> 1] = kN;
                            float aN = alum;
                            float bN = blum;
                            float hN = xatan2f(bN, aN);
                            float cN = sqrt(SQR(aN) + SQR(bN));

                            if (cN < 100.f) {
                                cN = 100.f;    //avoid divided by zero
                            }

                            noisevarchrom[i1 >> 1][j1 >> 1] = cN;
                            noisevarhue[i1 >> 1][j1 >> 1] = hN;
                        }

                        labdn->a[i1][j1] = a;
                        labdn->b[i1][j1] = b;
                    }
                }
            }

            int datalen = labdn->W * labdn->H;

            //now perform basic wavelet denoise
            //last two arguments of wavelet decomposition are max number of wavelet decomposition levels;
            //and whether to subsample the image after wavelet filtering.  Subsampling is coded as
            //binary 1 or 0 for each level, eg subsampling = 0 means no subsampling, 1 means subsample
            //the first level only, 7 means subsample the first three levels, etc.
            noisevarab_r = SQR(realred) + 0.01f;
            noisevarab_b = SQR(realblue) + 0.01f;

            wavelet_decomposition* adecomp;
            wavelet_decomposition* bdecomp;

            int schoice = 0;//shrink method

            if (dnparams.smethod == "shalbi") {
                schoice = 2;
            }

            const int levwav = 5;
#ifdef _RT_NESTED_OPENMP
            #pragma omp parallel sections if (multiThread)
#endif
            {
#ifdef _RT_NESTED_OPENMP
                #pragma omp section
#endif
                {
                    adecomp = new wavelet_decomposition (labdn->data + datalen, labdn->W, labdn->H, levwav, 1);
                }
#ifdef _RT_NESTED_OPENMP
                #pragma omp section
#endif
                {
                    bdecomp = new wavelet_decomposition (labdn->data + 2 * datalen, labdn->W, labdn->H, levwav, 1);
                }
            }
            bool autoch = dnparams.autochroma;

            if (comptlevel == 0) {
                WaveletDenoiseAll_info(
                    levwav,
                    *adecomp,
                    *bdecomp,
                    noisevarlum,
                    noisevarchrom,
                    noisevarhue,
                    width,
                    height,
                    noisevarab_r,
                    noisevarab_b,
                    labdn,
                    chaut,
                    Nb,
                    redaut,
                    blueaut,
                    maxredaut,
                    maxblueaut,
                    minredaut,
                    minblueaut,
                    schoice,
                    autoch,
                    chromina,
                    sigma,
                    lumema,
                    sigma_L,
                    redyel,
                    skinc,
                    nsknc,
                    maxchred,
                    maxchblue,
                    minchred,
                    minchblue,
                    nb,
                    chau,
                    chred,
                    chblue,
                    denoiseMethodRgb,
                    multiThread
                ); // Enhance mode
            }

            comptlevel += 1;
            delete adecomp;
            delete bdecomp;
            delete labdn;

            for (int i = 0; i < (height + 1) / 2; ++i) {
                delete[] noisevarlum[i];
            }

            delete[] noisevarlum;

            for (int i = 0; i < (height + 1) / 2; ++i) {
                delete[] noisevarchrom[i];
            }

            delete[] noisevarchrom;

            for (int i = 0; i < (height + 1) / 2; ++i) {
                delete[] noisevarhue[i];
            }

            delete[] noisevarhue;

        }//end of tile row
    }//end of tile loop

    for (int i = 0; i < hei; ++i) {
        delete[] lumcalc[i];
    }

    delete[] lumcalc;

    for (int i = 0; i < hei; ++i) {
        delete[] acalc[i];
    }

    delete[] acalc;

    for (int i = 0; i < hei; ++i) {
        delete[] bcalc[i];
    }

    delete[] bcalc;

#undef TS
#undef fTS
#undef offset
#undef epsilon

} // End of main RGB_denoise

}
