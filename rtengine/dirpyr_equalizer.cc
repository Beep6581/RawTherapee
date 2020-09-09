/*
 *  This file is part of RawTherapee.
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
 *
 *  (C) 2010 Emil Martinec <ejmartin@uchicago.edu>
 *
*/

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "array2D.h"
#include "cieimage.h"
#include "color.h"
#include "curves.h"
#include "improcfun.h"
#include "LUT.h"
#include "opthelper.h"
#include "rt_math.h"
#include "settings.h"

namespace {

float rangeFn(float i) {
    return 1.f / (i + 1000.f);
}


void dirpyr_channel(const float * const * data_fine, float ** data_coarse, int width, int height, int level, int scale)
{
    // scale is spacing of directional averaging weights
    // calculate weights, compute directionally weighted average

    if (level > 1) {
        //generate domain kernel
        //  multiplied each value of domker by 1000 to avoid multiplication by 1000 inside the loop
#ifdef __SSE2__
        const float domkerv[5][5][4] ALIGNED16 = {{{1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}},
                                                  {{1000, 1000, 1000, 1000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {1000, 1000, 1000, 1000}},
                                                  {{1000, 1000, 1000, 1000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {1000, 1000, 1000, 1000}},
                                                  {{1000, 1000, 1000, 1000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {1000, 1000, 1000, 1000}},
                                                  {{1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}}};
#endif
        const float domker[5][5] = {{1000, 1000, 1000, 1000, 1000},
                                    {1000, 2000, 2000, 2000, 1000},
                                    {1000, 2000, 2000, 2000, 1000},
                                    {1000, 2000, 2000, 2000, 1000},
                                    {1000, 1000, 1000, 1000, 1000}};
        constexpr int halfwin = 2;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            const int scalewin = halfwin * scale;
#ifdef __SSE2__
            const vfloat thousandv = F2V(1000.f);
#endif

#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < height; i++) {
                int j;
                for (j = 0; j < scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for (int inbr = std::max(0, i - scalewin); inbr <= std::min(height - 1, i + scalewin); inbr += scale) {
                        for (int jnbr = std::max(0, j - scalewin); jnbr <= j + scalewin; jnbr += scale) {
                            const float dirwt = domker[(inbr - i) / scale + halfwin][(jnbr - j)/ scale + halfwin] * rangeFn(fabsf(data_fine[inbr][jnbr] - data_fine[i][j])); 
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }

#ifdef __SSE2__

                for (; j < width - scalewin - 3; j += 4) {
                    vfloat valv = ZEROV;
                    vfloat normv = ZEROV;
                    const vfloat dftemp1v = LVFU(data_fine[i][j]);

                    for (int inbr = MAX(0, i - scalewin); inbr <= MIN(height - 1, i + scalewin); inbr += scale) {
                        const int indexihlp = (inbr - i) / scale + halfwin;
                        for (int jnbr = j - scalewin, indexjhlp = 0; jnbr <= j + scalewin; jnbr += scale, ++indexjhlp) {
                            const vfloat dftemp2v = LVFU(data_fine[inbr][jnbr]);
                            const vfloat dirwtv = LVF(domkerv[indexihlp][indexjhlp]) / (vabsf(dftemp1v - dftemp2v) + thousandv);
                            valv += dirwtv * dftemp2v;
                            normv += dirwtv;
                        }
                    }
                    STVFU(data_coarse[i][j], valv / normv); //low pass filter
                }
#endif
                for (; j < width - scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for (int inbr = std::max(0, i - scalewin); inbr <= std::min(height - 1, i + scalewin); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= j + scalewin; jnbr += scale) {
                            const float dirwt = domker[(inbr - i) / scale + halfwin][(jnbr - j)/ scale + halfwin] * rangeFn(fabsf(data_fine[inbr][jnbr] - data_fine[i][j])); 
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }
                    data_coarse[i][j] = val / norm; //low pass filter
                }

                for (; j < width; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for (int inbr = std::max(0, i - scalewin); inbr <= std::min(height - 1, i + scalewin); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= std::min(width - 1, j + scalewin); jnbr += scale) {
                            const float dirwt = domker[(inbr - i) / scale + halfwin][(jnbr - j)/ scale + halfwin] * rangeFn(fabsf(data_fine[inbr][jnbr] - data_fine[i][j])); 
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }
                    data_coarse[i][j] = val / norm; //low pass filter
                }
            }
        }
    } else {    // level <=1 means that all values of domker would be 1.0f, so no need for multiplication
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            const vfloat thousandv = F2V(1000.0f);
#endif
#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int i = 0; i < height; i++)
            {
                int j = 0;
                for (; j < scale; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for (int inbr = std::max(0, i - scale); inbr <= std::min(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = std::max(0, j - scale); jnbr <= j + scale; jnbr += scale) {
                            const float dirwt = rangeFn(fabsf(data_fine[inbr][jnbr] - data_fine[i][j]));
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }
                    data_coarse[i][j] = val / norm; //low pass filter
                }

#ifdef __SSE2__

                for (; j < width - scale - 3; j += 4) {
                    vfloat valv = ZEROV;
                    vfloat normv = ZEROV;
                    const vfloat dftemp1v = LVFU(data_fine[i][j]);

                    for (int inbr = MAX(0, i - scale); inbr <= MIN(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = j - scale; jnbr <= j + scale; jnbr += scale) {
                            const vfloat dftemp2v = LVFU(data_fine[inbr][jnbr]);
                            const vfloat dirwtv = thousandv / (vabsf(dftemp2v - dftemp1v) + thousandv);
                            valv += dirwtv * dftemp2v;
                            normv += dirwtv;
                        }
                    }
                    STVFU(data_coarse[i][j], valv / normv); //low pass filter
                }
#endif

                for (; j < width - scale; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for (int inbr = std::max(0, i - scale); inbr <= std::min(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = j - scale; jnbr <= j + scale; jnbr += scale) {
                            const float dirwt = rangeFn(fabsf(data_fine[inbr][jnbr] - data_fine[i][j]));
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }
                    data_coarse[i][j] = val / norm; //low pass filter
                }

                for (; j < width; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for (int inbr = std::max(0, i - scale); inbr <= std::min(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = j - scale; jnbr <= std::min(width - 1, j + scale); jnbr += scale) {
                            const float dirwt = rangeFn(fabsf(data_fine[inbr][jnbr] - data_fine[i][j]));
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }
                    data_coarse[i][j] = val / norm; //low pass filter
                }
            }
        }
    }
}

void fillLut(LUTf &irangefn, int level, double dirpyrThreshold, float mult, float skinprot) {

    float multbis;
    if (level == 4 && mult > 1.f) {
        multbis = 1.f + 0.65f * (mult - 1.f);
    } else if (level == 5 && mult > 1.f) {
        multbis = 1.f + 0.45f * (mult - 1.f);
    } else {
        multbis = mult; //multbis to reduce artifacts for high values mult
    }

    const float offs = skinprot == 0.f ? 0.f : -1.f;
    constexpr double noise = 2000.0;
    const float noisehi = 1.33 * noise * dirpyrThreshold / exp(level * log(3.0)), noiselo = 0.66 * noise * dirpyrThreshold / exp(level * log(3.0));

    for (int i = 0; i < 0x20000; i++) {
        if (abs(i - 0x10000) > noisehi || multbis < 1.f) {
            irangefn[i] = multbis + offs;
        } else {
            if (abs(i - 0x10000) < noiselo) {
                irangefn[i] = 1.f + offs;
            } else {
                irangefn[i] = 1.f + offs + (multbis - 1.f) * (noisehi - abs(i - 0x10000)) / (noisehi - noiselo + 0.01f);
            }
        }
    }
}

void idirpyr_eq_channel_loc(float ** data_coarse, float ** data_fine, float ** buffer, int width, int height, int level, float mult, const double dirpyrThreshold, float ** hue, float ** chrom, const double skinprot, const bool gamutlab, float b_l, float t_l, float t_r, float b_r, int choice, int scaleprev, bool multiThread)
{
    LUTf irangefn(0x20000);
    fillLut(irangefn, level, dirpyrThreshold, mult, skinprot);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            const float hipass = data_fine[i][j] - data_coarse[i][j];
            buffer[i][j] += irangefn[hipass + 0x10000] * hipass;
        }
    }
}

void idirpyr_eq_channel(const float * const * data_coarse, const float * const * data_fine, float ** buffer, int width, int height, int level, float mult, const double dirpyrThreshold, const float * const * hue, const float * const * chrom, const double skinprot, float b_l, float t_l, float t_r)
{
    const float skinprotneg = -skinprot;
    const float factorHard = (1.f - skinprotneg / 100.f);

    LUTf irangefn(0x20000);
    fillLut(irangefn, level, dirpyrThreshold, mult, skinprot);

    if (!skinprot) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                const float hipass = data_fine[i][j] - data_coarse[i][j];
                buffer[i][j] += irangefn[hipass + 0x10000] * hipass;
            }
        }
    } else if (skinprot > 0.0) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                float scale = 1.f;
                const float hipass = data_fine[i][j] - data_coarse[i][j];
                rtengine::Color::SkinSatCbdl(data_fine[i][j] / 327.68f, hue[i][j], chrom[i][j], skinprot, scale, true, b_l, t_l, t_r);
                buffer[i][j] += (1.f + (irangefn[hipass + 0x10000]) * scale) * hipass;
            }
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                float scale = 1.f;
                const float hipass = data_fine[i][j] - data_coarse[i][j];
                rtengine::Color::SkinSatCbdl(data_fine[i][j] / 327.68f, hue[i][j], chrom[i][j], skinprotneg, scale, false, b_l, t_l, t_r);
                const float correct = irangefn[hipass + 0x10000];

                if (scale == 1.f) {//image hard
                    buffer[i][j] += (1.f + correct * factorHard) * hipass;
                } else { //image soft with scale < 1 ==> skin
                    buffer[i][j] += (1.f + correct) * hipass;
                }
            }
        }
    }
}

void idirpyr_eq_channelcam(const float * const * data_coarse, const float * const * data_fine, float ** buffer, int width, int height, int level, float mult, const double dirpyrThreshold, const float * const * h_p, const float * const * C_p, const double skinprot, float b_l, float t_l, float t_r)
{

    const float skinprotneg = -skinprot;
    const float factorHard = 1.f - skinprotneg / 100.f;

    LUTf irangefn(0x20000);
    fillLut(irangefn, level, dirpyrThreshold, mult, skinprot);

    if (!skinprot) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                const float hipass = data_fine[i][j] - data_coarse[i][j];
                buffer[i][j] += irangefn[hipass + 0x10000] * hipass;
            }
        }
    } else if (skinprot > 0.0) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                const float hipass = data_fine[i][j] - data_coarse[i][j];
                float scale = 1.f;
                rtengine::Color::SkinSatCbdlCam(data_fine[i][j] / 327.68f, h_p[i][j] , C_p[i][j], skinprot, scale, true, b_l, t_l, t_r);
                buffer[i][j] += (1.f + (irangefn[hipass + 0x10000]) * scale) * hipass;
            }
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                const float hipass = data_fine[i][j] - data_coarse[i][j];
                float scale = 1.f;
                const float correct = irangefn[hipass + 0x10000];
                rtengine::Color::SkinSatCbdlCam(data_fine[i][j] / 327.68f, h_p[i][j], C_p[i][j], skinprotneg, scale, false, b_l, t_l, t_r);

                if (scale == 1.f) {//image hard
                    buffer[i][j] += (1.f + correct * factorHard) * hipass;
                } else { //image soft
                    buffer[i][j] += (1.f + correct) * hipass;
                }
            }
        }
    }
}

}

namespace rtengine
{

void ImProcFunctions::dirpyr_equalizer(const float * const * src, float ** dst, int srcwidth, int srcheight, const float * const * l_a, const float * const * l_b, const double * mult, const double dirpyrThreshold, const double skinprot, float b_l, float t_l, float t_r, int scaleprev)
{
    //sequence of scales
    constexpr int maxlevel = 6;
    constexpr int scales[maxlevel] = {1, 2, 4, 8, 16, 32};
    const float atten123 = rtengine::LIM<float>(settings->level123_cbdl, 0.f, 50.f);
    const float atten0 = rtengine::LIM<float>(settings->level0_cbdl, 0.f, 40.f);

    int lastlevel = maxlevel;
    while (lastlevel > 0 && fabs(mult[lastlevel - 1] - 1) < 0.001) {
        --lastlevel;
    }

    if (lastlevel == 0) {
        return;
    }

    float multi[maxlevel];

    for (int lv = 0; lv < maxlevel; ++lv) {
        if (scales[lv] < scaleprev) {
            const float factor = lv >= 1 ? atten123 : atten0;
            multi[lv] = (factor * ((float) mult[lv] - 1.f) / 100.f) + 1.f;    //modulate action if zoom < 100%
        } else {
            multi[lv] = mult[lv];
        }
    }

    multi_array2D<float, maxlevel> dirpyrlo (srcwidth, srcheight);

    dirpyr_channel(src, dirpyrlo[0], srcwidth, srcheight, 0, std::max(scales[0] / scaleprev, 1));

    for (int level = 1; level < lastlevel; ++level) {
        dirpyr_channel(dirpyrlo[level - 1], dirpyrlo[level], srcwidth, srcheight, level, std::max(scales[level] / scaleprev, 1));
    }

    array2D<float> tmpHue, tmpChr;

    if (skinprot) {
        // precalculate hue and chroma, use SSE, if available
        // by precalculating these values we can greatly reduce the number of calculations in idirpyr_eq_channel()
        // but we need two additional buffers for this preprocessing
        tmpHue(srcwidth, srcheight);
        tmpChr(srcwidth, srcheight);

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            const vfloat div = F2V(327.68f);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = 0; i < srcheight; i++) {
                int j = 0;
#ifdef __SSE2__
                for (; j < srcwidth - 3; j += 4) {
                    const vfloat lav = LVFU(l_a[i][j]);
                    const vfloat lbv = LVFU(l_b[i][j]);
                    STVFU(tmpHue[i][j], xatan2f(lbv, lav));
                    STVFU(tmpChr[i][j], vsqrtf(SQRV(lbv) + SQRV(lav)) / div);
                }
#endif
                for (; j < srcwidth; j++) {
                    tmpHue[i][j] = xatan2f(l_b[i][j], l_a[i][j]);
                    tmpChr[i][j] = sqrtf(SQR((l_b[i][j])) + SQR((l_a[i][j]))) / 327.68f;
                }
            }
        }
    }

    // with the current implementation of idirpyr_eq_channel we can safely use the buffer from last level as buffer, saves some memory
    float** buffer = dirpyrlo[lastlevel - 1];

    for (int level = lastlevel - 1; level > 0; --level) {
        idirpyr_eq_channel(dirpyrlo[level], dirpyrlo[level - 1], buffer, srcwidth, srcheight, level, multi[level], dirpyrThreshold, tmpHue, tmpChr, skinprot, b_l, t_l, t_r);
    }

    idirpyr_eq_channel(dirpyrlo[0], dst, buffer, srcwidth, srcheight, 0, multi[0], dirpyrThreshold, tmpHue, tmpChr, skinprot, b_l, t_l, t_r);

#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 0; i < srcheight; i++) {
        for (int j = 0; j < srcwidth; j++) {
            dst[i][j] = buffer[i][j];
        }
    }
}

void ImProcFunctions::dirpyr_equalizercam(const CieImage *ncie, float ** src, float ** dst, int srcwidth, int srcheight, const float * const * h_p, const float * const * C_p, const double * mult, const double dirpyrThreshold, const double skinprot, float b_l, float t_l, float t_r, int scaleprev)
{

    //sequence of scales
    constexpr int maxlevel = 6;
    constexpr int scales[maxlevel] = {1, 2, 4, 8, 16, 32};
    const float atten123 = rtengine::LIM<float>(settings->level123_cbdl, 0.f, 50.f);
    const float atten0 = rtengine::LIM<float>(settings->level0_cbdl, 0.f, 40.f);

    int lastlevel = maxlevel;
    while (fabs(mult[lastlevel - 1] - 1) < 0.001 && lastlevel > 0) {
        --lastlevel;
    }

    if (lastlevel == 0) {
        return;
    }

    float multi[maxlevel];

    for (int lv = 0; lv < maxlevel; lv++) {
        if (scales[lv] < scaleprev) {
            const float factor = lv >= 1 ? atten123 : atten0;
            multi[lv] = (factor * ((float) mult[lv] - 1.f) / 100.f) + 1.f;
        } else {
            multi[lv] = mult[lv];
        }
    }

    multi_array2D<float, maxlevel> dirpyrlo (srcwidth, srcheight);

    dirpyr_channel(src, dirpyrlo[0], srcwidth, srcheight, 0, std::max(scales[0] / scaleprev, 1));

    for (int level = 1; level < lastlevel; ++level) {
        dirpyr_channel(dirpyrlo[level - 1], dirpyrlo[level], srcwidth, srcheight, level, std::max(scales[level] / scaleprev, 1));
    }

    // with the current implementation of idirpyr_eq_channel we can safely use the buffer from last level as buffer, saves some memory
    float ** buffer = dirpyrlo[lastlevel - 1];

    for (int level = lastlevel - 1; level > 0; --level) {
        idirpyr_eq_channelcam(dirpyrlo[level], dirpyrlo[level - 1], buffer, srcwidth, srcheight, level, multi[level], dirpyrThreshold , h_p, C_p, skinprot, b_l, t_l, t_r);
    }

    idirpyr_eq_channelcam(dirpyrlo[0], dst, buffer, srcwidth, srcheight, 0, multi[0], dirpyrThreshold, h_p, C_p, skinprot, b_l, t_l, t_r);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int i = 0; i < srcheight; i++) {
        for (int j = 0; j < srcwidth; j++) {
            if (ncie->J_p[i][j] > 8.f && ncie->J_p[i][j] < 92.f) {
                dst[i][j] = buffer[i][j];
            } else {
                dst[i][j] = src[i][j];
            }
        }
    }
}

void ImProcFunctions::cbdl_local_temp(float ** src, float ** loctemp, int srcwidth, int srcheight, const float * mult, float kchro, const double dirpyrThreshold, const float mergeL, const float contres, const double skinprot, const bool gamutlab, float b_l, float t_l, float t_r, float b_r, int choice, int scaleprev, bool multiThread)
{
    constexpr int maxlevelloc = 6;
    constexpr int scalesloc[maxlevelloc] = {1, 2, 4, 8, 16, 32};
    const float atten123 = rtengine::LIM<float>(settings->level123_cbdl, 0.f, 50.f);
    const float atten0 = rtengine::LIM<float>(settings->level0_cbdl, 0.f, 40.f);
    int lastlevel = maxlevelloc;

    if (settings->verbose) { 
        printf("Dirpyr scaleprev=%i\n", scaleprev);
    }

    while (lastlevel > 0 && fabs(mult[lastlevel - 1] - 1) < 0.001) {

        lastlevel--;
        //printf("last level to process %d \n",lastlevel);
    }

    if (lastlevel == 0) {
        return;
    }

    float multi[6];

    for (int lv = 0; lv < 6; ++lv) {
        if (scalesloc[lv] < scaleprev) {
            const float factor = lv >= 1 ? atten123 : atten0;
            multi[lv] = (factor * ((float) mult[lv] - 1.f) / 100.f) + 1.f;    //modulate action if zoom < 100%
        } else {
            multi[lv] = mult[lv];
        }
    }

    if (settings->verbose) {
        printf("CbDL local mult0=%f  1=%f 2=%f 3=%f 4=%f 5%f\n", multi[0], multi[1], multi[2], multi[3], multi[4], multi[5]);
    }

    multi_array2D<float, maxlevelloc> dirpyrlo(srcwidth, srcheight);


    dirpyr_channel(src, dirpyrlo[0], srcwidth, srcheight, 0, std::max(scalesloc[0] / scaleprev, 1));


    for (int level = 1; level < lastlevel; ++level) {
        dirpyr_channel(dirpyrlo[level - 1], dirpyrlo[level], srcwidth, srcheight, level, std::max(scalesloc[level] / scaleprev, 1));
    }

    // with the current implementation of idirpyr_eq_channel we can safely use the buffer from last level as buffer, saves some memory
//    float ** buffer = dirpyrlo[lastlevel - 1];
        array2D<float> residbuff(srcwidth, srcheight);
        array2D<float> resid5(srcwidth, srcheight);
        
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < srcheight; i++)
            for (int j = 0; j < srcwidth; j++) {
                residbuff[i][j] = 0.f;
            }

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < srcheight; i++)
            for (int j = 0; j < srcwidth; j++) {
                residbuff[i][j] = dirpyrlo[lastlevel - 1][i][j];
                resid5[i][j] = dirpyrlo[lastlevel - 1][i][j];               
            }
   
    
    double avg = 0.f;
    if(contres != 0.f) {
        int ng = 0;
        
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:avg, ng)
#endif
        for (int i = 0; i < srcheight; i++) {
            for (int j = 0; j < srcwidth; j++) {
                avg += residbuff[i][j];
                ng++;
            }
        }
        avg /= ng;
        avg /= 32768.f;
        avg = LIM01(avg);
    }
    float contreal = 0.3f * contres;
    DiagonalCurve resid_contrast({
        DCT_NURBS,
            0, 0,
            avg - avg * (0.6 - contreal / 250.0), avg - avg * (0.6 + contreal / 250.0),
            avg + (1 - avg) * (0.6 - contreal / 250.0), avg + (1 - avg) * (0.6 + contreal / 250.0),
            1, 1
        });

    if(contres != 0.f) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < srcheight; i++)
            for (int j = 0; j < srcwidth; j++) {
                float buf = LIM01(residbuff[i][j] / 32768.f);
                buf = resid_contrast.getVal(buf);
                buf *= 32768.f;
                residbuff[i][j] = buf;
            }
    }
    
    
    for (int level = lastlevel - 1; level > 0; level--) {
        idirpyr_eq_channel_loc(dirpyrlo[level], dirpyrlo[level - 1], residbuff, srcwidth, srcheight, level, multi[level], dirpyrThreshold, nullptr, nullptr, skinprot, gamutlab, b_l, t_l, t_r, b_r, choice, scaleprev, multiThread);
    }

    scale = scalesloc[0];

    idirpyr_eq_channel_loc(dirpyrlo[0], src, residbuff, srcwidth, srcheight, 0, multi[0], dirpyrThreshold, nullptr, nullptr, skinprot, gamutlab, b_l, t_l, t_r, b_r, choice, scaleprev, multiThread);

    array2D<float> loct(srcwidth, srcheight);
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < srcheight; i++) {
        for (int j = 0; j < srcwidth; j++) {
            loct[i][j] = LIM(residbuff[i][j],0.f,32768.f);  // TODO: Really a clip necessary?
        }
    }

    float clar = 0.01f * mergeL;

/*
    if(clar == 0.f) {
        clar = 0.0f;
    }
//    printf("clar=%f \n", clar);
*/
    if(clar > 0.f) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < srcheight; i++) {
            for (int j = 0; j < srcwidth; j++) {
                 loctemp[i][j] = LIM((1.f + clar) * loct[i][j] - clar * resid5[i][j],0.f,32768.f);
            }
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < srcheight; i++) {
            for (int j = 0; j < srcwidth; j++) {
                loctemp[i][j] = LIM(loct[i][j],0.f,32768.f);
            }
        }
    }
}

}
