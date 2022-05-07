/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017-2018 Ingo Weyrich <heckflosse67@gmx.de>
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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "gauss.h"
#include "opthelper.h"
#include "rt_algo.h"
#include "rt_math.h"
#include "sleef.h"

namespace {

float calcBlendFactor(float val, float threshold) {
    // sigmoid function
    // result is in ]0;1] range
    // inflexion point is at (x, y) (threshold, 0.5)
    const float x = -16.f + (16.f / threshold) * val;
    return 0.5f * (1.f + x / std::sqrt(1.f + rtengine::SQR(x)));
}

#ifdef __SSE2__
vfloat calcBlendFactor(vfloat val, vfloat threshold) {
    // sigmoid function
    // result is in ]0;1] range
    // inflexion point is at (x, y) (threshold, 0.5)
    const vfloat x = -16.f + (16.f / threshold) * val;
    return 0.5f * (1.f + x * _mm_rsqrt_ps(1.f + rtengine::SQR(x)));
}
#endif

float tileAverage(const float * const *data, size_t tileY, size_t tileX, size_t tilesize) {

    float avg = 0.f;
#ifdef __SSE2__
    vfloat avgv = ZEROV;
#endif
    for (std::size_t y = tileY; y < tileY + tilesize; ++y) {
        std::size_t x = tileX;
#ifdef __SSE2__
        for (; x < tileX + tilesize - 3; x += 4) {
            avgv += LVFU(data[y][x]);
        }
#endif
        for (; x < tileX + tilesize; ++x) {
            avg += data[y][x];
        }
    }
#ifdef __SSE2__
    avg += vhadd(avgv);
#endif
    return avg / rtengine::SQR(tilesize);
}

float tileVariance(const float * const *data, size_t tileY, size_t tileX, size_t tilesize, float avg) {

    float var = 0.f;
#ifdef __SSE2__
    vfloat varv = ZEROV;
    const vfloat avgv = F2V(avg);
#endif
    for (std::size_t y = tileY; y < tileY + tilesize; ++y) {
        std::size_t x = tileX;
#ifdef __SSE2__
        for (; x < tileX + tilesize - 3; x += 4) {
            varv += SQRV(LVFU(data[y][x]) - avgv);
        }
#endif
        for (; x < tileX + tilesize; ++x) {
            var += rtengine::SQR(data[y][x] - avg);
        }
    }
#ifdef __SSE2__
    var += vhadd(varv);
#endif
    return var / (rtengine::SQR(tilesize) * avg);
}

float calcContrastThreshold(const float* const * luminance, int tileY, int tileX, int tilesize) {

    constexpr float scale = 0.0625f / 327.68f;
    std::vector<std::vector<float>> blend(tilesize - 4, std::vector<float>(tilesize - 4));

#ifdef __SSE2__
    const vfloat scalev = F2V(scale);
#endif

    for(int j = tileY + 2; j < tileY + tilesize - 2; ++j) {
        int i = tileX + 2;
#ifdef __SSE2__
        for(; i < tileX + tilesize - 5; i += 4) {
            vfloat contrastv = vsqrtf(SQRV(LVFU(luminance[j][i+1]) - LVFU(luminance[j][i-1])) + SQRV(LVFU(luminance[j+1][i]) - LVFU(luminance[j-1][i])) +
                                      SQRV(LVFU(luminance[j][i+2]) - LVFU(luminance[j][i-2])) + SQRV(LVFU(luminance[j+2][i]) - LVFU(luminance[j-2][i]))) * scalev;
            STVFU(blend[j - tileY - 2][i - tileX - 2], contrastv);
        }
#endif
        for(; i < tileX + tilesize - 2; ++i) {

            float contrast = sqrtf(rtengine::SQR(luminance[j][i+1] - luminance[j][i-1]) + rtengine::SQR(luminance[j+1][i] - luminance[j-1][i]) + 
                                   rtengine::SQR(luminance[j][i+2] - luminance[j][i-2]) + rtengine::SQR(luminance[j+2][i] - luminance[j-2][i])) * scale;

            blend[j - tileY - 2][i - tileX - 2] = contrast;
        }
    }

    const float limit = rtengine::SQR(tilesize - 4) / 100.f;

    int c;
    for (c = 1; c < 100; ++c) {
        const float contrastThreshold = c / 100.f;
        float sum = 0.f;
#ifdef __SSE2__
        const vfloat contrastThresholdv = F2V(contrastThreshold);
        vfloat sumv = ZEROV;
#endif

        for(int j = 0; j < tilesize - 4; ++j) {
            int i = 0;
#ifdef __SSE2__
            for(; i < tilesize - 7; i += 4) {
                sumv += calcBlendFactor(LVFU(blend[j][i]), contrastThresholdv);
            }
#endif
            for(; i < tilesize - 4; ++i) {
                sum += calcBlendFactor(blend[j][i], contrastThreshold);
            }
        }
#ifdef __SSE2__
        sum += vhadd(sumv);
#endif
        if (sum <= limit) {
            break;
        }
    }

    return (c + 1) / 100.f;
}
}

namespace rtengine
{

void findMinMaxPercentile(const float* data, size_t size, float minPrct, float& minOut, float maxPrct, float& maxOut, bool multithread)
{
    // Copyright (c) 2017 Ingo Weyrich <heckflosse67@gmx.de>
    // We need to find the (minPrct*size) smallest value and the (maxPrct*size) smallest value in data.
    // We use a histogram based search for speed and to reduce memory usage.
    // Memory usage of this method is histoSize * sizeof(uint32_t) * (t + 1) byte,
    // where t is the number of threads and histoSize is in [1;65536].
    // Processing time is O(n) where n is size of the input array.
    // It scales well with multiple threads if the size of the input array is large.
    // The current implementation is not guaranteed to work correctly if size > 2^32 (4294967296).

    assert(minPrct <= maxPrct);

    if (size == 0) {
        return;
    }

    size_t numThreads = 1;
#ifdef _OPENMP
    // Because we have an overhead in the critical region of the main loop for each thread
    // we make a rough calculation to reduce the number of threads for small data size.
    // This also works fine for the minmax loop.
    if (multithread) {
        const size_t maxThreads = omp_get_max_threads();
        while (size > numThreads * numThreads * 16384 && numThreads < maxThreads) {
            ++numThreads;
        }
    }
#endif

    // We need min and max value of data to calculate the scale factor for the histogram
    float minVal = data[0];
    float maxVal = data[0];
#ifdef _OPENMP
    #pragma omp parallel for reduction(min:minVal) reduction(max:maxVal) num_threads(numThreads)
#endif
    for (size_t i = 1; i < size; ++i) {
        minVal = std::min(minVal, data[i]);
        maxVal = std::max(maxVal, data[i]);
    }

    if (std::fabs(maxVal - minVal) == 0.f) { // fast exit, also avoids division by zero in calculation of scale factor
        minOut = maxOut = minVal;
        return;
    }

    // Caution: Currently this works correctly only for histoSize in range[1;65536].
    // For small data size (i.e. thumbnails) we reduce the size of the histogram to the size of data.
    const unsigned int histoSize = std::min<size_t>(65536, size);

    // calculate scale factor to use full range of histogram
    const float scale = (histoSize - 1) / (maxVal - minVal);

    // We need one main histogram
    std::vector<uint32_t> histo(histoSize, 0);

    if (numThreads == 1) {
        // just one thread => use main histogram
        for (size_t i = 0; i < size; ++i) {
            // we have to subtract minVal and multiply with scale to get the data in [0;histosize] range
            histo[static_cast<uint16_t>(scale * (data[i] - minVal))]++;
        }
    } else {
#ifdef _OPENMP
    #pragma omp parallel num_threads(numThreads)
#endif
        {
            // We need one histogram per thread
            std::vector<uint32_t> histothr(histoSize, 0);

#ifdef _OPENMP
            #pragma omp for nowait
#endif
            for (size_t i = 0; i < size; ++i) {
                // we have to subtract minVal and multiply with scale to get the data in [0;histosize] range
                histothr[static_cast<uint16_t>(scale * (data[i] - minVal))]++;
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                // add per thread histogram to main histogram
#ifdef _OPENMP
                #pragma omp simd
#endif

                for (size_t i = 0; i < histoSize; ++i) {
                    histo[i] += histothr[i];
                }
            }
        }
    }

    size_t k = 0;
    size_t count = 0;

    // find (minPrct*size) smallest value
    const float threshmin = minPrct * size;
    while (count < threshmin) {
        count += histo[k++];
    }

    if (k > 0) { // interpolate
        const size_t count_ = count - histo[k - 1];
        const float c0 = count - threshmin;
        const float c1 = threshmin - count_;
        minOut = (c1 * k + c0 * (k - 1)) / (c0 + c1);
    } else {
        minOut = k;
    }
    // go back to original range
    minOut /= scale;
    minOut += minVal;
    minOut = rtengine::LIM(minOut, minVal, maxVal);

    // find (maxPrct*size) smallest value
    const float threshmax = maxPrct * size;
    while (count < threshmax) {
        count += histo[k++];
    }

    if (k > 0) { // interpolate
        const size_t count_ = count - histo[k - 1];
        const float c0 = count - threshmax;
        const float c1 = threshmax - count_;
        maxOut = (c1 * k + c0 * (k - 1)) / (c0 + c1);
    } else {
        maxOut = k;
    }
    // go back to original range
    maxOut /= scale;
    maxOut += minVal;
    maxOut = rtengine::LIM(maxOut, minVal, maxVal);
}

void buildBlendMask(const float* const * luminance, float **blend, int W, int H, float &contrastThreshold, bool autoContrast, float ** clipMask) {

    if (autoContrast) {
        constexpr float minLuminance = 2000.f;
        constexpr float maxLuminance = 20000.f;
        constexpr float minTileVariance = 0.5f;
        for (int pass = 0; pass < 2; ++pass) {
            const int tilesize = 80 / (pass + 1);
            const int skip = pass == 0 ? tilesize : tilesize / 4;
            const int numTilesW = W / skip - 3 * pass;
            const int numTilesH = H / skip - 3 * pass;
            std::vector<std::vector<float>> variances(numTilesH, std::vector<float>(numTilesW));

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic)
#endif
            for (int i = 0; i < numTilesH; ++i) {
                const int tileY = i * skip;
                for (int j = 0; j < numTilesW; ++j) {
                    const int tileX = j * skip;
                    const float avg = tileAverage(luminance, tileY, tileX, tilesize);
                    if (avg < minLuminance || avg > maxLuminance) {
                        // too dark or too bright => skip the tile
                        variances[i][j] = RT_INFINITY_F;
                        continue;
                    } else {
                        variances[i][j] = tileVariance(luminance, tileY, tileX, tilesize, avg);
                        // exclude tiles with a variance less than minTileVariance
                        variances[i][j] = variances[i][j] < minTileVariance ? RT_INFINITY_F : variances[i][j];
                    }
                }
            }

            float minvar = RT_INFINITY_F;
            int minI = 0, minJ = 0;
            for (int i = 0; i < numTilesH; ++i) {
                for (int j = 0; j < numTilesW; ++j) {
                    if (variances[i][j] < minvar) {
                        minvar = variances[i][j];
                        minI = i;
                        minJ = j;
                    }
                }
            }

            if (minvar <= 1.f || pass == 1) {
                const int minY = skip * minI;
                const int minX = skip * minJ;
                if (pass == 0) {
                    // a variance <= 1 means we already found a flat region and can skip second pass
                    contrastThreshold = calcContrastThreshold(luminance, minY, minX, tilesize);
                    break;
                } else {
                    // in second pass we allow a variance of 8
                    // we additionally scan the tiles +-skip pixels around the best tile from pass 2
                    // Means we scan (2 * skip + 1)^2 tiles in this step to get a better hit rate
                    // fortunately the scan is quite fast, so we use only one core and don't parallelize
                    const int topLeftYStart = std::max(minY - skip, 0);
                    const int topLeftXStart = std::max(minX - skip, 0);
                    const int topLeftYEnd = std::min(minY + skip, H - tilesize);
                    const int topLeftXEnd = std::min(minX + skip, W - tilesize);
                    const int numTilesH = topLeftYEnd - topLeftYStart + 1;
                    const int numTilesW = topLeftXEnd - topLeftXStart + 1;

                    std::vector<std::vector<float>> variances(numTilesH, std::vector<float>(numTilesW));
                    for (int i = 0; i < numTilesH; ++i) {
                        const int tileY = topLeftYStart + i;
                        for (int j = 0; j < numTilesW; ++j) {
                            const int tileX = topLeftXStart + j;
                            const float avg = tileAverage(luminance, tileY, tileX, tilesize);

                            if (avg < minLuminance || avg > maxLuminance) {
                                // too dark or too bright => skip the tile
                                variances[i][j] = RT_INFINITY_F;
                                continue;
                            } else {
                                variances[i][j] = tileVariance(luminance, tileY, tileX, tilesize, avg);
                            // exclude tiles with a variance less than minTileVariance
                            variances[i][j] = variances[i][j] < minTileVariance ? RT_INFINITY_F : variances[i][j];
                            }
                        }
                    }

                    float minvar = RT_INFINITY_F;
                    int minI = 0, minJ = 0;
                    for (int i = 0; i < numTilesH; ++i) {
                        for (int j = 0; j < numTilesW; ++j) {
                            if (variances[i][j] < minvar) {
                                minvar = variances[i][j];
                                minI = i;
                                minJ = j;
                            }
                        }
                    }

                    contrastThreshold = minvar <= 8.f ? calcContrastThreshold(luminance, topLeftYStart + minI, topLeftXStart + minJ, tilesize) : 0.f;
                }
            }
        }
    }

    if(contrastThreshold == 0.f) {
        for(int j = 0; j < H; ++j) {
            for(int i = 0; i < W; ++i) {
                blend[j][i] = 1.f;
            }
        }
    } else {
        constexpr float scale = 0.0625f / 327.68f;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            const vfloat contrastThresholdv = F2V(contrastThreshold);
            const vfloat scalev = F2V(scale);
#endif
#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for(int j = 2; j < H - 2; ++j) {
                int i = 2;
#ifdef __SSE2__
                if (clipMask) {
                    for(; i < W - 5; i += 4) {
                        vfloat contrastv = vsqrtf(SQRV(LVFU(luminance[j][i+1]) - LVFU(luminance[j][i-1])) + SQRV(LVFU(luminance[j+1][i]) - LVFU(luminance[j-1][i])) +
                                                  SQRV(LVFU(luminance[j][i+2]) - LVFU(luminance[j][i-2])) + SQRV(LVFU(luminance[j+2][i]) - LVFU(luminance[j-2][i]))) * scalev;

                        STVFU(blend[j][i], LVFU(clipMask[j][i]) * calcBlendFactor(contrastv, contrastThresholdv));
                    }
                } else {
                    for(; i < W - 5; i += 4) {
                        vfloat contrastv = vsqrtf(SQRV(LVFU(luminance[j][i+1]) - LVFU(luminance[j][i-1])) + SQRV(LVFU(luminance[j+1][i]) - LVFU(luminance[j-1][i])) +
                                                  SQRV(LVFU(luminance[j][i+2]) - LVFU(luminance[j][i-2])) + SQRV(LVFU(luminance[j+2][i]) - LVFU(luminance[j-2][i]))) * scalev;

                        STVFU(blend[j][i], calcBlendFactor(contrastv, contrastThresholdv));
                    }
                }
#endif
                for(; i < W - 2; ++i) {

                    float contrast = sqrtf(rtengine::SQR(luminance[j][i+1] - luminance[j][i-1]) + rtengine::SQR(luminance[j+1][i] - luminance[j-1][i]) + 
                                           rtengine::SQR(luminance[j][i+2] - luminance[j][i-2]) + rtengine::SQR(luminance[j+2][i] - luminance[j-2][i])) * scale;

                    blend[j][i] = (clipMask ? clipMask[j][i] : 1.f) * calcBlendFactor(contrast, contrastThreshold);
                }
            }

#ifdef _OPENMP
            #pragma omp single
#endif
            {
                // upper border
                for(int j = 0; j < 2; ++j) {
                    for(int i = 2; i < W - 2; ++i) {
                        blend[j][i] = blend[2][i];
                    }
                }
                // lower border
                for(int j = H - 2; j < H; ++j) {
                    for(int i = 2; i < W - 2; ++i) {
                        blend[j][i] = blend[H-3][i];
                    }
                }
                for(int j = 0; j < H; ++j) {
                    // left border
                    blend[j][0] = blend[j][1] = blend[j][2];
                    // right border
                    blend[j][W - 2] = blend[j][W - 1] = blend[j][W - 3];
                }
            }

#ifdef __SSE2__
            // flush denormals to zero for gaussian blur to avoid performance penalty if there are a lot of zero values in the mask
            const auto oldMode = _MM_GET_FLUSH_ZERO_MODE();
            _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

            // blur blend mask to smooth transitions
            gaussianBlur(blend, blend, W, H, 2.0);

#ifdef __SSE2__
            _MM_SET_FLUSH_ZERO_MODE(oldMode);
#endif
        }
    }
}

double accumulateProduct(const float* data1, const float* data2, size_t n, bool multiThread) {
    if (n == 0) {
        return 0.0;
    }

    // use two accumulators to reduce dependencies (improves speed) and increase accuracy
    double acc1 = 0.0;
    double acc2 = 0.0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:acc1,acc2) if(multiThread)
#endif
    for (size_t i = 0; i < n - 1; i += 2) {
        acc1 += static_cast<double>(data1[i]) * static_cast<double>(data2[i]);
        acc2 += static_cast<double>(data1[i + 1]) * static_cast<double>(data2[i + 1]);
    }

    if (n & 1) {
        acc1 += static_cast<double>(data1[n -1]) * static_cast<double>(data2[n -1]);
    }
    return acc1 + acc2;
}
}
