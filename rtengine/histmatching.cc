/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

#include <iostream>

#include "color.h"
#include "curves.h"
#include "improcfun.h"
#include "procparams.h"
#include "rawimagesource.h"
#include "rt_math.h"
#include "rtthumbnail.h"
#include "settings.h"

//#define BENCHMARK
#include "StopWatch.h"

namespace rtengine
{

namespace {

struct CdfInfo {
    std::vector<int> cdf;
    int min_val;
    int max_val;

    CdfInfo(): cdf(256), min_val(-1), max_val(-1) {}
};


CdfInfo getCdf(const IImage8 &img)
{
    CdfInfo ret;

    for (int y = 0; y < img.getHeight(); ++y) {
        for (int x = 0; x < img.getWidth(); ++x) {
            int lum = LIM(int(Color::rgbLuminance(float(img.r(y, x)), float(img.g(y, x)), float(img.b(y, x)))), 0, 255);
            ++ret.cdf[lum];
        }
    }

    int sum = 0;
    for (size_t i = 0; i < ret.cdf.size(); ++i) {
        if (ret.cdf[i] > 0) {
            if (ret.min_val < 0) {
                ret.min_val = i;
            }
            ret.max_val = i;
        }
        sum += ret.cdf[i];
        ret.cdf[i] = sum;
    }

    return ret;
}


int findMatch(int val, const std::vector<int> &cdf, int j)
{
    if (cdf[j] <= val) {
        for (; j < int(cdf.size()); ++j) {
            if (cdf[j] == val) {
                return j;
            } else if (cdf[j] > val) {
                return (cdf[j] - val <= val - cdf[j-1] ? j : j-1);
            }
        }
        return 255;
    } else {
        for (; j >= 0; --j) {
            if (cdf[j] == val) {
                return j;
            } else if (cdf[j] < val) {
                return (val - cdf[j] <= cdf[j+1] - val ? j : j+1);
            }
        }
        return 0;
    }
}


void mappingToCurve(const std::vector<int> &mapping, std::vector<double> &curve)
{
    curve.clear();

    int idx = 15;
    for (; idx < int(mapping.size()); ++idx) {
        if (mapping[idx] >= idx) {
            break;
        }
    }
    if (idx == int(mapping.size())) {
        for (idx = 1; idx < int(mapping.size())-1; ++idx) {
            if (mapping[idx] >= idx) {
                break;
            }
        }
    }

    auto coord = [](int v) -> double { return double(v)/255.0; };
    auto doit =
        [&](int start, int stop, int step, bool addstart, int maxdelta=0) -> void
        {
            if (!maxdelta) maxdelta = step * 2;
            int prev = start;
            if (addstart && mapping[start] >= 0) {
                curve.push_back(coord(start));
                curve.push_back(coord(mapping[start]));
            }
            for (int i = start; i < stop; ++i) {
                int v = mapping[i];
                if (v < 0) {
                    continue;
                }
                bool change = i > 0 && v != mapping[i-1];
                int diff = i - prev;
                if ((change && std::abs(diff - step) <= 1) || diff > maxdelta) {
                    curve.push_back(coord(i));
                    curve.push_back(coord(v));
                    prev = i;
                }
            }
        };

    curve.push_back(0.0);
    curve.push_back(0.0);

    int start = 0;
    while (start < idx && (mapping[start] < 0 || start < idx / 2)) {
        ++start;
    }

    const int npoints = 8;
    int step = std::max(int(mapping.size())/npoints, 1);
    int end = mapping.size();
    if (idx <= end / 3) {
        doit(start, idx, idx / 2, true);
        step = (end - idx) / 4;
        doit(idx, end, step, false, step);
    } else {
        doit(start, idx, idx > step ? step : idx / 2, true);
        doit(idx, end, step, idx - step > step / 2 && std::abs(curve[curve.size()-2] - coord(idx)) > 0.01);
    }

    if (curve.size() > 2 && (1 - curve[curve.size()-2] <= coord(step) / 3)) {
        curve.pop_back();
        curve.pop_back();
    }

    curve.push_back(1.0);
    curve.push_back(1.0);

    // we assume we are matching an S-shaped curve, so try to avoid
    // concavities in the upper part of the S
    const auto getpos =
        [](float x, float xa, float ya, float xb, float yb)
        {
            // line equation:
            // (x - xa) / (xb - xa) = (y - ya) / (yb - ya)
            return (x - xa) / (xb - xa) * (yb - ya) + ya;
        };
    idx = -1;
    for (ssize_t i = curve.size()-1; i > 0; i -= 2) {
        if (curve[i] <= 0.0) {
            idx = i+1;
            break;
        }
    }
    if (idx >= 0 && size_t(idx) < curve.size()) {
        // idx is the position of the first point in the upper part of the S
        // for each 3 consecutive points (xa, ya), (x, y), (xb, yb) we check
        // that y is above the point at x of the line between the other two
        // if this is not the case, we remove (x, y) from the curve
        while (size_t(idx+5) < curve.size()) {
            float xa = curve[idx];
            float ya = curve[idx+1];
            float x = curve[idx+2];
            float y = curve[idx+3];
            float xb = curve[idx+4];
            float yb = curve[idx+5];
            float yy = getpos(x, xa, ya, xb, yb);
            if (yy > y) {
                // we have to remove (x, y) from the curve
                curve.erase(curve.begin()+(idx+2), curve.begin()+(idx+4));
            } else {
                // move on to the next point
                idx += 2;
            }
        }
    }

    if (curve.size() < 4) {
        curve = { DCT_Linear }; // not enough points, fall back to linear
    } else {
        curve.insert(curve.begin(), DCT_Spline);
        DiagonalCurve c(curve);
        double gap = 0.05;
        double x = 0.0;
        curve = { DCT_CatumullRom };
        while (x < 1.0) {
            curve.push_back(x);
            curve.push_back(c.getVal(x));
            x += gap;
            gap *= 1.4;
        }
        curve.push_back(1.0);
        curve.push_back(c.getVal(1.0));
    }
}

} // namespace


void RawImageSource::getAutoMatchedToneCurve(const ColorManagementParams &cp, const procparams::RAWParams &rawParams, StandardObserver observer, std::vector<double> &outCurve)
{
    BENCHFUN

    if (settings->verbose) {
        std::cout << "performing histogram matching for " << getFileName() << " on the embedded thumbnail" << std::endl;
    }

    const auto same_profile =
        [](const ColorManagementParams &a, const ColorManagementParams &b) -> bool
        {
            return (a.inputProfile == b.inputProfile
                    && a.toneCurve == b.toneCurve
                    && a.applyLookTable == b.applyLookTable
                    && a.applyBaselineExposureOffset == b.applyBaselineExposureOffset
                    && a.applyHueSatMap == b.applyHueSatMap
                    && a.dcpIlluminant == b.dcpIlluminant);
        };

    if (!histMatchingCache.empty() && same_profile(*histMatchingParams, cp)) {
        if (settings->verbose) {
            std::cout << "tone curve found in cache" << std::endl;
        }
        outCurve = histMatchingCache;
        return;
    }

    outCurve = { DCT_Linear };

    int fw, fh;
    getFullSize(fw, fh, TR_NONE);
    if (getRotateDegree() == 90 || getRotateDegree() == 270) {
        std::swap(fw, fh);
    }
    int skip = 3;

    if (settings->verbose) {
        std::cout << "histogram matching: full raw image size is " << fw << "x" << fh << std::endl;
    }

    ProcParams neutral;
    neutral.icm = cp;
    neutral.raw.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST);
    neutral.raw.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FAST);
    neutral.icm.outputProfile = ColorManagementParams::NoICMString;

    std::unique_ptr<IImage8> source;
    {
        eSensorType sensor_type;
        int w = 0, h = 0;
        const std::unique_ptr<Thumbnail> thumb(Thumbnail::loadQuickFromRaw(getFileName(), sensor_type, w, h, 1, false, true, true));
        if (!thumb) {
            if (settings->verbose) {
                std::cout << "histogram matching: no thumbnail found, generating a neutral curve" << std::endl;
            }
            histMatchingCache = outCurve;
            *histMatchingParams = cp;
            return;
        } else if (w * 33 < fw || w * h < 19200) {
             // Some cameras have extremely small thumbs, for example Canon PowerShot A3100 IS has 128x96 thumbs.
             // For them we skip histogram matching.
             // With 160x120 thumbs from RICOH GR DIGITAL 2 it works fine, so we use 19200 as limit.
            if (settings->verbose) {
                std::cout << "histogram matching: the embedded thumbnail is too small: " << w << "x" << h << std::endl;
            }
            histMatchingCache = outCurve;
            *histMatchingParams = cp;
            return;
        }
        skip = LIM(skip * fh / h, 6, 10); // adjust the skip factor -- the larger the thumbnail, the less we should skip to get a good match
        source.reset(thumb->quickProcessImage(neutral, fh / skip, TI_Nearest));

        if (settings->verbose) {
            std::cout << "histogram matching: extracted embedded thumbnail" << std::endl;
        }
    }

    std::unique_ptr<IImage8> target;
    {
        eSensorType sensor_type;
        double scale;
        int w = fw / skip, h = fh / skip;
        const std::unique_ptr<Thumbnail> thumb(Thumbnail::loadFromRaw(getFileName(), sensor_type, w, h, 1, false, observer, false, &rawParams, true));
        if (!thumb) {
            if (settings->verbose) {
                std::cout << "histogram matching: raw decoding failed, generating a neutral curve" << std::endl;
            }
            histMatchingCache = outCurve;
            *histMatchingParams = cp;
            return;
        }
        target.reset(thumb->processImage(neutral, sensor_type, fh / skip, TI_Nearest, getMetaData(), scale, false, true));

        int sw = source->getWidth(), sh = source->getHeight();
        int tw = target->getWidth(), th = target->getHeight();
        float thumb_ratio = float(std::max(sw, sh)) / float(std::min(sw, sh));
        float target_ratio = float(std::max(tw, th)) / float(std::min(tw, th));
        if (std::abs(thumb_ratio - target_ratio) > 0.01f) {
            int cx = 0, cy = 0;
            if (thumb_ratio > target_ratio) {
                // crop the height
                int ch = th - (tw * float(sh) / float(sw));
                cy += ch / 2;
                th -= ch;
            } else {
                // crop the width
                int cw = tw - (th * float(sw) / float(sh));
                cx += cw / 2;
                tw -= cw;
            }
            if (settings->verbose) {
                std::cout << "histogram matching: cropping target to get an aspect ratio of " << round(thumb_ratio * 100)/100.f << ":1, new size is " << tw << "x" << th << std::endl;
            }

            if (cx || cy) {
                Image8 *tmp = new Image8(tw, th);
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for (int y = 0; y < th; ++y) {
                    for (int x = 0; x < tw; ++x) {
                        tmp->r(y, x) = target->r(y+cy, x+cx);
                        tmp->g(y, x) = target->g(y+cy, x+cx);
                        tmp->b(y, x) = target->b(y+cy, x+cx);
                    }
                }
                target.reset(tmp);
            }
        }

        if (settings->verbose) {
            std::cout << "histogram matching: generated neutral rendering" << std::endl;
        }
    }
    if (target->getWidth() != source->getWidth() || target->getHeight() != source->getHeight()) {
        Image8 *tmp = new Image8(source->getWidth(), source->getHeight());
        target->resizeImgTo(source->getWidth(), source->getHeight(), TI_Nearest, tmp);
        target.reset(tmp);
    }
    CdfInfo scdf = getCdf(*source);
    CdfInfo tcdf = getCdf(*target);

    std::vector<int> mapping;
    int j = 0;
    for (int i = 0; i < int(tcdf.cdf.size()); ++i) {
        j = findMatch(tcdf.cdf[i], scdf.cdf, j);
        if (i >= tcdf.min_val && i <= tcdf.max_val && j >= scdf.min_val && j <= scdf.max_val) {
            mapping.push_back(j);
        } else {
            mapping.push_back(-1);
        }
    }

    mappingToCurve(mapping, outCurve);

    if (settings->verbose) {
        std::cout << "histogram matching: generated curve with " << outCurve.size()/2 << " control points" << std::endl;
    }

    histMatchingCache = outCurve;
    *histMatchingParams = cp;
}

} // namespace rtengine
