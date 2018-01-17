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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "rawimagesource.h"
#include "rtthumbnail.h"
#include "curves.h"
#include "color.h"
#include "rt_math.h"
#include "iccstore.h"
#include "../rtgui/mydiagonalcurve.h"
#include "improcfun.h"
#define BENCHMARK
#include "StopWatch.h"
#include <iostream>


namespace rtengine {

extern const Settings *settings;

namespace {

std::vector<int> getCdf(const IImage8 &img)
{
    std::vector<int> ret(256);
    for (int y = 0; y < img.getHeight(); ++y) {
        for (int x = 0; x < img.getWidth(); ++x) {
            int lum = LIM(0, int(Color::rgbLuminance(float(img.r(y, x)), float(img.g(y, x)), float(img.b(y, x)))), 255);
            ++ret[lum];
        }
    }

    int sum = 0;
    for (size_t i = 0; i < ret.size(); ++i) {
        sum += ret[i];
        ret[i] = sum;
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
    
    const int npoints = 8;
    int idx = 1;
    for (; idx < int(mapping.size()); ++idx) {
        if (mapping[idx] >= idx) {
            break;
        }
    }
    int step = max(int(mapping.size())/npoints, 1);

    auto coord = [](int v) -> double { return double(v)/255.0; };
    auto doit =
        [&](int start, int stop, int step) -> void
        {
            int prev = start;
            for (int i = start; i < stop; ++i) {
                int v = mapping[i];
                bool change = i > 0 && v != mapping[i-1];
                int diff = i - prev;
                if (change && std::abs(diff - step) <= 1) {
                    curve.emplace_back(coord(i));
                    curve.emplace_back(coord(v));
                    prev = i;
                }
            }
        };
    doit(0, idx, idx > step ? step : idx / 2);
    doit(idx, int(mapping.size()), step);
    if (curve[1] > 0.01) {
        curve.insert(curve.begin(), 0.0);
        curve.insert(curve.begin(), 0.0);
    }
    if (curve.back() < 0.99 || (1 - curve[curve.size()-2] > step / 512.0 && curve.back() < coord(mapping.back()))) {
        curve.emplace_back(1.0);
        curve.emplace_back(coord(mapping.back()));
    }
    curve.insert(curve.begin(), DCT_Spline);
}

} // namespace


void RawImageSource::getAutoMatchedToneCurve(std::vector<double> &outCurve)
{
    BENCHFUN
        
    if (settings->verbose) {
        std::cout << "performing histogram matching for " << getFileName() << " on the embedded thumbnail" << std::endl;
    }
    
    const int rheight = 200;
    ProcParams neutral;
    std::unique_ptr<IImage8> target;
    { 
        int tr = TR_NONE;
        int fw, fh;
        getFullSize(fw, fh, tr);
        int skip = fh / rheight;
        PreviewProps pp(0, 0, fw, fh, skip);
        ColorTemp currWB = getWB();
        std::unique_ptr<Imagefloat> image(new Imagefloat(int(fw / skip), int(fh / skip)));
        getImage(currWB, tr, image.get(), pp, neutral.toneCurve, neutral.raw);

        // this could probably be made faster -- ideally we would need to just
        // perform the transformation from camera space to the output space
        // (taking gamma into account), but I couldn't find anything
        // ready-made, so for now this will do. Remember the famous quote:
        // "premature optimization is the root of all evil" :-)
        convertColorSpace(image.get(), neutral.icm, currWB);
        ImProcFunctions ipf(&neutral);
        LabImage tmplab(image->getWidth(), image->getHeight());
        ipf.rgb2lab(*image, tmplab, neutral.icm.working);
        image.reset(ipf.lab2rgbOut(&tmplab, 0, 0, tmplab.W, tmplab.H, neutral.icm));
        target.reset(image->to8());

        if (settings->verbose) {
            std::cout << "histogram matching: generated neutral rendering" << std::endl;
        }
    }
    std::unique_ptr<IImage8> source;
    {
        RawMetaDataLocation rml;
        eSensorType sensor_type;
        int w, h;
        std::unique_ptr<Thumbnail> thumb(Thumbnail::loadQuickFromRaw(getFileName(), rml, sensor_type, w, h, 1, false, true));
        source.reset(thumb->quickProcessImage(neutral, target->getHeight(), TI_Nearest));

        if (settings->verbose) {
            std::cout << "histogram matching: extracted embedded thumbnail" << std::endl;
        }
    }
    if (target->getWidth() != source->getWidth() || target->getHeight() != source->getHeight()) {
        Image8 *tmp = new Image8(source->getWidth(), source->getHeight());
        target->resizeImgTo(source->getWidth(), source->getHeight(), TI_Nearest, tmp);
        target.reset(tmp);
    }
    std::vector<int> scdf = getCdf(*source);
    std::vector<int> tcdf = getCdf(*target);

    std::vector<int> mapping;
    int j = 0;
    for (size_t i = 0; i < tcdf.size(); ++i) {
        j = findMatch(tcdf[i], scdf, j);
        mapping.emplace_back(j);
    }

    mappingToCurve(mapping, outCurve);

    if (settings->verbose) {
        std::cout << "histogram matching: generated curve with " << outCurve.size()/2 << " control points" << std::endl;
    }
}

} // namespace rtengine
