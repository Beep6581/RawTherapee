/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2024 Rawtherapee developers
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
#include <array>
#include <iostream>

#include "lensmetadata.h"

namespace rtengine
{

CenterRadiusMetadataLensCorrection::CenterRadiusMetadataLensCorrection(const FramesMetaData *meta) :
    swap_xy(false)
{
    metadata = Exiv2Metadata(meta->getFileName());
    metadata.load();
}

void CenterRadiusMetadataLensCorrection::initCorrections(int width, int height, const procparams::CoarseTransformParams &coarse, int rawRotationDeg)
{
    if (rawRotationDeg >= 0) {
        int rot = (coarse.rotate + rawRotationDeg) % 360;
        swap_xy = (rot == 90 || rot == 270);
        if (swap_xy) {
            std::swap(width, height);
        }
    }

    w2 = width * 0.5f;
    h2 = height * 0.5f;
    rf = 1 / std::sqrt(SQR(w2) + SQR(h2));
}

void CenterRadiusMetadataLensCorrection::process(double &x, double &y, int cx, int cy, int channel, bool dist, bool ca) const
{
    double xx = x + cx;
    double yy = y + cy;
    if (swap_xy) {
        std::swap(xx, yy);
    }

    double xc = xx - w2;
    double yc = yy - h2;

    double rout = rf * std::sqrt(SQR(xc) + SQR(yc));
    double cf = 1;
    if (dist && ca) {
        cf = distortionAndCACorrectionFactor(rout, channel);
    } else if (dist) {
        cf = distortionCorrectionFactor(rout);
    } else if (ca) {
        cf = caCorrectionFactor(rout, channel);
    }

    x = cf * xc + w2;
    y = cf * yc + h2;

    if (swap_xy) {
        std::swap(x, y);
    }
    x -= cx;
    y -= cy;
}

void CenterRadiusMetadataLensCorrection::correctDistortionAndCA(double &x, double &y, int cx, int cy, int channel) const
{
    if (!hasDistortionCorrection() || !hasCACorrection()) {
        return;
    }

    process(x, y, cx, cy, channel, true, true);
}

void CenterRadiusMetadataLensCorrection::correctDistortion(double &x, double &y, int cx, int cy) const
{
    if (!hasDistortionCorrection()) {
        return;
    }

    process(x, y, cx, cy, 1, true, false);
}

void CenterRadiusMetadataLensCorrection::correctCA(double &x, double &y, int cx, int cy, int channel) const
{
    if (!hasCACorrection()) {
        return;
    }

    process(x, y, cx, cy, channel, false, true);
}

void CenterRadiusMetadataLensCorrection::processVignetteNChannels(int width, int height, float **rawData, int channels) const
{
    if (!hasVignettingCorrection()) {
        return;
    }

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double xc = x - w2;
            double yc = y - h2;
            double sf = vignettingCorrectionFactor(rf * std::sqrt(SQR(xc) + SQR(yc)));
            for (int c = 0; c < channels; c++) {
                rawData[y][x + c] *= sf;
            }
        }
    }
}

void CenterRadiusMetadataLensCorrection::processVignette(int width, int height, float **rawData) const
{
    return processVignetteNChannels(width, height, rawData, 1);
}

void CenterRadiusMetadataLensCorrection::processVignette3Channels(int width, int height, float **rawData) const
{
    return processVignetteNChannels(width, height, rawData, 3);
}

std::unique_ptr<MetadataLensCorrection> MetadataLensCorrectionFinder::findCorrection(const FramesMetaData *meta)
{
    static const std::unordered_set<std::string> makers = {};

    std::string make = Glib::ustring(meta->getMake()).uppercase();

    if (makers.find(make) == makers.end()) {
        return nullptr;
    }

    std::unique_ptr<MetadataLensCorrection> correction;

    return correction;
}

} // namespace rtengine
