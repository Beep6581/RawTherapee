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

namespace
{

/* interpolateLinearSpline does a simple linear spline interpolation. Values
 * outside the external knots will return the value of the nearest knot without
 * any additional interpolation. */
double interpolateLinearSpline(const std::vector<double> &xi, const std::vector<double> &yi, double x)
{
    if (x < xi[0]) {
        return yi[0];
    }

    for (size_t i = 1; i < xi.size(); i++) {
        if (x >= xi[i - 1] && x <= xi[i]) {
            double dydx = (yi[i] - yi[i - 1]) / (xi[i] - xi[i - 1]);

            return yi[i - 1] + (x - xi[i - 1]) * dydx;
        }
    }

    return yi[yi.size() - 1];
}

} // namespace

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

/* Fuji, Sony, Olympus metadata handling and algorithms adapted from
 * - src/iop/lens.cc
 * - src/common/exif.cc
 * in darktable 4.6 */
/*
    This file is part of darktable,
    Copyright (C) 2019-2024 darktable developers.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

class SonyMetadataLensCorrection : public CenterRadiusMetadataLensCorrection
{
public:
    SonyMetadataLensCorrection(const FramesMetaData *meta) :
        CenterRadiusMetadataLensCorrection(meta)
    {
        parse();
        setup();
    }

private:
    int nc;
    std::array<short, 16> distortion;
    std::array<short, 16> ca_r;
    std::array<short, 16> ca_b;
    std::array<short, 16> vignetting;

    std::vector<double> knots;
    std::vector<double> dist;
    std::array<std::vector<double>, 3> ca;
    std::vector<double> vig;

    void parse()
    {
        if (Exiv2::versionNumber() < EXIV2_MAKE_VERSION(0, 27, 4)) {
            throw std::runtime_error("cannot get Sony correction data, too old exiv2 version " + Exiv2::versionString());
        }

        auto &exif = metadata.exifData();

        auto posd = exif.findKey(Exiv2::ExifKey("Exif.SubImage1.DistortionCorrParams"));
        auto posc = exif.findKey(Exiv2::ExifKey("Exif.SubImage1.ChromaticAberrationCorrParams"));
        auto posv = exif.findKey(Exiv2::ExifKey("Exif.SubImage1.VignettingCorrParams"));

        /* Sony metadata corrections parameters define some splines with N knots */
        if (posd == exif.end() || posc == exif.end() || posv == exif.end()) {
            throw std::runtime_error("cannot get Sony correction data");
        }

        const int nc = to_long(posd);
        if (nc <= 16 && 2 * nc == to_long(posc) && nc == to_long(posv)) {
            this->nc = nc;
            for (int i = 0; i < nc; i++) {
                distortion[i] = to_long(posd, i + 1);
                ca_r[i] = to_long(posc, i + 1);
                ca_b[i] = to_long(posc, nc + i + 1);
                vignetting[i] = to_long(posv, i + 1);
            }
        } else {
            throw std::runtime_error("cannot get Sony correction data");
        }
    }

    void setup()
    {
        knots.resize(nc);
        dist.resize(nc);
        vig.resize(nc);
        for (int i = 0; i < 3; ++i) {
            ca[i].resize(nc);
        }

        for (int i = 0; i < this->nc; i++) {
            knots[i] = (i + 0.5) / (nc - 1);

            dist[i] = distortion[i] * powf(2, -14) + 1;

            ca[0][i] = ca[1][i] = ca[2][i] = 1.f;
            ca[0][i] *= ca_r[i] * powf(2, -21) + 1;
            ca[2][i] *= ca_b[i] * powf(2, -21) + 1;

            vig[i] = 1 / powf(2, 0.5f - powf(2, vignetting[i] * powf(2, -13) - 1));
        }
    }

    double distortionCorrectionFactor(double rout) const override
    {
        return interpolateLinearSpline(knots, dist, rout);
    }

    double caCorrectionFactor(double rout, int channel) const override
    {
        return interpolateLinearSpline(knots, ca[channel], rout);
    }

    double distortionAndCACorrectionFactor(double rout, int channel) const override
    {
        return distortionCorrectionFactor(rout) * caCorrectionFactor(rout, channel);
    }

    double vignettingCorrectionFactor(double r) const override
    {
        return interpolateLinearSpline(knots, vig, r);
    }

    bool hasDistortionCorrection() const override { return true; }
    bool hasVignettingCorrection() const override { return true; }
    bool hasCACorrection() const override { return true; }
};

std::unique_ptr<MetadataLensCorrection> MetadataLensCorrectionFinder::findCorrection(const FramesMetaData *meta)
{
    static const std::unordered_set<std::string> makers = {
        "SONY",
    };

    std::string make = Glib::ustring(meta->getMake()).uppercase();

    if (makers.find(make) == makers.end()) {
        return nullptr;
    }

    std::unique_ptr<MetadataLensCorrection> correction;

    try {
        if (make == "SONY") {
            correction.reset(new SonyMetadataLensCorrection(meta));
        }
    } catch (std::exception &exc) {
        if (settings->verbose) {
            std::cerr << "error parsing lens metadata: " << exc.what() << std::endl;
        }

        correction.reset(nullptr);
    }

    return correction;
}

} // namespace rtengine
