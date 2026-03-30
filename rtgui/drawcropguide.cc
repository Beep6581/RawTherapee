/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2025 Daniel Gao <daniel.gao.work@gmail.com>
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

#include "drawcropguide.h"

#include "options.h"

#include "rtengine/aspectratios.h"
#include "rtengine/procparams.h"
#include "rtengine/rt_math.h"

using namespace rtengine;
using namespace rtengine::procparams;

namespace
{

using PresetIndex = CropGuideParams::PresetIndex;

constexpr double GUIDE_ALPHA = 0.618;
constexpr double GOLDEN_RATIO = 1.618033988749895;
constexpr double GOLDEN_RATIO_RECIPROCAL = 0.618033988749895;

struct CropRect {
    double x0;
    double y0;
    double x1;
    double y1;
};

enum class Dir { UP, RIGHT, DOWN, LEFT };

struct GuideDrawer {
    const Cairo::RefPtr<Cairo::Context>& cr;
    const CropRect& rect;
    const CropGuideParams& params;

    double red = 1.0;
    double green = 1.0;
    double blue = 1.0;

    GuideDrawer(const Cairo::RefPtr<Cairo::Context>& context,
                const CropRect& bounds,
                const CropGuideParams& guideParams)
        : cr(context), rect(bounds), params(guideParams) {}

    void drawRuleOfThirds();
    void drawHarmonicMeans();
    void drawCrosshair();
    void drawGrid();
    void drawEpassport();
    void drawCenteredSquare();

    void drawDiagonals();
    void drawGoldenTriangle();

    void drawGoldenRatio(double device_scale);
    void drawAspectRatios();

    void drawHorizontal(double ratio);
    void drawVertical(double ratio);

    void drawDashedLine(double x0, double y0, double x1, double y1);
    void drawDashedArc(double xc, double yc, double r, double rad0, double rad1);

private:
    void drawGoldenRatioRecursive(double x, double y, double length, Dir dir,
                                  double limit, bool clockwise);
};

void drawFrame(const Cairo::RefPtr<Cairo::Context>& cr,
               const CropRect& rect)
{
    cr->set_line_width(1.0);
    cr->set_source_rgba(1.0, 1.0, 1.0, GUIDE_ALPHA);
    cr->move_to(rect.x0, rect.y0);
    cr->line_to(rect.x1, rect.y0);
    cr->line_to(rect.x1, rect.y1);
    cr->line_to(rect.x0, rect.y1);
    cr->line_to(rect.x0, rect.y0);
    cr->stroke();

    cr->set_source_rgba(0.0, 0.0, 0.0, GUIDE_ALPHA);
    cr->set_dash(std::valarray<double>({4}), 0);
    cr->move_to(rect.x0, rect.y0);
    cr->line_to(rect.x1, rect.y0);
    cr->line_to(rect.x1, rect.y1);
    cr->line_to(rect.x0, rect.y1);
    cr->line_to(rect.x0, rect.y0);
    cr->stroke();
    cr->set_dash(std::valarray<double>(), 0);
}

void drawCropGuides(const Cairo::RefPtr<Cairo::Context>& cr,
                    const CropRect& crop_rect,
                    const CropRect& bleed_rect,
                    const CropGuideParams& params,
                    CropGuideOverride override,
                    double device_scale)
{
    switch (override) {
        case CropGuideOverride::FRAME:
            break;
        case CropGuideOverride::DONT_TOUCH:
            if (!params.enabled) {
                return; // Skip drawing guides
            } else {
                break;
            }
        case CropGuideOverride::NO_GUIDES:
        default:
            return; // Skip drawing guides
    }

    drawFrame(cr, crop_rect);

    if (params.bleed > 0) {
        drawFrame(cr, bleed_rect);
    }

    if (override == CropGuideOverride::FRAME) return;

    GuideDrawer util{cr, bleed_rect, params};

    if (params.presets.any()) {
        // Horizontal/vertical lines only
        util.drawRuleOfThirds();
        util.drawHarmonicMeans();
        util.drawCrosshair();
        util.drawGrid();
        util.drawEpassport();
        util.drawCenteredSquare();

        // Diagonals
        util.drawDiagonals();
        util.drawGoldenTriangle();

        util.drawGoldenRatio(device_scale);
    }

    util.drawAspectRatios();
}

CropRect calculateScaleBleedRect(const CropRect& crop_rect, double bleed)
{
    // outer = inner + 2 * bleed * inner = (1 + 2 * bleed) * inner
    double scale = (1.0 + 2.0 * bleed / 100.0);

    double w = crop_rect.x1 - crop_rect.x0;
    double h = crop_rect.y1 - crop_rect.y0;
    double new_w = w / scale;
    double new_h = h / scale;

    double dx = (w - new_w) / 2.0;
    double dy = (h - new_h) / 2.0;

    CropRect bleed_rect = crop_rect;
    bleed_rect.x0 += dx;
    bleed_rect.y0 += dy;
    bleed_rect.x1 -= dx;
    bleed_rect.y1 -= dy;

    return bleed_rect;
}

CropRect calculateUniformBleedRect(const CropRect& crop_rect, double bleed, bool use_width)
{
    // outer = inner + 2 * bleed * inner = (1 + 2 * bleed) * inner
    double scale = (1.0 + 2.0 * bleed / 100.0);

    double w = crop_rect.x1 - crop_rect.x0;
    double h = crop_rect.y1 - crop_rect.y0;

    double new_w = 0.0;
    double new_h = 0.0;
    if (use_width) {
        new_w = w / scale;
        new_h = std::max(0.0, h - (w - new_w));
    } else {
        new_h = h / scale;
        new_w = std::max(0.0, w - (h - new_h));
    }

    double dx = (w - new_w) / 2.0;
    double dy = (h - new_h) / 2.0;

    CropRect bleed_rect = crop_rect;
    bleed_rect.x0 += dx;
    bleed_rect.y0 += dy;
    bleed_rect.x1 -= dx;
    bleed_rect.y1 -= dy;

    return bleed_rect;
}

CropRect calculateBleedRect(const CropRect& crop_rect, const CropGuideParams& params)
{
    if (params.bleed <= 0) return crop_rect;

    double w = crop_rect.x1 - crop_rect.x0;
    double h = crop_rect.y1 - crop_rect.y0;

    switch (params.basis) {
        case CropGuideParams::Basis::WIDTH:
            return calculateUniformBleedRect(crop_rect, params.bleed, true);
        case CropGuideParams::Basis::HEIGHT:
            return calculateUniformBleedRect(crop_rect, params.bleed, false);
        case CropGuideParams::Basis::LONG:
            return calculateUniformBleedRect(crop_rect, params.bleed, w > h);
        case CropGuideParams::Basis::SHORT:
            return calculateUniformBleedRect(crop_rect, params.bleed, w < h);
        case CropGuideParams::Basis::SCALE:
        default:
            return calculateScaleBleedRect(crop_rect, params.bleed);
    }
}

}  // namespace

void GuideDrawer::drawRuleOfThirds()
{
    if (!params.presets[PresetIndex::RULE_OF_THIRDS]) return;

    drawHorizontal(1.0 / 3.0);
    drawHorizontal(2.0 / 3.0);
    drawVertical(1.0 / 3.0);
    drawVertical(2.0 / 3.0);
}

void GuideDrawer::drawHarmonicMeans()
{
    if (!params.presets[PresetIndex::HARMONIC_MEANS]) return;

    drawHorizontal(GOLDEN_RATIO_RECIPROCAL);
    drawHorizontal(1.0 - GOLDEN_RATIO_RECIPROCAL);
    drawVertical(GOLDEN_RATIO_RECIPROCAL);
    drawVertical(1.0 - GOLDEN_RATIO_RECIPROCAL);
}

void GuideDrawer::drawCrosshair()
{
    if (!params.presets[PresetIndex::CROSSHAIR]) return;

    drawHorizontal(0.5);
    drawVertical(0.5);
}

void GuideDrawer::drawGrid()
{
    if (!params.presets[PresetIndex::GRID]) return;

    // To have even distribution, normalize it a bit
    const int longSideNumLines = 10;

    const double w = rect.x1 - rect.x0;
    const double h = rect.y1 - rect.y0;

    if (w > longSideNumLines && h > longSideNumLines) {
        if (w > h) {
            for (int i = 1; i < longSideNumLines; i++) {
                drawVertical(static_cast<double>(i) / longSideNumLines);
            }

            int shortSideNumLines = static_cast<int>(std::round(h * longSideNumLines / w));
            for (int i = 1; i < shortSideNumLines; i++) {
                drawHorizontal(static_cast<double>(i) / shortSideNumLines);
            }
        } else {
            for (int i = 1; i < longSideNumLines; i++) {
                drawHorizontal(static_cast<double>(i) / longSideNumLines);
            }

            int shortSideNumLines = static_cast<int>(std::round(w * longSideNumLines / h));
            for (int i = 1; i < shortSideNumLines; i++) {
                drawVertical(static_cast<double>(i) / shortSideNumLines);
            }
        }
    }
}

/**
 * Official measurements do not specify exact ratios, just min/max measurements
 * within which the eyes and chin-crown distance must lie. I averaged those
 * measurements to produce these guides.
 *
 * The first horizontal guide is for the crown, the second is roughly for the
 * nostrils, the third is for the chin.
 *
 * http://www.homeoffice.gov.uk/agencies-public-bodies/ips/passports/information-photographers/
 * "(...) the measurement of the face from the bottom of the chin to the crown
 * (ie the top of the head, not the top of the hair) is between 29mm and 34mm."
 */
void GuideDrawer::drawEpassport()
{
    if (!params.presets[PresetIndex::EPASSPORT]) return;

    drawHorizontal(7.0 / 45.0);
    drawHorizontal(26.0 / 45.0);
    drawHorizontal(37.0 / 45.0);
    drawVertical(0.5);
}

void GuideDrawer::drawCenteredSquare()
{
    if (!params.presets[PresetIndex::CENTERED_SQUARE]) return;

    const double w = rect.x1 - rect.x0;
    const double h = rect.y1 - rect.y0;
    double ratio = w / h;

    if (ratio < 1.0) {
        ratio = 1.0 / ratio;
        const double pos = (ratio - 1.0) / (2.0 * ratio);
        drawHorizontal(pos);
        drawHorizontal(1.0 - pos);
    } else {
        const double pos = (ratio - 1.0) / (2.0 * ratio);
        drawVertical(pos);
        drawVertical(1.0 - pos);
    }
}

void GuideDrawer::drawDiagonals()
{
    if (!params.presets[PresetIndex::RULE_OF_DIAGONALS]) return;

    double corners_from[4][2];
    double corners_to[4][2];

    double mindim = std::min(rect.x1 - rect.x0, rect.y1 - rect.y0);

    // clang-format off
    corners_from[0][0] = rect.x0;
    corners_from[0][1] = rect.y0;
    corners_to[0][0]   = rect.x0 + mindim;
    corners_to[0][1]   = rect.y0 + mindim;
    corners_from[1][0] = rect.x0;
    corners_from[1][1] = rect.y1;
    corners_to[1][0]   = rect.x0 + mindim;
    corners_to[1][1]   = rect.y1 - mindim;
    corners_from[2][0] = rect.x1;
    corners_from[2][1] = rect.y0;
    corners_to[2][0]   = rect.x1 - mindim;
    corners_to[2][1]   = rect.y0 + mindim;
    corners_from[3][0] = rect.x1;
    corners_from[3][1] = rect.y1;
    corners_to[3][0]   = rect.x1 - mindim;
    corners_to[3][1]   = rect.y1 - mindim;
    // clang-format on

    for (int i = 0; i < 4; i++) {
        drawDashedLine(corners_from[i][0], corners_from[i][1],
                       corners_to[i][0], corners_to[i][1]);
    }
}

void GuideDrawer::drawGoldenTriangle()
{
    if (!params.presets[PresetIndex::GOLDEN_TRIANGLE]) return;

    double x0 = rect.x0;
    double x1 = rect.x1;
    const double y0 = rect.y0;
    const double y1 = rect.y1;

    if (params.mirror_golden_triangle) {
        std::swap(x0, x1);
    }

    drawDashedLine(x0, y0, x1, y1);

    const double height = y1 - y0;
    const double width = x1 - x0;
    const double d = std::sqrt(height * height + width * width);
    const double alpha = std::asin(width / d);
    const double beta = std::asin(height / d);
    const double a = std::sin(beta) * height;
    const double b = std::sin(alpha) * height;

    double x = (a * b) / height;
    double y = height - (b * (d - a)) / width;
    drawDashedLine(x0, y1, x0 + x, y0 + y);

    x = width - (a * b) / height;
    y = (b * (d - a)) / width;
    drawDashedLine(x1, y0, x0 + x, y0 + y);
}

void GuideDrawer::drawGoldenRatio(double device_scale)
{
    if (!params.presets[PresetIndex::GOLDEN_RATIO]) return;

    const double w = rect.x1 - rect.x0;
    const double h = rect.y1 - rect.y0;
    const double aspect_ratio = w / h;

    // If the current ratio doesn't match the golden ratio, use a fitted and
    // centered guide that does match the golden ratio.
    double fitted_w = w;
    double fitted_h = h;
    if (w >= h) {
        if (aspect_ratio >= GOLDEN_RATIO) {
            fitted_w = h * GOLDEN_RATIO;
            fitted_h = h;
        } else {
            fitted_w = w;
            fitted_h = w * GOLDEN_RATIO_RECIPROCAL;
        }
    } else {
        if (aspect_ratio <= GOLDEN_RATIO_RECIPROCAL) {
            fitted_w = w;
            fitted_h = w * GOLDEN_RATIO;
        } else {
            fitted_w = h * GOLDEN_RATIO_RECIPROCAL;
            fitted_h = h;
        }
    }

    const double delta_w = w - fitted_w;
    const double delta_h = h - fitted_h;

    const CropRect fitted_rect {
        rect.x0 + (delta_w / 2.0),
        rect.y0 + (delta_h / 2.0),
        rect.x1 - (delta_w / 2.0),
        rect.y1 - (delta_h / 2.0)
    };

    // Don't show new fitted rect bounds if close enough
    if (delta_w >= 4.0) {
        drawDashedLine(fitted_rect.x0, fitted_rect.y0, fitted_rect.x0, fitted_rect.y1);
        drawDashedLine(fitted_rect.x1, fitted_rect.y0, fitted_rect.x1, fitted_rect.y1);
    }
    if (delta_h >= 4.0) {
        drawDashedLine(fitted_rect.x0, fitted_rect.y0, fitted_rect.x1, fitted_rect.y0);
        drawDashedLine(fitted_rect.x0, fitted_rect.y1, fitted_rect.x1, fitted_rect.y1);
    }

    double start_x = fitted_rect.x0;
    double start_y = fitted_rect.y1;
    Dir dir = Dir::RIGHT;

    if (fitted_w >= fitted_h) {
        if (params.rotate_golden_ratio && params.mirror_golden_ratio) {
            start_x = fitted_rect.x0;
            start_y = fitted_rect.y0;
            dir = Dir::RIGHT;
        } else if (params.rotate_golden_ratio) {
            start_x = fitted_rect.x1;
            start_y = fitted_rect.y0;
            dir = Dir::LEFT;
        } else if (params.mirror_golden_ratio) {
            start_x = fitted_rect.x1;
            start_y = fitted_rect.y1;
            dir = Dir::LEFT;
        } else {
            start_x = fitted_rect.x0;
            start_y = fitted_rect.y1;
            dir = Dir::RIGHT;
        }
    } else {
        if (params.rotate_golden_ratio && params.mirror_golden_ratio) {
            start_x = fitted_rect.x0;
            start_y = fitted_rect.y1;
            dir = Dir::UP;
        } else if (params.rotate_golden_ratio) {
            start_x = fitted_rect.x1;
            start_y = fitted_rect.y1;
            dir = Dir::UP;
        } else if (params.mirror_golden_ratio) {
            start_x = fitted_rect.x1;
            start_y = fitted_rect.y0;
            dir = Dir::DOWN;
        } else {
            start_x = fitted_rect.x0;
            start_y = fitted_rect.y0;
            dir = Dir::DOWN;
        }
    }

    constexpr double LIMIT_THRESHOLD = 0.05;
    double length = fitted_w;
    double limit = fitted_h * LIMIT_THRESHOLD;
    if (fitted_w < fitted_h) {
        length = fitted_h;
        limit = fitted_w * LIMIT_THRESHOLD;
    }
    limit = std::max(limit, 10.0 * device_scale);

    bool clockwise = !params.mirror_golden_ratio;
    drawGoldenRatioRecursive(start_x, start_y, length, dir, limit, clockwise);
}

/**
 * @param x Current spiral arc start point x
 * @param y Current spiral arc start point y
 * @param length Current length along dir that will be split into golden ratio
 * @param dir Current direction of split (i.e. where smaller section is positioned)
 * @param limit Min length before stopping recursion
 * @param clockwise Direction of spiral while travelling inwards
 */
void GuideDrawer::drawGoldenRatioRecursive(double x, double y, double length,
                                           Dir dir, double limit, bool clockwise)
{
    if (length <= limit) return;

    double offset = length * GOLDEN_RATIO_RECIPROCAL;

    // Position along square line segment of split.
    // This is also the center of the spiral's circular arc segment.
    double arc_x = x;
    double arc_y = y;
    // In radians, the angles of the spiral's circular arc segment.
    double arc_start = 2 * RT_PI;
    double arc_end = 0;
    // The end point of the spiral after this iteration.
    double next_x = x;
    double next_y = y;
    Dir next_dir = dir;

    // These angles are if you are viewing the image on a normal x/y graph.
    // (i.e. they don't match cairo angles)
    constexpr double ARC_0 = 0.0;
    constexpr double ARC_90 = 3.0 * RT_PI / 2.0;
    constexpr double ARC_180 = RT_PI;
    constexpr double ARC_270 = RT_PI / 2.0;

    if (clockwise) {
        switch (dir) {
            case Dir::UP:
                arc_y -= offset;
                next_x -= offset;
                next_y -= offset;
                arc_start = ARC_270;
                arc_end = ARC_180;
                next_dir = Dir::RIGHT;
                break;
            case Dir::RIGHT:
                arc_x += offset;
                next_x += offset;
                next_y -= offset;
                arc_start = ARC_180;
                arc_end = ARC_90;
                next_dir = Dir::DOWN;
                break;
            case Dir::DOWN:
                arc_y += offset;
                next_x += offset;
                next_y += offset;
                arc_start = ARC_90;
                arc_end = ARC_0;
                next_dir = Dir::LEFT;
                break;
            case Dir::LEFT:
                arc_x -= offset;
                next_x -= offset;
                next_y += offset;
                arc_start = ARC_0;
                arc_end = ARC_270;
                next_dir = Dir::UP;
                break;
        }
    } else {
        switch (dir) {
            case Dir::UP:
                arc_y -= offset;
                next_x += offset;
                next_y -= offset;
                arc_start = ARC_0;
                arc_end = ARC_270;
                next_dir = Dir::LEFT;
                break;
            case Dir::RIGHT:
                arc_x += offset;
                next_x += offset;
                next_y += offset;
                arc_start = ARC_270;
                arc_end = ARC_180;
                next_dir = Dir::UP;
                break;
            case Dir::DOWN:
                arc_y += offset;
                next_x -= offset;
                next_y += offset;
                arc_start = ARC_180;
                arc_end = ARC_90;
                next_dir = Dir::RIGHT;
                break;
            case Dir::LEFT:
                arc_x -= offset;
                next_x -= offset;
                next_y -= offset;
                arc_start = ARC_90;
                arc_end = ARC_0;
                next_dir = Dir::DOWN;
                break;
        }
    }

    drawDashedArc(arc_x, arc_y, offset, arc_start, arc_end);
    drawDashedLine(arc_x, arc_y, next_x, next_y);

    drawGoldenRatioRecursive(next_x, next_y, offset, next_dir, limit, clockwise);
}

void GuideDrawer::drawAspectRatios()
{
    const double w = rect.x1 - rect.x0;
    const double h = rect.y1 - rect.y0;
    const double rect_ratio = w / h;

    for (const auto& entry : params.aspect_ratios) {
        if (!entry.enabled) continue;

        red = entry.red;
        green = entry.green;
        blue = entry.blue;

        const double aspect_ratio = [&]() {
            double result = getAspectRatioValue(entry.preset_index);
            if ((entry.is_portrait && result > 1.0)
                || (!entry.is_portrait && result < 1.0))
            {
                result = 1.0 / result;
            }
            return result;
        }();

        // If the current ratio doesn't match the desired aspect ratio, use a
        // fitted and centered guide that does match the aspect ratio.
        double fitted_w = w;
        double fitted_h = h;
        if (rect_ratio >= aspect_ratio) {
            fitted_w = h * aspect_ratio;
            fitted_h = h;
        } else {
            fitted_w = w;
            fitted_h = w / aspect_ratio;
        }

        const double delta_w = w - fitted_w;
        const double delta_h = h - fitted_h;

        const CropRect fitted_rect {
            rect.x0 + (delta_w / 2.0),
            rect.y0 + (delta_h / 2.0),
            rect.x1 - (delta_w / 2.0),
            rect.y1 - (delta_h / 2.0)
        };

        // Don't show new fitted rect bounds if close enough
        if (delta_w >= 2.0) {
            drawDashedLine(fitted_rect.x0, fitted_rect.y0, fitted_rect.x0, fitted_rect.y1);
            drawDashedLine(fitted_rect.x1, fitted_rect.y0, fitted_rect.x1, fitted_rect.y1);
        }
        if (delta_h >= 2.0) {
            drawDashedLine(fitted_rect.x0, fitted_rect.y0, fitted_rect.x1, fitted_rect.y0);
            drawDashedLine(fitted_rect.x0, fitted_rect.y1, fitted_rect.x1, fitted_rect.y1);
        }
    }

    // Restore default crop guide color
    red = 1.0;
    green = 1.0;
    blue = 1.0;
}

void GuideDrawer::drawHorizontal(double ratio)
{
    double y = rect.y0 + std::round((rect.y1 - rect.y0) * ratio);
    drawDashedLine(rect.x0, y, rect.x1, y);
}

void GuideDrawer::drawVertical(double ratio)
{
    double x = rect.x0 + std::round((rect.x1 - rect.x0) * ratio);
    drawDashedLine(x, rect.y0, x, rect.y1);
}

void GuideDrawer::drawDashedLine(double x0, double y0, double x1, double y1)
{
    cr->set_source_rgba(red, green, blue, GUIDE_ALPHA);
    cr->move_to(x0, y0);
    cr->line_to(x1, y1);
    cr->stroke();

    cr->set_source_rgba(0.0, 0.0, 0.0, GUIDE_ALPHA);
    cr->set_dash(std::valarray<double>({4}), 0);
    cr->move_to(x0, y0);
    cr->line_to(x1, y1);
    cr->stroke();

    // Reset state
    cr->set_dash(std::valarray<double>(), 0);
}

void GuideDrawer::drawDashedArc(double xc, double yc, double r, double rad0, double rad1)
{
    cr->set_source_rgba(1.0, 1.0, 1.0, GUIDE_ALPHA);
    cr->arc(xc, yc, r, rad0, rad1);
    cr->stroke();

    cr->set_source_rgba(0.0, 0.0, 0.0, GUIDE_ALPHA);
    cr->set_dash(std::valarray<double>({4}), 0);
    cr->arc(xc, yc, r, rad0, rad1);
    cr->stroke();

    // Reset state
    cr->set_dash(std::valarray<double>(), 0);
}

void drawCrop(const Cairo::RefPtr<Cairo::Context>& cr,
              double imx, double imy, double imw, double imh,
              double clipWidth, double clipHeight,
              double startx, double starty,
              double scale, double deviceScale,
              const CropParams& cropParams,
              const CropGuideParams& cropGuideParams,
              CropGuideOverride cropGuideOverride,
              bool useBgColor, bool fullImageVisible)
{
    cr->save();

    cr->set_line_width(0.0);
    cr->rectangle(imx, imy, clipWidth, clipHeight);
    cr->clip();

    double c1x = (cropParams.x - startx) * scale;
    double c1y = (cropParams.y - starty) * scale;
    double c2x = (cropParams.x + cropParams.w - startx) * scale - (fullImageVisible ? 0.0 : 1.0);
    double c2y = (cropParams.y + cropParams.h - starty) * scale - (fullImageVisible ? 0.0 : 1.0);

    // crop overlay color, linked with crop windows background
    const auto& options = App::get().options();
    if (options.bgcolor == 0 || !useBgColor) {
        cr->set_source_rgba (options.cutOverlayBrush[0], options.cutOverlayBrush[1], options.cutOverlayBrush[2], options.cutOverlayBrush[3]);
    } else if (options.bgcolor == 1) {
        cr->set_source_rgb (0, 0, 0);
    } else if (options.bgcolor == 2) {
        cr->set_source_rgb (1, 1, 1);
    } else if (options.bgcolor == 3) {
        cr->set_source_rgb (0.467, 0.467, 0.467);
    }

    cr->rectangle(imx, imy, imw + 0.5, std::round(c1y) + 0.5);
    cr->rectangle(imx, std::round(imy + c2y) + 0.5,
                  imw + 0.5, std::round(imh - c2y) + 0.5);
    cr->rectangle(imx, std::round(imy + c1y) + 0.5,
                  std::round(c1x) + 0.5, std::round(c2y - c1y + 1) + 0.5);
    cr->rectangle(std::round(imx + c2x) + 0.5, std::round(imy + c1y) + 0.5,
                  std::round(imw - c2x) + 0.5, std::round(c2y - c1y + 1) + 0.5);
    cr->fill();

    cr->restore();

    if (cropGuideOverride == CropGuideOverride::NO_GUIDES) return;

    // Rectangle around the cropped area for drawing crop guide frame
    CropRect cropRect {
        /*.x0 =*/ std::round(c1x) + imx + 0.5,
        /*.y0 =*/ std::round(c1y) + imy + 0.5,
        /*.x1 =*/ std::round(c2x) + imx + 0.5,
        /*.y1 =*/ std::round(c2y) + imy + 0.5
    };
    // If bleed is enabled, shrink the crop rect to simulate the bleed margin
    CropRect bleedRect = calculateBleedRect(cropRect, cropGuideParams);

    if (fullImageVisible) {
        cropRect.x1 = std::min(cropRect.x1, imx + imw - 0.5);
        cropRect.y1 = std::min(cropRect.y1, imy + imh - 0.5);

        bleedRect.x1 = std::min(bleedRect.x1, imx + imw - 0.5);
        bleedRect.y1 = std::min(bleedRect.y1, imy + imh - 0.5);
    }

    drawCropGuides(cr, cropRect, bleedRect, cropGuideParams, cropGuideOverride,
                   deviceScale);
}
