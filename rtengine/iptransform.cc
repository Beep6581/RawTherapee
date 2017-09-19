/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include "rtengine.h"
#include "improcfun.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mytime.h"
#include "rt_math.h"
#include "sleef.c"
#include "rtlensfun.h"


using namespace std;

namespace
{

float pow3 (float x)
{
    return x * x * x;
}

float pow4 (float x)
{
    return (x * x) * (x * x);
}

float pown (float x, int n)
{

    switch (n) {
        case 0:
            return 1;

        case 2:
            return x * x;

        case 4:
            return pow4 (x);

        case 6:
            return (x * x) * pow4 (x);

        case 8:
            return pow4 (x) * pow4 (x);

        default:
            return pow_F (x, n);
    }
}

float normn (float a, float b, int n)
{
    switch (n) {
        case 2:
            return sqrtf (a * a + b * b);

        case 4:
            return sqrtf (sqrtf (pow4 (a) + pow4 (b)));

        case 6:
            return sqrtf (xcbrtf (pow3 (a) * pow3 (a) + pow3 (b) * pow3 (b)));

        case 8:
            return sqrtf (sqrtf (sqrtf (pow4 (a) * pow4 (a) + pow4 (b) * pow4 (b))));

        default:
            return pow_F (pown (a, n) + pown (b, n), 1.f / n);
    }
}


}

namespace rtengine
{
#undef CLIPTOC

#define CLIPTOC(a,b,c,d) ((a)>=(b)?((a)<=(c)?(a):(d=true,(c))):(d=true,(b)))

bool ImProcFunctions::transCoord (int W, int H, const std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef,
                                  const LensCorrection *pLCPMap)
{

    bool clipped = false;

    red.clear ();
    green.clear ();
    blue.clear ();

    if (!needsCA() && !needsDistortion() && !needsRotation() && !needsPerspective() && (!params->lensProf.useDist || pLCPMap == nullptr)) {
        for (size_t i = 0; i < src.size(); i++) {
            red.push_back   (Coord2D (src[i].x, src[i].y));
            green.push_back (Coord2D (src[i].x, src[i].y));
            blue.push_back  (Coord2D (src[i].x, src[i].y));
        }

        return clipped;
    }

    double oW = W, oH = H;
    double w2 = (double) oW  / 2.0 - 0.5;
    double h2 = (double) oH  / 2.0 - 0.5;
    double maxRadius = sqrt ( (double) ( oW * oW + oH * oH ) ) / 2;

    // auxiliary variables for distortion correction
    bool needsDist = needsDistortion();  // for performance
    double distAmount = params->distortion.amount;

    // auxiliary variables for rotation
    double cost = cos (params->rotate.degree * rtengine::RT_PI / 180.0);
    double sint = sin (params->rotate.degree * rtengine::RT_PI / 180.0);

    // auxiliary variables for vertical perspective correction
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * rtengine::RT_PI;
    double vpteta  = fabs (vpalpha - rtengine::RT_PI / 2) < 3e-4 ? 0.0 : acos ((vpdeg > 0 ? 1.0 : -1.0) * sqrt ((-oW * oW * tan (vpalpha) * tan (vpalpha) + (vpdeg > 0 ? 1.0 : -1.0) * oW * tan (vpalpha) * sqrt (16 * maxRadius * maxRadius + oW * oW * tan (vpalpha) * tan (vpalpha))) / (maxRadius * maxRadius * 8)));
    double vpcospt = (vpdeg >= 0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

    // auxiliary variables for horizontal perspective correction
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * rtengine::RT_PI;
    double hpteta  = fabs (hpalpha - rtengine::RT_PI / 2) < 3e-4 ? 0.0 : acos ((hpdeg > 0 ? 1.0 : -1.0) * sqrt ((-oH * oH * tan (hpalpha) * tan (hpalpha) + (hpdeg > 0 ? 1.0 : -1.0) * oH * tan (hpalpha) * sqrt (16 * maxRadius * maxRadius + oH * oH * tan (hpalpha) * tan (hpalpha))) / (maxRadius * maxRadius * 8)));
    double hpcospt = (hpdeg >= 0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

    double ascale = ascaleDef > 0 ? ascaleDef : (params->commonTrans.autofill ? getTransformAutoFill (oW, oH, pLCPMap) : 1.0);

    for (size_t i = 0; i < src.size(); i++) {
        double x_d = src[i].x, y_d = src[i].y;

        if (pLCPMap && params->lensProf.useDist) {
            pLCPMap->correctDistortion(x_d, y_d, 0, 0, ascale);
        } else {
            x_d *= ascale;
            y_d *= ascale;
        }

        x_d += ascale * (0 - w2);     // centering x coord & scale
        y_d += ascale * (0 - h2);     // centering y coord & scale

        if (needsPerspective()) {
            // horizontal perspective transformation
            y_d *= maxRadius / (maxRadius + x_d * hptanpt);
            x_d *= maxRadius * hpcospt / (maxRadius + x_d * hptanpt);

            // vertical perspective transformation
            x_d *= maxRadius / (maxRadius - y_d * vptanpt);
            y_d *= maxRadius * vpcospt / (maxRadius - y_d * vptanpt);
        }

        // rotate
        double Dx = x_d * cost - y_d * sint;
        double Dy = x_d * sint + y_d * cost;

        // distortion correction
        double s = 1;

        if (needsDist) {
            double r = sqrt (Dx * Dx + Dy * Dy) / maxRadius; // sqrt is slow
            s = 1.0 - distAmount + distAmount * r ;
        }

        // LCP CA is not reflected in preview (and very small), so don't add it here

        red.push_back (Coord2D (Dx * (s + params->cacorrection.red) + w2, Dy * (s + params->cacorrection.red) + h2));
        green.push_back (Coord2D (Dx * s + w2, Dy * s + h2));
        blue.push_back (Coord2D (Dx * (s + params->cacorrection.blue) + w2, Dy * (s + params->cacorrection.blue) + h2));
    }

    // Clip all points and track if they were any corrections
    for (size_t i = 0; i < src.size(); i++) {
        red[i].x = CLIPTOC (red[i].x, 0, W - 1, clipped);
        red[i].y = CLIPTOC (red[i].y, 0, H - 1, clipped);
        green[i].x = CLIPTOC (green[i].x, 0, W - 1, clipped);
        green[i].y = CLIPTOC (green[i].y, 0, H - 1, clipped);
        blue[i].x = CLIPTOC (blue[i].x, 0, W - 1, clipped);
        blue[i].y = CLIPTOC (blue[i].y, 0, H - 1, clipped);
    }

    return clipped;
}

// Transform all corners and critical sidelines of an image
bool ImProcFunctions::transCoord (int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef, const LensCorrection *pLCPMap)
{
    const int DivisionsPerBorder = 32;

    int x1 = x, y1 = y;
    int x2 = x1 + w - 1;
    int y2 = y1 + h - 1;

    // Build all edge points and half-way points
    std::vector<Coord2D> corners (8);
    corners[0].set (x1, y1);
    corners[1].set (x1, y2);
    corners[2].set (x2, y2);
    corners[3].set (x2, y1);
    corners[4].set ((x1 + x2) / 2, y1);
    corners[5].set ((x1 + x2) / 2, y2);
    corners[6].set (x1, (y1 + y2) / 2);
    corners[7].set (x2, (y1 + y2) / 2);

    // Add several steps inbetween
    int xstep = (x2 - x1) / DivisionsPerBorder;

    if (xstep < 1) {
        xstep = 1;
    }

    for (int i = x1 + xstep; i <= x2 - xstep; i += xstep) {
        corners.push_back (Coord2D (i, y1));
        corners.push_back (Coord2D (i, y2));
    }

    int ystep = (y2 - y1) / DivisionsPerBorder;

    if (ystep < 1) {
        ystep = 1;
    }

    for (int i = y1 + ystep; i <= y2 - ystep; i += ystep) {
        corners.push_back (Coord2D (x1, i));
        corners.push_back (Coord2D (x2, i));
    }

    std::vector<Coord2D> r, g, b;

    bool clipped = transCoord (W, H, corners, r, g, b, ascaleDef, pLCPMap);

    // Merge all R G Bs into one X/Y pool
    std::vector<Coord2D> transCorners;
    transCorners.insert (transCorners.end(), r.begin(), r.end());
    transCorners.insert (transCorners.end(), g.begin(), g.end());
    transCorners.insert (transCorners.end(), b.begin(), b.end());

    // find the min/max of all coordinates, so the borders
    double x1d = transCorners[0].x;

    for (size_t i = 1; i < transCorners.size(); i++)
        if (transCorners[i].x < x1d) {
            x1d = transCorners[i].x;
        }

    int x1v = (int) (x1d);

    double y1d = transCorners[0].y;

    for (size_t i = 1; i < transCorners.size(); i++)
        if (transCorners[i].y < y1d) {
            y1d = transCorners[i].y;
        }

    int y1v = (int) (y1d);

    double x2d = transCorners[0].x;

    for (size_t i = 1; i < transCorners.size(); i++)
        if (transCorners[i].x > x2d) {
            x2d = transCorners[i].x;
        }

    int x2v = (int)ceil (x2d);

    double y2d = transCorners[0].y;

    for (size_t i = 1; i < transCorners.size(); i++)
        if (transCorners[i].y > y2d) {
            y2d = transCorners[i].y;
        }

    int y2v = (int)ceil (y2d);

    xv = x1v;
    yv = y1v;
    wv = x2v - x1v + 1;
    hv = y2v - y1v + 1;

    return clipped;
}

void ImProcFunctions::transform (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH,
                                 const ImageMetaData *metadata,
                                 int rawRotationDeg, bool fullImage)
{
    double focalLen = metadata->getFocalLen();
    double focalLen35mm = metadata->getFocalLen35mm();
    float focusDist = metadata->getFocusDist();
    double fNumber = metadata->getFNumber();

    std::unique_ptr<const LensCorrection> pLCPMap;

    if (needsLensfun()) {
        pLCPMap = LFDatabase::findModifier(params->lensProf, metadata, oW, oH, params->coarse, rawRotationDeg);
    } else if (needsLCP()) { // don't check focal length to allow distortion correction for lenses without chip
        const std::shared_ptr<LCPProfile> pLCPProf = LCPStore::getInstance()->getProfile (params->lensProf.lcpFile);

        if (pLCPProf) {
            pLCPMap.reset(
                new LCPMapper (pLCPProf, focalLen, focalLen35mm,
                               focusDist, fNumber, false,
                               params->lensProf.useDist,
                               oW, oH, params->coarse, rawRotationDeg
                )
            );
        }
    }

    if (! (needsCA() || needsDistortion() || needsRotation() || needsPerspective() || needsLCP() || needsLensfun()) && (needsVignetting() || needsPCVignetting() || needsGradient())) {
        transformLuminanceOnly (original, transformed, cx, cy, oW, oH, fW, fH);
    } else {
        TransformMode mode;
        if (!needsCA() && scale != 1) {
            mode = TRANSFORM_PREVIEW;
        } else if (!fullImage) {
            mode = TRANSFORM_HIGH_QUALITY;
        } else {
            mode = TRANSFORM_HIGH_QUALITY_FULLIMAGE;
        }
        transformGeneral(mode, original, transformed, cx, cy, sx, sy, oW, oH, fW, fH, pLCPMap.get());
    }
}

// helper function
void ImProcFunctions::calcVignettingParams (int oW, int oH, const VignettingParams& vignetting, double &w2, double &h2, double& maxRadius, double &v, double &b, double &mul)
{
    // vignette center is a point with coordinates between -1 and +1
    double x = vignetting.centerX / 100.0;
    double y = vignetting.centerY / 100.0;

    // calculate vignette center in pixels
    w2 = (double) oW  / 2.0 - 0.5 + x * oW;
    h2 = (double) oH  / 2.0 - 0.5 + y * oH;

    // max vignette radius in pixels
    maxRadius = sqrt ( (double) ( oW * oW + oH * oH ) ) / 2.;

    // vignette variables with applied strength
    v = 1.0 + vignetting.strength * fabs (vignetting.amount) * 3.0 / 400.0;
    b = 1.0 + vignetting.radius * 7.0 / 100.0;
    mul = (1.0 - v) / tanh (b);
}

struct grad_params {
    bool angle_is_zero, transpose, bright_top;
    float ta, yc, xc;
    float ys, ys_inv;
    float scale, botmul, topmul;
    float top_edge_0;
    int h;
};
static void calcGradientParams (int oW, int oH, const GradientParams& gradient, struct grad_params& gp)
{
    int w = oW;
    int h = oH;
    double gradient_stops = gradient.strength;
    double gradient_span = gradient.feather / 100.0;
    double gradient_center_x = gradient.centerX / 200.0 + 0.5;
    double gradient_center_y = gradient.centerY / 200.0 + 0.5;
    double gradient_angle = gradient.degree / 180.0 * rtengine::RT_PI;
    //fprintf(stderr, "%f %f %f %f %f %d %d\n", gradient_stops, gradient_span, gradient_center_x, gradient_center_y, gradient_angle, w, h);

    // make 0.0 <= gradient_angle < 2 * rtengine::RT_PI
    gradient_angle = fmod (gradient_angle, 2 * rtengine::RT_PI);

    if (gradient_angle < 0.0) {
        gradient_angle += 2.0 * rtengine::RT_PI;
    }

    gp.bright_top = false;
    gp.transpose = false;
    gp.angle_is_zero = false;
    gp.h = h;
    double cosgrad = cos (gradient_angle);

    if (fabs (cosgrad) < 0.707) {
        // we transpose to avoid division by zero at 90 degrees
        // (actually we could transpose only for 90 degrees, but this way we avoid
        // division with extremely small numbers
        gp.transpose = true;
        gradient_angle += 0.5 * rtengine::RT_PI;
        double gxc = gradient_center_x;
        gradient_center_x = 1.0 - gradient_center_y;
        gradient_center_y = gxc;
    }

    gradient_angle = fmod (gradient_angle, 2 * rtengine::RT_PI);

    if (gradient_angle > 0.5 * rtengine::RT_PI && gradient_angle < rtengine::RT_PI) {
        gradient_angle += rtengine::RT_PI;
        gp.bright_top = true;
    } else if (gradient_angle >= rtengine::RT_PI && gradient_angle < 1.5 * rtengine::RT_PI) {
        gradient_angle -= rtengine::RT_PI;
        gp.bright_top = true;
    }

    if (fabs (gradient_angle) < 0.001 || fabs (gradient_angle - 2 * rtengine::RT_PI) < 0.001) {
        gradient_angle = 0;
        gp.angle_is_zero = true;
    }

    if (gp.transpose) {
        gp.bright_top = !gp.bright_top;
    }

    if (gp.transpose) {
        int tmp = w;
        w = h;
        h = tmp;
    }

    gp.scale = 1.0 / pow (2, gradient_stops);

    if (gp.bright_top) {
        gp.topmul = 1.0;
        gp.botmul = gp.scale;
    } else {
        gp.topmul = gp.scale;
        gp.botmul = 1.0;
    }

    gp.ta = tan (gradient_angle);
    gp.xc = w * gradient_center_x;
    gp.yc = h * gradient_center_y;
    gp.ys = sqrt ((float)h * h + (float)w * w) * (gradient_span / cos (gradient_angle));
    gp.ys_inv = 1.0 / gp.ys;
    gp.top_edge_0 = gp.yc - gp.ys / 2.0;

    if (gp.ys < 1.0 / h) {
        gp.ys_inv = 0;
        gp.ys = 0;
    }
}

static float calcGradientFactor (const struct grad_params& gp, int x, int y)
{
    if (gp.angle_is_zero) {
        int gy = gp.transpose ? x : y;

        if (gy < gp.top_edge_0) {
            return gp.topmul;
        } else if (gy >= gp.top_edge_0 + gp.ys) {
            return gp.botmul;
        } else {
            float val = ((float) (gy - gp.top_edge_0) * gp.ys_inv);

            if (gp.bright_top) {
                val = 1.f - val;
            }

            val *= rtengine::RT_PI_F_2;

            if (gp.scale < 1.f) {
                val = pow3 (xsinf (val));
            } else {
                val = 1.f - pow3 (xcosf (val));
            }

            return gp.scale + val * (1.0 - gp.scale);
        }
    } else {
        int gy = gp.transpose ? x : y;
        int gx = gp.transpose ? gp.h - y - 1 : x;
        float top_edge = gp.top_edge_0 - gp.ta * (gx - gp.xc);

        if (gy < top_edge) {
            return gp.topmul;
        } else if (gy >= top_edge + gp.ys) {
            return gp.botmul;
        } else {
            float val = ((float) (gy - top_edge) * gp.ys_inv);

            val = gp.bright_top ? 1.f - val : val;

            val *= rtengine::RT_PI_F_2;

            if (gp.scale < 1.f) {
                val = pow3 (xsinf (val));
            } else {
                val = 1.f - pow3 (xcosf (val));
            }

            return gp.scale + val * (1.0 - gp.scale);
        }
    }
}

struct pcv_params {
    float oe_a, oe_b, oe1_a, oe1_b, oe2_a, oe2_b;
    float ie_mul, ie1_mul, ie2_mul;
    float sepmix, feather;
    int w, h, x1, x2, y1, y2;
    int sep;
    bool is_super_ellipse_mode, is_portrait;
    float scale;
    float fadeout_mul;
};
static void calcPCVignetteParams (int fW, int fH, int oW, int oH, const PCVignetteParams& pcvignette, const CropParams &crop, struct pcv_params& pcv)
{

    // ellipse formula: (x/a)^2 + (y/b)^2 = 1
    double roundness = pcvignette.roundness / 100.0;
    pcv.feather = pcvignette.feather / 100.0;

    if (crop.enabled) {
        pcv.w = (crop.w * oW) / fW;
        pcv.h = (crop.h * oH) / fH;
        pcv.x1 = (crop.x * oW) / fW;
        pcv.y1 = (crop.y * oH) / fH;
        pcv.x2 = pcv.x1 + pcv.w;
        pcv.y2 = pcv.y1 + pcv.h;
    } else {
        pcv.x1 = 0, pcv.y1 = 0;
        pcv.x2 = oW, pcv.y2 = oH;
        pcv.w = oW;
        pcv.h = oH;
    }

    pcv.fadeout_mul = 1.0 / (0.05 * sqrtf (oW * oW + oH * oH));
    float short_side = (pcv.w < pcv.h) ? pcv.w : pcv.h;
    float long_side =  (pcv.w > pcv.h) ? pcv.w : pcv.h;

    pcv.sep = 2;
    pcv.sepmix = 0;
    pcv.oe_a = sqrt (2.0) * long_side * 0.5;
    pcv.oe_b = pcv.oe_a * short_side / long_side;
    pcv.ie_mul = (1.0 / sqrt (2.0)) * (1.0 - pcv.feather);
    pcv.is_super_ellipse_mode = false;
    pcv.is_portrait = (pcv.w < pcv.h);

    if (roundness < 0.5) {
        // make super-ellipse of higher and higher degree
        pcv.is_super_ellipse_mode = true;
        float sepf = 2 + 4 * powf (1.0 - 2 * roundness, 1.3); // gamma 1.3 used to balance the effect in the 0.0...0.5 roundness range
        pcv.sep = ((int)sepf) & ~0x1;
        pcv.sepmix = (sepf - pcv.sep) * 0.5; // 0.0 to 1.0
        pcv.oe1_a = powf (2.0, 1.0 / pcv.sep) * long_side * 0.5;
        pcv.oe1_b = pcv.oe1_a * short_side / long_side;
        pcv.ie1_mul = (1.0 / powf (2.0, 1.0 / pcv.sep)) * (1.0 - pcv.feather);
        pcv.oe2_a = powf (2.0, 1.0 / (pcv.sep + 2)) * long_side * 0.5;
        pcv.oe2_b = pcv.oe2_a * short_side / long_side;
        pcv.ie2_mul = (1.0 / powf (2.0, 1.0 / (pcv.sep + 2))) * (1.0 - pcv.feather);
    }

    if (roundness > 0.5) {
        // scale from fitted ellipse towards circle
        float rad = sqrtf (pcv.w * pcv.w + pcv.h * pcv.h) / 2.0;
        float diff_a = rad - pcv.oe_a;
        float diff_b = rad - pcv.oe_b;
        pcv.oe_a = pcv.oe_a + diff_a * 2 * (roundness - 0.5);
        pcv.oe_b = pcv.oe_b + diff_b * 2 * (roundness - 0.5);
    }

    pcv.scale = powf (2, -pcvignette.strength);

    if (pcvignette.strength >= 6.0) {
        pcv.scale = 0.0;
    }
}

static float calcPCVignetteFactor (const struct pcv_params& pcv, int x, int y)
{

    float fo = 1.f;

    if (x < pcv.x1 || x > pcv.x2 || y < pcv.y1 || y > pcv.y2) {
        /*
          The initial plan was to have 1.0 directly outside the crop box (ie no fading), but due to
          rounding/trunction here and there I didn't succeed matching up exactly on the pixel with
          the crop box. To hide that mismatch I made a fade.
         */
        int dist_x = (x < pcv.x1) ? pcv.x1 - x : x - pcv.x2;
        int dist_y = (y < pcv.y1) ? pcv.y1 - y : y - pcv.y2;

        if (dist_x < 0) {
            dist_x = 0;
        }

        if (dist_y < 0) {
            dist_y = 0;
        }

        fo = sqrtf (dist_x * dist_x + dist_y * dist_y) * pcv.fadeout_mul;

        if (fo >= 1.f) {
            return 1.f;
        }
    }

    float a = fabs ((x - pcv.x1) - pcv.w * 0.5f);
    float b = fabs ((y - pcv.y1) - pcv.h * 0.5f);

    if (pcv.is_portrait) {
        std::swap (a, b);
    }

    float dist = normn (a, b, 2);
    float dist_oe, dist_ie;
    float2 sincosval;

    if (dist == 0.0f) {
        sincosval.y = 1.0f;         // cos
        sincosval.x = 0.0f;         // sin
    } else {
        sincosval.y = a / dist;     // cos
        sincosval.x = b / dist;     // sin
    }

    if (pcv.is_super_ellipse_mode) {
        float dist_oe1 = pcv.oe1_a * pcv.oe1_b / normn (pcv.oe1_b * sincosval.y, pcv.oe1_a * sincosval.x, pcv.sep);
        float dist_oe2 = pcv.oe2_a * pcv.oe2_b / normn (pcv.oe2_b * sincosval.y, pcv.oe2_a * sincosval.x, pcv.sep + 2);
        float dist_ie1 = pcv.ie1_mul * dist_oe1;
        float dist_ie2 = pcv.ie2_mul * dist_oe2;
        dist_oe = dist_oe1 * (1.f - pcv.sepmix) + dist_oe2 * pcv.sepmix;
        dist_ie = dist_ie1 * (1.f - pcv.sepmix) + dist_ie2 * pcv.sepmix;
    } else {
        dist_oe = pcv.oe_a * pcv.oe_b / sqrtf (SQR (pcv.oe_b * sincosval.y) + SQR (pcv.oe_a * sincosval.x));
        dist_ie = pcv.ie_mul * dist_oe;
    }

    if (dist <= dist_ie) {
        return 1.f;
    }

    float val;

    if (dist >= dist_oe) {
        val = pcv.scale;
    } else {
        val = rtengine::RT_PI_F_2 * (dist - dist_ie) / (dist_oe - dist_ie);

        if (pcv.scale < 1.f) {
            val = pow4 (xcosf (val));
        } else {
            val = 1 - pow4 (xsinf (val));
        }

        val = pcv.scale + val * (1.f - pcv.scale);
    }

    if (fo < 1.f) {
        val = fo + val * (1.f - fo);
    }

    return val;
}

void ImProcFunctions::transformLuminanceOnly (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int oW, int oH, int fW, int fH)
{

    const bool applyVignetting = needsVignetting();
    const bool applyGradient = needsGradient();
    const bool applyPCVignetting = needsPCVignetting();

    double vig_w2, vig_h2, maxRadius, v, b, mul;

    if (applyVignetting) {
        calcVignettingParams (oW, oH, params->vignetting, vig_w2, vig_h2, maxRadius, v, b, mul);
    }

    struct grad_params gp;

    if (applyGradient) {
        calcGradientParams (oW, oH, params->gradient, gp);
    }

    struct pcv_params pcv;

    if (applyPCVignetting) {
        //fprintf(stderr, "%d %d | %d %d | %d %d | %d %d [%d %d]\n", fW, fH, oW, oH, transformed->getWidth(), transformed->getHeight(), cx, cy, params->crop.w, params->crop.h);
        calcPCVignetteParams (fW, fH, oW, oH, params->pcvignette, params->crop, pcv);
    }

    bool darkening = (params->vignetting.amount <= 0.0);
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->getHeight(); y++) {
        double vig_y_d = applyVignetting ? (double) (y + cy) - vig_h2 : 0.0;

        for (int x = 0; x < transformed->getWidth(); x++) {
            double factor = 1.0;

            if (applyVignetting) {
                double vig_x_d = (double) (x + cx) - vig_w2 ;
                double r = sqrt (vig_x_d * vig_x_d + vig_y_d * vig_y_d);

                if (darkening) {
                    factor /= std::max (v + mul * tanh (b * (maxRadius - r) / maxRadius), 0.001);
                } else {
                    factor = v + mul * tanh (b * (maxRadius - r) / maxRadius);
                }
            }

            if (applyGradient) {
                factor *= calcGradientFactor (gp, cx + x, cy + y);
            }

            if (applyPCVignetting) {
                factor *= calcPCVignetteFactor (pcv, cx + x, cy + y);
            }

            transformed->r (y, x) = original->r (y, x) * factor;
            transformed->g (y, x) = original->g (y, x) * factor;
            transformed->b (y, x) = original->b (y, x) * factor;
        }
    }
}


void ImProcFunctions::transformGeneral(ImProcFunctions::TransformMode mode, Imagefloat *original, Imagefloat *transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH, const LensCorrection *pLCPMap)
{
    double w2 = (double) oW  / 2.0 - 0.5;
    double h2 = (double) oH  / 2.0 - 0.5;

    double vig_w2, vig_h2, maxRadius, v, b, mul;
    calcVignettingParams (oW, oH, params->vignetting, vig_w2, vig_h2, maxRadius, v, b, mul);

    struct grad_params gp;

    if (needsGradient()) {
        calcGradientParams (oW, oH, params->gradient, gp);
    }

    struct pcv_params pcv;

    if (needsPCVignetting()) {
        calcPCVignetteParams (fW, fH, oW, oH, params->pcvignette, params->crop, pcv);
    }

    float** chOrig[3];
    chOrig[0] = original->r.ptrs;
    chOrig[1] = original->g.ptrs;
    chOrig[2] = original->b.ptrs;

    float** chTrans[3];
    chTrans[0] = transformed->r.ptrs;
    chTrans[1] = transformed->g.ptrs;
    chTrans[2] = transformed->b.ptrs;

    // auxiliary variables for c/a correction
    double chDist[3];
    chDist[0] = params->cacorrection.red;
    chDist[1] = 0.0;
    chDist[2] = params->cacorrection.blue;

    // auxiliary variables for distortion correction
    bool needsDist = needsDistortion();  // for performance
    double distAmount = params->distortion.amount;

    // auxiliary variables for rotation
    double cost = cos (params->rotate.degree * rtengine::RT_PI / 180.0);
    double sint = sin (params->rotate.degree * rtengine::RT_PI / 180.0);

    // auxiliary variables for vertical perspective correction
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * rtengine::RT_PI;
    double vpteta  = fabs (vpalpha - rtengine::RT_PI / 2) < 3e-4 ? 0.0 : acos ((vpdeg > 0 ? 1.0 : -1.0) * sqrt ((-SQR (oW * tan (vpalpha)) + (vpdeg > 0 ? 1.0 : -1.0) *
                     oW * tan (vpalpha) * sqrt (SQR (4 * maxRadius) + SQR (oW * tan (vpalpha)))) / (SQR (maxRadius) * 8)));
    double vpcospt = (vpdeg >= 0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

    // auxiliary variables for horizontal perspective correction
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * rtengine::RT_PI;
    double hpteta  = fabs (hpalpha - rtengine::RT_PI / 2) < 3e-4 ? 0.0 : acos ((hpdeg > 0 ? 1.0 : -1.0) * sqrt ((-SQR (oH * tan (hpalpha)) + (hpdeg > 0 ? 1.0 : -1.0) *
                     oH * tan (hpalpha) * sqrt (SQR (4 * maxRadius) + SQR (oH * tan (hpalpha)))) / (SQR (maxRadius) * 8)));
    double hpcospt = (hpdeg >= 0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

    double ascale = params->commonTrans.autofill ? getTransformAutoFill (oW, oH, pLCPMap) : 1.0;

    // smaller crop images are a problem, so only when processing fully
    bool enableLCPCA = false;
    bool enableLCPDist = false;
    bool enableCA = false;

    switch (mode) {
        case ImProcFunctions::TRANSFORM_HIGH_QUALITY_FULLIMAGE: {
            enableLCPCA = pLCPMap && params->lensProf.useCA && pLCPMap->isCACorrectionAvailable();
        }
        //no break on purpose

        case ImProcFunctions::TRANSFORM_HIGH_QUALITY: {
            enableLCPDist = pLCPMap && params->lensProf.useDist;
            if (enableLCPCA) {
                enableLCPDist = false;
            }
            enableCA = enableLCPCA || needsCA();
        }
        //no break on purpose

        default:
        case ImProcFunctions::TRANSFORM_PREVIEW: {
            enableLCPDist = pLCPMap && params->lensProf.useDist;
            break;
        }
    }

    if (!enableCA) {
        chDist[0] = 0.0;
    }

    // main cycle
    bool darkening = (params->vignetting.amount <= 0.0);
    #pragma omp parallel for if (multiThread)

    for (int y = 0; y < transformed->getHeight(); y++) {
        for (int x = 0; x < transformed->getWidth(); x++) {
            double x_d = x, y_d = y;

            if (enableLCPDist) {
                pLCPMap->correctDistortion(x_d, y_d, cx, cy, ascale); // must be first transform
            } else {
                x_d *= ascale;
                y_d *= ascale;
            }

            x_d += ascale * (cx - w2);     // centering x coord & scale
            y_d += ascale * (cy - h2);     // centering y coord & scale

            double vig_x_d = 0., vig_y_d = 0.;

            if (needsVignetting()) {
                vig_x_d = ascale * (x + cx - vig_w2);       // centering x coord & scale
                vig_y_d = ascale * (y + cy - vig_h2);       // centering y coord & scale
            }

            if (needsPerspective()) {
                // horizontal perspective transformation
                y_d *= maxRadius / (maxRadius + x_d * hptanpt);
                x_d *= maxRadius * hpcospt / (maxRadius + x_d * hptanpt);

                // vertical perspective transformation
                x_d *= maxRadius / (maxRadius - y_d * vptanpt);
                y_d *= maxRadius * vpcospt / (maxRadius - y_d * vptanpt);
            }

            // rotate
            double Dxc = x_d * cost - y_d * sint;
            double Dyc = x_d * sint + y_d * cost;

            // distortion correction
            double s = 1;

            if (needsDist) {
                double r = sqrt (Dxc * Dxc + Dyc * Dyc) / maxRadius; // sqrt is slow
                s = 1.0 - distAmount + distAmount * r ;
            }

            double r2 = 0.;

            if (needsVignetting()) {
                double vig_Dx = vig_x_d * cost - vig_y_d * sint;
                double vig_Dy = vig_x_d * sint + vig_y_d * cost;
                r2 = sqrt (vig_Dx * vig_Dx + vig_Dy * vig_Dy);
            }

            for (int c = 0; c < (enableCA ? 3 : 1); c++) {
                double Dx = Dxc * (s + chDist[c]);
                double Dy = Dyc * (s + chDist[c]);

                // de-center
                Dx += w2;
                Dy += h2;

                // LCP CA
                if (enableLCPCA) {
                    pLCPMap->correctCA (Dx, Dy, c);
                }

                // Extract integer and fractions of source screen coordinates
                int xc = (int)Dx;
                Dx -= (double)xc;
                xc -= sx;
                int yc = (int)Dy;
                Dy -= (double)yc;
                yc -= sy;

                // Convert only valid pixels
                if (yc >= 0 && yc < original->getHeight() && xc >= 0 && xc < original->getWidth()) {

                    // multiplier for vignetting correction
                    double vignmul = 1.0;

                    if (needsVignetting()) {
                        if (darkening) {
                            vignmul /= std::max (v + mul * tanh (b * (maxRadius - s * r2) / maxRadius), 0.001);
                        } else {
                            vignmul *= (v + mul * tanh (b * (maxRadius - s * r2) / maxRadius));
                        }
                    }

                    if (needsGradient()) {
                        vignmul *= calcGradientFactor (gp, cx + x, cy + y);
                    }

                    if (needsPCVignetting()) {
                        vignmul *= calcPCVignetteFactor (pcv, cx + x, cy + y);
                    }

                    if (yc > 0 && yc < original->getHeight() - 2 && xc > 0 && xc < original->getWidth() - 2) {
                        // all interpolation pixels inside image
                        if (enableCA) {
                            interpolateTransformChannelsCubic (chOrig[c], xc - 1, yc - 1, Dx, Dy, & (chTrans[c][y][x]), vignmul);
                        } else if (mode == ImProcFunctions::TRANSFORM_PREVIEW) {
                            transformed->r (y, x) = vignmul * (original->r (yc, xc) * (1.0 - Dx) * (1.0 - Dy) + original->r (yc, xc + 1) * Dx * (1.0 - Dy) + original->r (yc + 1, xc) * (1.0 - Dx) * Dy + original->r (yc + 1, xc + 1) * Dx * Dy);
                            transformed->g (y, x) = vignmul * (original->g (yc, xc) * (1.0 - Dx) * (1.0 - Dy) + original->g (yc, xc + 1) * Dx * (1.0 - Dy) + original->g (yc + 1, xc) * (1.0 - Dx) * Dy + original->g (yc + 1, xc + 1) * Dx * Dy);
                            transformed->b (y, x) = vignmul * (original->b (yc, xc) * (1.0 - Dx) * (1.0 - Dy) + original->b (yc, xc + 1) * Dx * (1.0 - Dy) + original->b (yc + 1, xc) * (1.0 - Dx) * Dy + original->b (yc + 1, xc + 1) * Dx * Dy);
                        } else {
                            interpolateTransformCubic (original, xc - 1, yc - 1, Dx, Dy, & (transformed->r (y, x)), & (transformed->g (y, x)), & (transformed->b (y, x)), vignmul);
                        }
                    } else {
                        // edge pixels
                        int y1 = LIM (yc,   0, original->getHeight() - 1);
                        int y2 = LIM (yc + 1, 0, original->getHeight() - 1);
                        int x1 = LIM (xc,   0, original->getWidth() - 1);
                        int x2 = LIM (xc + 1, 0, original->getWidth() - 1);

                        if (enableCA) {
                            chTrans[c][y][x] = vignmul * (chOrig[c][y1][x1] * (1.0 - Dx) * (1.0 - Dy) + chOrig[c][y1][x2] * Dx * (1.0 - Dy) + chOrig[c][y2][x1] * (1.0 - Dx) * Dy + chOrig[c][y2][x2] * Dx * Dy);
                        } else {
                            transformed->r (y, x) = vignmul * (original->r (y1, x1) * (1.0 - Dx) * (1.0 - Dy) + original->r (y1, x2) * Dx * (1.0 - Dy) + original->r (y2, x1) * (1.0 - Dx) * Dy + original->r (y2, x2) * Dx * Dy);
                            transformed->g (y, x) = vignmul * (original->g (y1, x1) * (1.0 - Dx) * (1.0 - Dy) + original->g (y1, x2) * Dx * (1.0 - Dy) + original->g (y2, x1) * (1.0 - Dx) * Dy + original->g (y2, x2) * Dx * Dy);
                            transformed->b (y, x) = vignmul * (original->b (y1, x1) * (1.0 - Dx) * (1.0 - Dy) + original->b (y1, x2) * Dx * (1.0 - Dy) + original->b (y2, x1) * (1.0 - Dx) * Dy + original->b (y2, x2) * Dx * Dy);
                        }
                    }
                } else {
                    if (enableCA) {
                        // not valid (source pixel x,y not inside source image, etc.)
                        chTrans[c][y][x] = 0;
                    } else {
                        transformed->r (y, x) = 0;
                        transformed->g (y, x) = 0;
                        transformed->b (y, x) = 0;
                    }
                }
            }
        }
    }
}


double ImProcFunctions::getTransformAutoFill (int oW, int oH, const LensCorrection *pLCPMap)
{
    if (!needsCA() && !needsDistortion() && !needsRotation() && !needsPerspective() && (!params->lensProf.useDist || pLCPMap == nullptr)) {
        return 1;
    }

    double scaleU = 2, scaleL = 0.001;  // upper and lower border, iterate inbetween

    do {
        double scale = (scaleU + scaleL) * 0.5;

        int orx, ory, orw, orh;
        bool clipped = transCoord (oW, oH, 0, 0, oW, oH, orx, ory, orw, orh, scale, pLCPMap);

        if (clipped) {
            scaleU = scale;
        } else {
            scaleL = scale;
        }
    } while (scaleU - scaleL > 0.001);

    return scaleL;
}

bool ImProcFunctions::needsCA ()
{
    return fabs (params->cacorrection.red) > 1e-15 || fabs (params->cacorrection.blue) > 1e-15;
}

bool ImProcFunctions::needsDistortion ()
{
    return fabs (params->distortion.amount) > 1e-15;
}

bool ImProcFunctions::needsRotation ()
{
    return fabs (params->rotate.degree) > 1e-15;
}

bool ImProcFunctions::needsPerspective ()
{
    return params->perspective.horizontal || params->perspective.vertical;
}

bool ImProcFunctions::needsGradient ()
{
    return params->gradient.enabled && fabs (params->gradient.strength) > 1e-15;
}

bool ImProcFunctions::needsPCVignetting ()
{
    return params->pcvignette.enabled && fabs (params->pcvignette.strength) > 1e-15;
}

bool ImProcFunctions::needsVignetting ()
{
    return params->vignetting.amount;
}

bool ImProcFunctions::needsLCP ()
{
    return params->lensProf.useLcp();
}

bool ImProcFunctions::needsLensfun()
{
    return params->lensProf.useLensfun();
}

bool ImProcFunctions::needsTransform ()
{
    return needsCA () || needsDistortion () || needsRotation () || needsPerspective () || needsGradient () || needsPCVignetting () || needsVignetting () || needsLCP() || needsLensfun();
}


}

