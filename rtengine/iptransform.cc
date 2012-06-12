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

using namespace std;

namespace rtengine {
#undef CLIPTOC

#define CLIPTOC(a,b,c,d) ((a)>=(b)?((a)<=(c)?(a):(d=true,(c))):(d=true,(b)))
#define RT_PI 3.141592653589

bool ImProcFunctions::transCoord (int W, int H, const std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef,
    const LCPMapper *pLCPMap) {

    bool clipresize = true, clipped = false;

    red.clear (); green.clear (); blue.clear ();

    if (!needsCA() && !needsDistortion() && !needsRotation() && !needsPerspective() && (!params->lensProf.useDist || pLCPMap==NULL)) {
        if (clipresize) {
                for (size_t i=0; i<src.size(); i++) {
                    red.push_back   (Coord2D (src[i].x, src[i].y));
                    green.push_back (Coord2D (src[i].x, src[i].y));
                    blue.push_back  (Coord2D (src[i].x, src[i].y));
                }
        }
        return clipped;
    }

    double oW = W, oH = H;
	double w2 = (double) oW  / 2.0 - 0.5;
	double h2 = (double) oH  / 2.0 - 0.5;
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;

    // auxiliary variables for distortion correction
    bool needsDist = needsDistortion();  // for performance
	double distAmount = params->distortion.amount;

    // auxiliary variables for rotation
	double cost = cos(params->rotate.degree * RT_PI/180.0);
	double sint = sin(params->rotate.degree * RT_PI/180.0);

    // auxiliary variables for vertical perspective correction
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * RT_PI;
    double vpteta  = fabs(vpalpha-RT_PI/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

    // auxiliary variables for horizontal perspective correction
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * RT_PI;
    double hpteta  = fabs(hpalpha-RT_PI/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

	double ascale = ascaleDef>0 ? ascaleDef : (params->commonTrans.autofill ? getTransformAutoFill (oW, oH, pLCPMap) : 1.0);

	for (size_t i=0; i<src.size(); i++) {
        double x_d=src[i].x, y_d=src[i].y;
        if (pLCPMap && params->lensProf.useDist) pLCPMap->correctDistortion(x_d,y_d);  // must be first transform

		y_d = ascale * (y_d - h2);
		x_d = ascale * (x_d - w2);

        if (needsPerspective()) {
            // horizontal perspective transformation
            y_d *= maxRadius / (maxRadius + x_d*hptanpt);
            x_d *= maxRadius * hpcospt / (maxRadius + x_d*hptanpt);

            // vertical perspective transformation
            x_d *= maxRadius / (maxRadius - y_d*vptanpt);
            y_d *= maxRadius * vpcospt / (maxRadius - y_d*vptanpt);
        }

        // rotate
		double Dx = x_d * cost - y_d * sint;
		double Dy = x_d * sint + y_d * cost;

        // distortion correction
        double s = 1;
        if (needsDist) {
            double r = sqrt(Dx*Dx + Dy*Dy) / maxRadius;  // sqrt is slow
            s = 1.0 - distAmount + distAmount * r ;
        }

        // LCP CA is not reflected in preview (and very small), so don't add it here

		red.push_back (Coord2D(Dx*(s+params->cacorrection.red)+w2, Dy*(s+params->cacorrection.red)+h2));
		green.push_back (Coord2D(Dx*s+w2, Dy*s+h2));
		blue.push_back (Coord2D(Dx*(s+params->cacorrection.blue)+w2, Dy*(s+params->cacorrection.blue)+h2));
	}

	if (clipresize) {
        for (size_t i=0; i<src.size(); i++) {
            red[i].x = CLIPTOC(red[i].x,0,W-1,clipped);
            red[i].y = CLIPTOC(red[i].y,0,H-1,clipped);
            green[i].x = CLIPTOC(green[i].x,0,W-1,clipped);
            green[i].y = CLIPTOC(green[i].y,0,H-1,clipped);
            blue[i].x = CLIPTOC(blue[i].x,0,W-1,clipped);
            blue[i].y = CLIPTOC(blue[i].y,0,H-1,clipped);
        }
    }
    return clipped;
}

// Transform all corners and critical sidelines of an image
bool ImProcFunctions::transCoord (int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef, const LCPMapper *pLCPMap) {

    int x1 = x, y1 = y;
    int x2 = x1 + w - 1;
    int y2 = y1 + h - 1;

    // Build all edge points and half-way points
    std::vector<Coord2D> corners (8);
    corners[0].set (x1, y1);
    corners[1].set (x1, y2);
    corners[2].set (x2, y2);
    corners[3].set (x2, y1);
    corners[4].set ((x1+x2)/2, y1);
    corners[5].set ((x1+x2)/2, y2);
    corners[6].set (x1, (y1+y2)/2);
    corners[7].set (x2, (y1+y2)/2);

    // Add several steps inbetween
    int xstep = (x2-x1)/32;
    if (xstep<1) xstep = 1;
    for (int i=x1+xstep; i<=x2-xstep; i+=xstep) {
        corners.push_back (Coord2D (i, y1));
        corners.push_back (Coord2D (i, y2));
    }
    int ystep = (y2-y1)/32;
    if (ystep<1) ystep = 1;
    for (int i=y1+ystep; i<=y2-ystep; i+=ystep) {
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
    for (size_t i=1; i<transCorners.size(); i++)
        if (transCorners[i].x<x1d)
            x1d = transCorners[i].x;
   int x1v = (int)(x1d);

    double y1d = transCorners[0].y;
    for (size_t i=1; i<transCorners.size(); i++)
        if (transCorners[i].y<y1d)
            y1d = transCorners[i].y;
    int y1v = (int)(y1d);

    double x2d = transCorners[0].x;
    for (size_t i=1; i<transCorners.size(); i++)
        if (transCorners[i].x>x2d)
            x2d = transCorners[i].x;
    int x2v = (int)ceil(x2d);

    double y2d = transCorners[0].y;
    for (size_t i=1; i<transCorners.size(); i++)
        if (transCorners[i].y>y2d)
            y2d = transCorners[i].y;
    int y2v = (int)ceil(y2d);

    xv = x1v;
    yv = y1v;
    wv = x2v - x1v + 1;
    hv = y2v - y1v + 1;

    return clipped;
}

void ImProcFunctions::transform (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, 
    double focalLen, double focalLen35mm, float focusDist, int rawRotationDeg, bool fullImage) {

    LCPMapper *pLCPMap=NULL;
    if (needsLCP() && focalLen>0) {
        LCPProfile *pLCPProf=lcpStore->getProfile(params->lensProf.lcpFile);
        if (pLCPProf) pLCPMap=new LCPMapper(pLCPProf, focalLen, focalLen35mm, focusDist, 0, false, params->lensProf.useDist,
            original->width, original->height, params->coarse, rawRotationDeg);
    }

	if (!(needsCA() || needsDistortion() || needsRotation() || needsPerspective() || needsLCP()) && needsVignetting())
		transformVignetteOnly (original, transformed, cx, cy, oW, oH);
	else if (!needsCA() && scale!=1)
		transformPreview (original, transformed, cx, cy, sx, sy, oW, oH, pLCPMap);
		else
		transformHighQuality (original, transformed, cx, cy, sx, sy, oW, oH, pLCPMap, fullImage);

    if (pLCPMap) delete pLCPMap;
}

// helper function
void ImProcFunctions::calcVignettingParams(int oW, int oH, const VignettingParams& vignetting, double &w2, double &h2, double& maxRadius, double &v, double &b, double &mul)
{
	// vignette center is a point with coordinates between -1 and +1
	double x = vignetting.centerX / 100.0;
	double y = vignetting.centerY / 100.0;

	// calculate vignette center in pixels 
	w2 = (double) oW  / 2.0 - 0.5 + x * oW;
	h2 = (double) oH  / 2.0 - 0.5 + y * oH;

	// max vignette radius in pixels
	maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2.;

	// vignette variables with applied strength
	v = 1.0 - vignetting.strength * vignetting.amount * 3.0 / 400.0;
	b = 1.0 + vignetting.radius * 7.0 / 100.0;
	mul = (1.0-v) / tanh(b);
}

// Transform vignetting only
void ImProcFunctions::transformVignetteOnly (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int oW, int oH) {

	double vig_w2, vig_h2, maxRadius, v, b, mul;
	calcVignettingParams(oW, oH, params->vignetting, vig_w2, vig_h2, maxRadius, v, b, mul);

	#pragma omp parallel for if (multiThread)
	for (int y=0; y<transformed->height; y++) {
		double vig_y_d = (double) (y + cy) - vig_h2 ;
		for (int x=0; x<transformed->width; x++) {
			double vig_x_d = (double) (x + cx) - vig_w2 ;
			double r = sqrt(vig_x_d*vig_x_d + vig_y_d*vig_y_d);
			double vign = v + mul * tanh (b*(maxRadius-r) / maxRadius);
			transformed->r[y][x] = original->r[y][x] / vign;
			transformed->g[y][x] = original->g[y][x] / vign;
			transformed->b[y][x] = original->b[y][x] / vign;
		}
	}
}

// Transform WITH scaling (opt.) and CA, cubic interpolation
#include "cubintch.cc"
#include "cubint.cc"

void ImProcFunctions::transformHighQuality (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, 
    const LCPMapper *pLCPMap, bool fullImage) {

	double w2 = (double) oW  / 2.0 - 0.5;
	double h2 = (double) oH  / 2.0 - 0.5;

	double vig_w2,vig_h2,maxRadius,v,b,mul;
	calcVignettingParams(oW, oH, params->vignetting, vig_w2, vig_h2, maxRadius, v, b, mul);

    float** chOrig[3];
    chOrig[0] = original->r;
    chOrig[1] = original->g;
    chOrig[2] = original->b;

    float** chTrans[3];
    chTrans[0] = transformed->r;
    chTrans[1] = transformed->g;
    chTrans[2] = transformed->b;

    // auxiliary variables for c/a correction
    double chDist[3];
    chDist[0] = params->cacorrection.red;
    chDist[1] = 0.0;
    chDist[2] = params->cacorrection.blue;

	// auxiliary variables for distortion correction
    bool needsDist = needsDistortion();  // for performance
	double distAmount = params->distortion.amount;

	// auxiliary variables for rotation
	double cost = cos(params->rotate.degree * RT_PI/180.0);
	double sint = sin(params->rotate.degree * RT_PI/180.0);

	// auxiliary variables for vertical perspective correction
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * RT_PI;
    double vpteta  = fabs(vpalpha-RT_PI/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-SQR(oW*tan(vpalpha)) + (vpdeg>0 ? 1.0 : -1.0) *
																oW*tan(vpalpha)*sqrt(SQR(4*maxRadius)+SQR(oW*tan(vpalpha))))/(SQR(maxRadius)*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

	// auxiliary variables for horizontal perspective correction
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * RT_PI;
    double hpteta  = fabs(hpalpha-RT_PI/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-SQR(oH*tan(hpalpha)) + (hpdeg>0 ? 1.0 : -1.0) *
																oH*tan(hpalpha)*sqrt(SQR(4*maxRadius)+SQR(oH*tan(hpalpha))))/(SQR(maxRadius)*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

	double ascale = params->commonTrans.autofill ? getTransformAutoFill (oW, oH, fullImage ? pLCPMap : NULL) : 1.0;

    // smaller crop images are a problem, so only when processing fully
            bool enableLCPCA   = pLCPMap && params->lensProf.useCA && fullImage && pLCPMap->enableCA; 
    bool enableLCPDist = pLCPMap && params->lensProf.useDist && fullImage;
    if (enableLCPCA) enableLCPDist=false;
            bool enableCA = enableLCPCA || needsCA();

	// main cycle
	#pragma omp parallel for if (multiThread)
    for (int y=0; y<transformed->height; y++) {
        for (int x=0; x<transformed->width; x++) {
            double x_d=x,y_d=y;
            if (enableLCPDist) pLCPMap->correctDistortion(x_d,y_d);  // must be first transform

            x_d = ascale * (x_d + cx - w2);		// centering x coord & scale
            y_d = ascale * (y_d + cy - h2);		// centering y coord & scale

            double vig_x_d, vig_y_d;
            if (needsVignetting()) {
                vig_x_d = ascale * (x + cx - vig_w2);		// centering x coord & scale
                vig_y_d = ascale * (y + cy - vig_h2);		// centering y coord & scale
            }

            if (needsPerspective()) {
            // horizontal perspective transformation
                y_d *= maxRadius / (maxRadius + x_d*hptanpt);
                x_d *= maxRadius * hpcospt / (maxRadius + x_d*hptanpt);

            // vertical perspective transformation
                x_d *= maxRadius / (maxRadius - y_d*vptanpt);
                y_d *= maxRadius * vpcospt / (maxRadius - y_d*vptanpt);
            }

            // rotate
            double Dxc = x_d * cost - y_d * sint;
            double Dyc = x_d * sint + y_d * cost;

            // distortion correction
            double s = 1;
            if (needsDist) {
                double r = sqrt(Dxc*Dxc + Dyc*Dyc) / maxRadius;  // sqrt is slow
                s = 1.0 - distAmount + distAmount * r ;
            }

            double r2;
            if (needsVignetting()) {  
            double vig_Dx = vig_x_d * cost - vig_y_d * sint;
            double vig_Dy = vig_x_d * sint + vig_y_d * cost;
                r2=sqrt(vig_Dx*vig_Dx + vig_Dy*vig_Dy);
            }

                    for (int c=0; c < (enableCA ? 3 : 1); c++) {
                double Dx = Dxc * (s + chDist[c]);
                double Dy = Dyc * (s + chDist[c]);

                // de-center
                Dx += w2; Dy += h2;

                // LCP CA
                if (enableLCPCA) pLCPMap->correctCA(Dx,Dy,c);

				// Extract integer and fractions of source screen coordinates
				int xc = (int)Dx; Dx -= (double)xc; xc -= sx;
				int yc = (int)Dy; Dy -= (double)yc; yc -= sy;

				// Convert only valid pixels
				if (yc>=0 && yc<original->height && xc>=0 && xc<original->width) {

					// multiplier for vignetting correction
					double vignmul = 1.0;
					if (needsVignetting())
						vignmul /= (v + mul * tanh (b*(maxRadius-s*r2) / maxRadius));

					if (yc > 0 && yc < original->height-2 && xc > 0 && xc < original->width-2) {
                        // all interpolation pixels inside image
                                if (enableCA) 
                        interpolateTransformChannelsCubic (chOrig[c], xc-1, yc-1, Dx, Dy, &(chTrans[c][y][x]), vignmul);
                                else
                                    interpolateTransformCubic (original, xc-1, yc-1, Dx, Dy, &(transformed->r[y][x]), &(transformed->g[y][x]), &(transformed->b[y][x]), vignmul);
                    } else { 
                        // edge pixels
						int y1 = LIM(yc,   0, original->height-1);
						int y2 = LIM(yc+1, 0, original->height-1);
						int x1 = LIM(xc,   0, original->width-1);
						int x2 = LIM(xc+1, 0, original->width-1);

                                if (enableCA) {
                        chTrans[c][y][x] = vignmul * (chOrig[c][y1][x1]*(1.0-Dx)*(1.0-Dy) + chOrig[c][y1][x2]*Dx*(1.0-Dy) + chOrig[c][y2][x1]*(1.0-Dx)*Dy + chOrig[c][y2][x2]*Dx*Dy);
                                } else {
                                    transformed->r[y][x] = vignmul*(original->r[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->r[y1][x2]*Dx*(1.0-Dy) + original->r[y2][x1]*(1.0-Dx)*Dy + original->r[y2][x2]*Dx*Dy);
                                    transformed->g[y][x] = vignmul*(original->g[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->g[y1][x2]*Dx*(1.0-Dy) + original->g[y2][x1]*(1.0-Dx)*Dy + original->g[y2][x2]*Dx*Dy);
                                    transformed->b[y][x] = vignmul*(original->b[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->b[y1][x2]*Dx*(1.0-Dy) + original->b[y2][x1]*(1.0-Dx)*Dy + original->b[y2][x2]*Dx*Dy);
					}
				}
                        }
                        else {
                            if (enableCA) {
					// not valid (source pixel x,y not inside source image, etc.)
					chTrans[c][y][x] = 0;
                            } else {
                                transformed->r[y][x] = 0;
                                transformed->g[y][x] = 0;
                                transformed->b[y][x] = 0;
                            }
                        }
			}
        }
    }
}

// Transform WITH scaling, WITHOUT CA, simple (and fast) interpolation. Used for preview
void ImProcFunctions::transformPreview (Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, const LCPMapper *pLCPMap) {

	double w2 = (double) oW  / 2.0 - 0.5;
	double h2 = (double) oH  / 2.0 - 0.5;

	double vig_w2, vig_h2, maxRadius, v, b, mul;
	calcVignettingParams(oW, oH, params->vignetting, vig_w2, vig_h2, maxRadius, v, b, mul);

	// auxiliary variables for distortion correction
    bool needsDist = needsDistortion();  // for performance
	double distAmount = params->distortion.amount;

	// auxiliary variables for rotation
	double cost = cos(params->rotate.degree * RT_PI/180.0);
	double sint = sin(params->rotate.degree * RT_PI/180.0);

	// auxiliary variables for vertical perspective correction
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90 - vpdeg) / 180.0 * RT_PI;
    double vpteta  = fabs(vpalpha-RT_PI/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

	// auxiliary variables for horizontal perspective correction
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90 - hpdeg) / 180.0 * RT_PI;
    double hpteta  = fabs(hpalpha-RT_PI/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

	double ascale = params->commonTrans.autofill ? getTransformAutoFill (oW, oH, pLCPMap) : 1.0;

    // main cycle
	#pragma omp parallel for if (multiThread)
    for (int y=0; y<transformed->height; y++) {
        for (int x=0; x<transformed->width; x++) {
            double x_d=x,y_d=y;
            if (pLCPMap && params->lensProf.useDist) pLCPMap->correctDistortion(x_d,y_d);  // must be first transform

            y_d = ascale * (y_d + cy - h2);		// centering y coord & scale
            x_d = ascale * (x_d + cx - w2);		// centering x coord & scale

            double vig_x_d, vig_y_d;
            if (needsVignetting()) {
                vig_x_d = ascale * (x + cx - vig_w2);		// centering x coord & scale
                vig_y_d = ascale * (y + cy - vig_h2);		// centering y coord & scale
            }

            if (needsPerspective()) {
            // horizontal perspective transformation
                y_d *= maxRadius / (maxRadius + x_d*hptanpt);
                x_d *= maxRadius * hpcospt / (maxRadius + x_d*hptanpt);

            // vertical perspective transformation
                x_d *= maxRadius / (maxRadius - y_d*vptanpt);
                y_d *= maxRadius * vpcospt / (maxRadius - y_d*vptanpt);
            }

			// rotate
            double Dx = x_d * cost - y_d * sint;
            double Dy = x_d * sint + y_d * cost;

            // distortion correction
            double s = 1;
            if (needsDist) {
                double r = sqrt(Dx*Dx + Dy*Dy) / maxRadius;  // sqrt is slow
                s = 1.0 - distAmount + distAmount * r ;
	        Dx *= s;
            Dy *= s;
            }

            double r2;
            if (needsVignetting()) {  
            double vig_Dx = vig_x_d * cost - vig_y_d * sint;
            double vig_Dy = vig_x_d * sint + vig_y_d * cost;
                r2=sqrt(vig_Dx*vig_Dx + vig_Dy*vig_Dy);
            }

            // de-center
            Dx += w2; Dy += h2;

            // Extract integer and fractions of source screen coordinates
            int xc = (int)Dx; Dx -= (double)xc; xc -= sx;
            int yc = (int)Dy; Dy -= (double)yc; yc -= sy;

            // Convert only valid pixels
            if (yc>=0 && yc<original->height && xc>=0 && xc<original->width) {

                // multiplier for vignetting correction
            	double vignmul = 1.0;
                if (needsVignetting())
                	vignmul /= (v + mul * tanh (b*(maxRadius-s*r2) / maxRadius));

                if (yc < original->height-1 && xc < original->width-1) {  
                    // all interpolation pixels inside image
                    transformed->r[y][x] = vignmul*(original->r[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->r[yc][xc+1]*Dx*(1.0-Dy) + original->r[yc+1][xc]*(1.0-Dx)*Dy + original->r[yc+1][xc+1]*Dx*Dy);
                    transformed->g[y][x] = vignmul*(original->g[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->g[yc][xc+1]*Dx*(1.0-Dy) + original->g[yc+1][xc]*(1.0-Dx)*Dy + original->g[yc+1][xc+1]*Dx*Dy);
                    transformed->b[y][x] = vignmul*(original->b[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->b[yc][xc+1]*Dx*(1.0-Dy) + original->b[yc+1][xc]*(1.0-Dx)*Dy + original->b[yc+1][xc+1]*Dx*Dy);
                }
                else { 
                    // edge pixels
                	int y1 = LIM(yc,   0, original->height-1);
                	int y2 = LIM(yc+1, 0, original->height-1);
                	int x1 = LIM(xc,   0, original->width-1);
                	int x2 = LIM(xc+1, 0, original->width-1);
                    transformed->r[y][x] = vignmul*(original->r[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->r[y1][x2]*Dx*(1.0-Dy) + original->r[y2][x1]*(1.0-Dx)*Dy + original->r[y2][x2]*Dx*Dy);
                    transformed->g[y][x] = vignmul*(original->g[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->g[y1][x2]*Dx*(1.0-Dy) + original->g[y2][x1]*(1.0-Dx)*Dy + original->g[y2][x2]*Dx*Dy);
                    transformed->b[y][x] = vignmul*(original->b[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->b[y1][x2]*Dx*(1.0-Dy) + original->b[y2][x1]*(1.0-Dx)*Dy + original->b[y2][x2]*Dx*Dy);
                }
            }
            else {
                // not valid (source pixel x,y not inside source image, etc.)
                transformed->r[y][x] = 0;
                transformed->g[y][x] = 0;
                transformed->b[y][x] = 0;
            }
        }
    }
}

double ImProcFunctions::getTransformAutoFill (int oW, int oH, const LCPMapper *pLCPMap) {

	double scaleU = 1.0;
	double scaleL = 0.001;
	while (scaleU - scaleL > 0.001) {
		double scale = (scaleU + scaleL) / 2.0;

        int orx, ory, orw, orh;
        bool clipped = transCoord (oW, oH, 0, 0, oW, oH, orx, ory, orw, orh, scale, pLCPMap);

        if (clipped)
        	scaleU = scale;
        else
        	scaleL = scale;
	}
	return scaleL;
}

bool ImProcFunctions::needsCA () {
	return fabs (params->cacorrection.red) > 1e-15 || fabs (params->cacorrection.blue) > 1e-15;
}

bool ImProcFunctions::needsDistortion () {
	return fabs (params->distortion.amount) > 1e-15;
}

bool ImProcFunctions::needsRotation	() {
	return fabs (params->rotate.degree) > 1e-15;
}

bool ImProcFunctions::needsPerspective () {
	return params->perspective.horizontal || params->perspective.vertical;
}

bool ImProcFunctions::needsVignetting () {
	return params->vignetting.amount;
}

bool ImProcFunctions::needsLCP () {
	return params->lensProf.lcpFile.length()>0;
}

bool ImProcFunctions::needsTransform () {
	return needsCA () || needsDistortion () || needsRotation () || needsPerspective () || needsVignetting () || needsLCP();
}


}

