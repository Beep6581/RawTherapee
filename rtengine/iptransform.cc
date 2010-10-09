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
#include <rtengine.h>
#include <improcfun.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <mytime.h>

namespace rtengine {

#undef CMAXVAL
#undef MAX
#undef MIN
#undef CLIP
#undef CLIPTO
#undef CLIPTOC

#define CMAXVAL 0xffff
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIPTOC(a,b,c,d) ((a)>=(b)?((a)<=(c)?(a):(d=true,(c))):(d=true,(b)))

bool ImProcFunctions::transCoord (int W, int H, std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef) {

    bool clipresize = true;
    bool clipped = false;

    red.clear ();
    green.clear ();
    blue.clear ();

    if (!needsCA() && !needsDistortion() && !needsRotation() && !needsPerspective()) {
        if (clipresize) {
            // Apply resizing
            if (fabs(params->resize.scale-1.0)>=1e-7) {
                for (int i=0; i<src.size(); i++) {
                    red.push_back   (Coord2D (src[i].x / params->resize.scale, src[i].y / params->resize.scale));
                    green.push_back (Coord2D (src[i].x / params->resize.scale, src[i].y / params->resize.scale));
                    blue.push_back  (Coord2D (src[i].x / params->resize.scale, src[i].y / params->resize.scale));
                }
                for (int i=0; i<src.size(); i++) {
                    red[i].x = CLIPTOC(red[i].x,0,W-1,clipped);
                    red[i].y = CLIPTOC(red[i].y,0,H-1,clipped);
                    green[i].x = CLIPTOC(green[i].x,0,W-1,clipped);
                    green[i].y = CLIPTOC(green[i].y,0,H-1,clipped);
                    blue[i].x = CLIPTOC(blue[i].x,0,W-1,clipped);
                    blue[i].y = CLIPTOC(blue[i].y,0,H-1,clipped);
                }
            }
            else
                for (int i=0; i<src.size(); i++) {
                    red.push_back   (Coord2D (src[i].x, src[i].y));
                    green.push_back (Coord2D (src[i].x, src[i].y));
                    blue.push_back  (Coord2D (src[i].x, src[i].y));
                }
        }
        return clipped;
    }

    double oW = W*params->resize.scale;
    double oH = H*params->resize.scale;
	double w2 = (double) oW  / 2.0 - 0.5;
	double h2 = (double) oH  / 2.0 - 0.5;
	double a = params->distortion.amount;
	double cost = cos(params->rotate.degree * 3.14/180.0);
	double sint = sin(params->rotate.degree * 3.14/180.0);
	double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * 3.14;
    double vpteta  = fabs(vpalpha-3.14/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * 3.14;
    double hpteta  = fabs(hpalpha-3.14/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

	double ascale = ascaleDef>0 ? ascaleDef : (params->commonTrans.autofill ? getTransformAutoFill (oW, oH) : 1.0);

	for (int i=0; i<src.size(); i++) {

		double y_d = ascale * (src[i].y - h2);
		double x_d = ascale * (src[i].x - w2);

        y_d = y_d * maxRadius / (maxRadius + x_d*hptanpt);
        x_d = x_d * maxRadius * hpcospt / (maxRadius + x_d*hptanpt);

        x_d = x_d * maxRadius / (maxRadius - y_d*vptanpt);
        y_d = y_d * maxRadius * vpcospt / (maxRadius - y_d*vptanpt);

		double Dx = x_d * cost - y_d * sint;
		double Dy = x_d * sint + y_d * cost;

		double r = sqrt(Dx*Dx + Dy*Dy) / maxRadius;
		double s = 1.0 - a + a * r ;

		red.push_back (Coord2D(Dx*(s+params->cacorrection.red)+w2, Dy*(s+params->cacorrection.red)+h2));
		green.push_back (Coord2D(Dx*s+w2, Dy*s+h2));
		blue.push_back (Coord2D(Dx*(s+params->cacorrection.blue)+w2, Dy*(s+params->cacorrection.blue)+h2));
	}

	if (clipresize) {
        if (fabs(params->resize.scale-1.0)>=1e-7) {
            for (int i=0; i<src.size(); i++) {
                red[i].x /= params->resize.scale;
                red[i].y /= params->resize.scale;
                green[i].x /= params->resize.scale;
                green[i].y /= params->resize.scale;
                blue[i].x /= params->resize.scale;
                blue[i].y /= params->resize.scale;
            }
        }
        for (int i=0; i<src.size(); i++) {
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

bool ImProcFunctions::transCoord (int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef) {

    int x1 = x, y1 = y;
    int x2 = x1 + w - 1;
    int y2 = y1 + h - 1;

    std::vector<Coord2D> corners (8);
    corners[0].set (x1, y1);
    corners[1].set (x1, y2);
    corners[2].set (x2, y2);
    corners[3].set (x2, y1);
    corners[4].set ((x1+x2)/2, y1);
    corners[5].set ((x1+x2)/2, y2);
    corners[6].set (x1, (y1+y2)/2);
    corners[7].set (x2, (y1+y2)/2);
    int xstep = (x2-x1)/16;
    if (xstep<1) xstep = 1;
    for (int i=x1+xstep; i<=x2-xstep; i+=xstep) {
        corners.push_back (Coord2D (i, y1));
        corners.push_back (Coord2D (i, y2));
    }
    int ystep = (y2-y1)/16;
    if (ystep<1) ystep = 1;
    for (int i=y1+ystep; i<=y2-ystep; i+=ystep) {
        corners.push_back (Coord2D (x1, i));
        corners.push_back (Coord2D (x2, i));
    }

    std::vector<Coord2D> r, g, b;

    bool result = transCoord (W, H, corners, r, g, b, ascaleDef);

    std::vector<Coord2D> transCorners;
    transCorners.insert (transCorners.end(), r.begin(), r.end());
    transCorners.insert (transCorners.end(), g.begin(), g.end());
    transCorners.insert (transCorners.end(), b.begin(), b.end());

    double x1d = transCorners[0].x;
    for (int i=1; i<transCorners.size(); i++)
        if (transCorners[i].x<x1d)
            x1d = transCorners[i].x;
   int x1v = (int)(x1d);

    double y1d = transCorners[0].y;
    for (int i=1; i<transCorners.size(); i++)
        if (transCorners[i].y<y1d)
            y1d = transCorners[i].y;
    int y1v = (int)(y1d);

    double x2d = transCorners[0].x;
    for (int i=1; i<transCorners.size(); i++)
        if (transCorners[i].x>x2d)
            x2d = transCorners[i].x;
    int x2v = (int)ceil(x2d);

    double y2d = transCorners[0].y;
    for (int i=1; i<transCorners.size(); i++)
        if (transCorners[i].y>y2d)
            y2d = transCorners[i].y;
    int y2v = (int)ceil(y2d);

    xv = x1v;
    yv = y1v;
    wv = x2v - x1v + 1;
    hv = y2v - y1v + 1;

    return result;
}

void ImProcFunctions::transform (Image16* original, Image16* transformed, int cx, int cy, int sx, int sy, int oW, int oH) {

	if (!(needsCA() || needsDistortion() || needsRotation() || needsPerspective()) && needsVignetting())
		vignetting (original, transformed, cx, cy, oW, oH);
	else if (!needsCA()) {
		if (scale==1)
			transformNonSep (original, transformed, cx, cy, sx, sy, oW, oH);
		else
			simpltransform (original, transformed, cx, cy, sx, sy, oW, oH);
	}
	else
		transformSep (original, transformed, cx, cy, sx, sy, oW, oH);
}

void ImProcFunctions::vignetting (Image16* original, Image16* transformed, int cx, int cy, int oW, int oH) {

	double  w2 = (double) oW  / 2.0 - 0.5;
	double  h2 = (double) oH  / 2.0 - 0.5;

	double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2.;

	double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
	double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;

	double mul = (1.0-v) / tanh(b);

	#pragma omp parallel for if (multiThread)
	for (int y=0; y<transformed->height; y++) {
		double y_d = (double) (y + cy) - h2 ;
		int val;
		for (int x=0; x<transformed->width; x++) {
			double x_d = (double) (x + cx) - w2 ;
			double r = sqrt(x_d*x_d + y_d*y_d);
			double vign = v + mul * tanh (b*(maxRadius-r) / maxRadius);
			val = original->r[y][x] / vign;
			transformed->r[y][x] = CLIP(val);
			val =  original->g[y][x] / vign;
			transformed->g[y][x] = CLIP(val);
			val = original->b[y][x] / vign;
			transformed->b[y][x] = CLIP(val);
		}
	}
}

#include "cubint.cc"
void ImProcFunctions::transformNonSep (Image16* original, Image16* transformed, int cx, int cy, int sx, int sy, int oW, int oH) {

	double  w2 = (double) oW  / 2.0 - 0.5;
	double  h2 = (double) oH  / 2.0 - 0.5;

	// auxiliary variables for distortion correction
	double a = params->distortion.amount;

	// auxiliary variables for rotation
	double cost = cos(params->rotate.degree * 3.14/180.0);
	double sint = sin(params->rotate.degree * 3.14/180.0);

	// auxiliary variables for vignetting
	double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2.;
	double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
	double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;
	double mul = (1.0-v) / tanh(b);
	bool dovign = params->vignetting.amount != 0;

	// auxiliary variables for vertical perspective correction
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * 3.14;
    double vpteta  = fabs(vpalpha-3.14/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

	// auxiliary variables for horizontal perspective correction
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * 3.14;
    double hpteta  = fabs(hpalpha-3.14/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

	double ascale = params->commonTrans.autofill ? getTransformAutoFill (oW, oH) : 1.0;

	// main cycle
	#pragma omp parallel for if (multiThread)
    for (int y=0; y<transformed->height; y++) {
        for (int x=0; x<transformed->width; x++) {
            double x_d = ascale * (x + cx - w2);		// centering x coord & scale
            double y_d = ascale * (y + cy - h2);		// centering y coord & scale

            // horizontal perspective transformation
            y_d = y_d * maxRadius / (maxRadius + x_d*hptanpt);
            x_d = x_d * maxRadius * hpcospt / (maxRadius + x_d*hptanpt);

            // vertical perspective transformation
            x_d = x_d * maxRadius / (maxRadius - y_d*vptanpt);
            y_d = y_d * maxRadius * vpcospt / (maxRadius - y_d*vptanpt);

			// rotate
            double Dx = x_d * cost - y_d * sint;
            double Dy = x_d * sint + y_d * cost;

            // distortion correction
            double r = sqrt(Dx*Dx + Dy*Dy);
            double s = 1.0 - a + a * r ;
	        Dx *= s;
            Dy *= s;

            // de-center
            Dx += w2;
            Dy += h2;

            // Extract integer and fractions of source screen coordinates
            int xc = (int)Dx; Dx -= (double)xc; xc -= sx;
            int yc = (int)Dy; Dy -= (double)yc; yc -= sy;

            // Convert only valid pixels
            if (yc>=0 && yc<original->height && xc>=0 && xc<original->width) {

                // multiplier for vignetting correction
            	double vignmul = 1.0;
                if (dovign)
                	vignmul /= (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                if (yc > 0 && yc < original->height-2 && xc > 0 && xc < original->width-2)   // all interpolation pixels inside image
                    cubint (original, xc-1, yc-1, Dx, Dy, &(transformed->r[y][x]), &(transformed->g[y][x]), &(transformed->b[y][x]), vignmul);
                else { // edge pixels
                	int y1 = CLIPTO(yc,   0, original->height-1);
                	int y2 = CLIPTO(yc+1, 0, original->height-1);
                	int x1 = CLIPTO(xc,   0, original->width-1);
                	int x2 = CLIPTO(xc+1, 0, original->width-1);
                    int r = vignmul*(original->r[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->r[y1][x2]*Dx*(1.0-Dy) + original->r[y2][x1]*(1.0-Dx)*Dy + original->r[y2][x2]*Dx*Dy);
                    int g = vignmul*(original->g[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->g[y1][x2]*Dx*(1.0-Dy) + original->g[y2][x1]*(1.0-Dx)*Dy + original->g[y2][x2]*Dx*Dy);
                    int b = vignmul*(original->b[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->b[y1][x2]*Dx*(1.0-Dy) + original->b[y2][x1]*(1.0-Dx)*Dy + original->b[y2][x2]*Dx*Dy);
                    transformed->r[y][x] = CLIP(r);
                    transformed->g[y][x] = CLIP(g);
                    transformed->b[y][x] = CLIP(b);
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

#include "cubintch.cc"
void ImProcFunctions::transformSep (Image16* original, Image16* transformed, int cx, int cy, int sx, int sy, int oW, int oH) {

	double  w2 = (double) oW  / 2.0 - 0.5;
	double  h2 = (double) oH  / 2.0 - 0.5;

	// auxiliary variables for c/a correction
    double cdist[3];
    cdist[0] = params->cacorrection.red;
    cdist[1] = 0.0;
    cdist[2] = params->cacorrection.blue;
    unsigned short** chorig[3];
    chorig[0] = original->r;
    chorig[1] = original->g;
    chorig[2] = original->b;
    unsigned short** chtrans[3];
    chtrans[0] = transformed->r;
    chtrans[1] = transformed->g;
    chtrans[2] = transformed->b;

	// auxiliary variables for distortion correction
	double a = params->distortion.amount;

	// auxiliary variables for rotation
	double cost = cos(params->rotate.degree * 3.14/180.0);
	double sint = sin(params->rotate.degree * 3.14/180.0);

	// auxiliary variables for vignetting
	double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2.;
	double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
	double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;
	double mul = (1.0-v) / tanh(b);
	bool dovign = params->vignetting.amount != 0;

	// auxiliary variables for vertical perspective correction
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * 3.14;
    double vpteta  = fabs(vpalpha-3.14/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

	// auxiliary variables for horizontal perspective correction
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * 3.14;
    double hpteta  = fabs(hpalpha-3.14/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

	double ascale = params->commonTrans.autofill ? getTransformAutoFill (oW, oH) : 1.0;

	// main cycle
	#pragma omp parallel for if (multiThread)
    for (int y=0; y<transformed->height; y++) {
        for (int x=0; x<transformed->width; x++) {
            double x_d = ascale * (x + cx - w2);		// centering x coord & scale
            double y_d = ascale * (y + cy - h2);		// centering y coord & scale

            // horizontal perspective transformation
            y_d = y_d * maxRadius / (maxRadius + x_d*hptanpt);
            x_d = x_d * maxRadius * hpcospt / (maxRadius + x_d*hptanpt);

            // vertical perspective transformation
            x_d = x_d * maxRadius / (maxRadius - y_d*vptanpt);
            y_d = y_d * maxRadius * vpcospt / (maxRadius - y_d*vptanpt);

            // rotate
            double Dxc = x_d * cost - y_d * sint;
            double Dyc = x_d * sint + y_d * cost;

            // distortion correction
            double r = sqrt(Dxc*Dxc + Dyc*Dyc);
            double s = 1.0 - a + a * r ;

            for (int c=0; c<3; c++) {

                double Dx = Dxc * (s + cdist[c]);
                double Dy = Dyc * (s + cdist[c]);

                // de-center
                Dx += w2;
                Dy += h2;

				// Extract integer and fractions of source screen coordinates
				int xc = (int)Dx; Dx -= (double)xc; xc -= sx;
				int yc = (int)Dy; Dy -= (double)yc; yc -= sy;

				// Convert only valid pixels
				if (yc>=0 && yc<original->height && xc>=0 && xc<original->width) {

					// multiplier for vignetting correction
					double vignmul = 1.0;
					if (dovign)
						vignmul /= (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

					if (yc > 0 && yc < original->height-2 && xc > 0 && xc < original->width-2)   // all interpolation pixels inside image
                        cubintch (chorig[c], xc-1, yc-1, Dx, Dy, &(chtrans[c][y][x]), vignmul);
					else { // edge pixels
						int y1 = CLIPTO(yc,   0, original->height-1);
						int y2 = CLIPTO(yc+1, 0, original->height-1);
						int x1 = CLIPTO(xc,   0, original->width-1);
						int x2 = CLIPTO(xc+1, 0, original->width-1);
                        int val = vignmul*(chorig[c][y1][x1]*(1.0-Dx)*(1.0-Dy) + chorig[c][y1][x2]*Dx*(1.0-Dy) + chorig[c][y2][x1]*(1.0-Dx)*Dy + chorig[c][y2][x2]*Dx*Dy);
                        chtrans[c][y][x] = CLIP(val);
					}
				}
				else
					// not valid (source pixel x,y not inside source image, etc.)
					chtrans[c][y][x] = 0;
			}
        }
    }
}

void ImProcFunctions::simpltransform (Image16* original, Image16* transformed, int cx, int cy, int sx, int sy, int oW, int oH) {

	double  w2 = (double) oW  / 2.0 - 0.5;
	double  h2 = (double) oH  / 2.0 - 0.5;

	// auxiliary variables for distortion correction
	double a = params->distortion.amount;

	// auxiliary variables for rotation
	double cost = cos(params->rotate.degree * 3.14/180.0);
	double sint = sin(params->rotate.degree * 3.14/180.0);

	// auxiliary variables for vignetting
	double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2.;
	double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
	double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;
	double mul = (1.0-v) / tanh(b);
	bool dovign = params->vignetting.amount != 0;

	// auxiliary variables for vertical perspective correction
    double vpdeg = params->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90 - vpdeg) / 180.0 * 3.14;
    double vpteta  = fabs(vpalpha-3.14/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

	// auxiliary variables for horizontal perspective correction
    double hpdeg = params->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90 - hpdeg) / 180.0 * 3.14;
    double hpteta  = fabs(hpalpha-3.14/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

	double ascale = params->commonTrans.autofill ? getTransformAutoFill (oW, oH) : 1.0;

    // main cycle
	#pragma omp parallel for if (multiThread)
    for (int y=0; y<transformed->height; y++) {
        for (int x=0; x<transformed->width; x++) {
            double y_d = ascale * (y + cy - h2);		// centering y coord & scale
            double x_d = ascale * (x + cx - w2);		// centering x coord & scale

            // horizontal perspective transformation
            y_d = y_d * maxRadius / (maxRadius + x_d*hptanpt);
            x_d = x_d * maxRadius * hpcospt / (maxRadius + x_d*hptanpt);

            // vertical perspective transformation
            x_d = x_d * maxRadius / (maxRadius - y_d*vptanpt);
            y_d = y_d * maxRadius * vpcospt / (maxRadius - y_d*vptanpt);

			// rotate
            double Dx = x_d * cost - y_d * sint;
            double Dy = x_d * sint + y_d * cost;

            // distortion correction
            double r = sqrt(Dx*Dx + Dy*Dy);
            double s = 1.0 - a + a * r ;
	        Dx *= s;
            Dy *= s;

            // de-center
            Dx += w2;
            Dy += h2;

            // Extract integer and fractions of source screen coordinates
            int xc = (int)Dx; Dx -= (double)xc; xc -= sx;
            int yc = (int)Dy; Dy -= (double)yc; yc -= sy;

            // Convert only valid pixels
            if (yc>=0 && yc<original->height && xc>=0 && xc<original->width) {

                // multiplier for vignetting correction
            	double vignmul = 1.0;
                if (dovign)
                	vignmul /= (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                if (yc < original->height-1 && xc < original->width-1) {  // all interpolation pixels inside image
                    int r = vignmul*(original->r[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->r[yc][xc+1]*Dx*(1.0-Dy) + original->r[yc+1][xc]*(1.0-Dx)*Dy + original->r[yc+1][xc+1]*Dx*Dy);
                    int g = vignmul*(original->g[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->g[yc][xc+1]*Dx*(1.0-Dy) + original->g[yc+1][xc]*(1.0-Dx)*Dy + original->g[yc+1][xc+1]*Dx*Dy);
                    int b = vignmul*(original->b[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->b[yc][xc+1]*Dx*(1.0-Dy) + original->b[yc+1][xc]*(1.0-Dx)*Dy + original->b[yc+1][xc+1]*Dx*Dy);
                    transformed->r[y][x] = CLIP(r);
                    transformed->g[y][x] = CLIP(g);
                    transformed->b[y][x] = CLIP(b);
                }
                else { // edge pixels
                	int y1 = CLIPTO(yc,   0, original->height-1);
                	int y2 = CLIPTO(yc+1, 0, original->height-1);
                	int x1 = CLIPTO(xc,   0, original->width-1);
                	int x2 = CLIPTO(xc+1, 0, original->width-1);
                    int r = vignmul*(original->r[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->r[y1][x2]*Dx*(1.0-Dy) + original->r[y2][x1]*(1.0-Dx)*Dy + original->r[y2][x2]*Dx*Dy);
                    int g = vignmul*(original->g[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->g[y1][x2]*Dx*(1.0-Dy) + original->g[y2][x1]*(1.0-Dx)*Dy + original->g[y2][x2]*Dx*Dy);
                    int b = vignmul*(original->b[y1][x1]*(1.0-Dx)*(1.0-Dy) + original->b[y1][x2]*Dx*(1.0-Dy) + original->b[y2][x1]*(1.0-Dx)*Dy + original->b[y2][x2]*Dx*Dy);
                    transformed->r[y][x] = CLIP(r);
                    transformed->g[y][x] = CLIP(g);
                    transformed->b[y][x] = CLIP(b);
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

double ImProcFunctions::getTransformAutoFill (int oW, int oH) {

	double scaleU = 1.0;
	double scaleL = 0.001;
	while (scaleU - scaleL > 0.001) {
		double scale = (scaleU + scaleL) / 2.0;

        int orx, ory, orw, orh;
        bool clipped = transCoord (oW, oH, 0, 0, oW, oH, orx, ory, orw, orh, scale);

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
bool ImProcFunctions::needsTransform () {

	return needsCA () || needsDistortion () || needsRotation () || needsPerspective () || needsVignetting ();
}


}

