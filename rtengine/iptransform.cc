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
#include <omp.h>

namespace rtengine {

#undef CMAXVAL
#undef MAX
#undef MIN
#undef CLIP
#undef CLIPTOC

#define CMAXVAL 0xffff
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define CLIPTOC(a,b,c,d) ((a)>=(b)?((a)<=(c)?(a):((c),d=true)):((b),d=true))

extern const Settings* settings;

void ImProcFunctions::vignetting_ (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to) {

  int oW = sizes.oW;
  int oH = sizes.oH;
  int cx = sizes.cx;
  int cy = sizes.cy;

  double  w2 = (double) oW  / 2.0 - 0.5;
  double  h2 = (double) oH  / 2.0 - 0.5;

  double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;

  double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
  double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;

  double mul = (1.0-v) / tanh(b);

  int val;
  for (int y=row_from; y<row_to; y++) {
      double y_d = (double) (y + cy) - h2 ;
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

void ImProcFunctions::vignetting (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int oW, int oH) {

    STemp sizes;
    sizes.cx = cx;
    sizes.cy = cy;
    sizes.oW = oW;
    sizes.oH = oH;

    if (settings->dualThreadEnabled) {
        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::vignetting_), original, transformed, params, sizes, 0, transformed->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::vignetting_), original, transformed, params, sizes, transformed->height/2, transformed->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else
        vignetting_ (original, transformed, params, sizes, 0, transformed->height);
}

#include "cubint.cc"
void ImProcFunctions::transform_ (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to) {

  int oW = sizes.oW;
  int oH = sizes.oH;
  int cx = sizes.cx;
  int cy = sizes.cy;
  int sx = sizes.sx;
  int sy = sizes.sy;

  double  w2 = (double) oW  / 2.0 - 0.5;
  double  h2 = (double) oH  / 2.0 - 0.5;

  double cost = cos(params->rotate.degree * 3.14/180.0);
  double sint = sin(params->rotate.degree * 3.14/180.0);

  double  max_x = (double) (sx + original->width - 1);
  double  max_y = (double) (sy + original->height - 1);
  double  min_x = (double) sx;
  double  min_y = (double) sy;

  const int n2 = 2;
  const int n = 4;

  int mix  = original->width - 1; // maximum x-index src
  int miy  = original->height - 1;// maximum y-index src
  int mix2 = mix +1 - n;
  int miy2 = miy +1 - n;

  double scale = (oW>oH) ? (double)oW / 2.0 : (double)oH / 2.0 ;
  double radius = sqrt( (double)( oW*oW + oH*oH ) );
  radius /= (oW<oH) ? oW : oH;

  double a = params->distortion.amount;

  double d = 1.0 - a;

    // magnify image to keep size
    double rotmagn = 1.0;
    if (params->rotate.fill) {
        double beta = atan((double)MIN(oH,oW)/MAX(oW,oH));
        rotmagn = sin(beta) / sin(fabs(params->rotate.degree) * 3.14/180.0 + beta);
    }
    // 1. check upper and lower border
    double d1 = rotmagn - a*h2/scale;
    double d2 = rotmagn - a*w2/scale;
    double d3 = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;
    d = MIN(d,MIN(d1,MIN(d2,d3)));

    // auxilary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2 / scale;

    double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;

    double mul = (1.0-v) / tanh(b);

    // main cycle
    double eps = 1e-10;
    bool calc_r=( (fabs(a)>eps) || (fabs(1.0-v)>eps) );
    bool do_vign = (fabs(1.0-v)>eps);

    for (int y=row_from; y<row_to; y++) {
        double y_d = (double) (y + cy) - h2 ;
        for (int x=0; x<transformed->width; x++) {
            double x_d = (double) (x + cx) - w2 ;

            double r=0.0;
            double s = d;//10000.0;
	    if (calc_r)
	    {
	            r=(sqrt(x_d*x_d + y_d*y_d)) / scale;
	            if (r<radius)
	            s += a * r ;
	    }

            double Dx = s*(x_d * cost - y_d * sint) + w2;
            double Dy = s*(x_d * sint + y_d * cost) + h2;

            if (fabs(Dx)<eps) Dx = 0;
            if (fabs(Dy)<eps) Dy = 0;
            if (fabs(Dx-max_x)<eps) Dx = nextafter(max_x,0);
            if (fabs(Dy-max_y)<eps) Dy = nextafter(max_y,0);

            bool valid = !((Dx >= max_x)   || (Dy >= max_y) || (Dx < min_x) || (Dy < min_y));

            // Convert only valid pixels
            if (valid) {
                // Extract integer and fractions of source screen coordinates
                int xc  =  (int) (Dx); Dx -= (double)xc;
                int yc  =  (int) (Dy); Dy -= (double)yc;
                int ys = yc +1 - n2 - sy; // smallest y-index used for interpolation
                int xs = xc +1 - n2 - sx; // smallest x-index used for interpolation

                double vignmul = 1.0;
		if (do_vign) vignmul /= (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                if (ys >= 0 && ys <= miy2 && xs >= 0 && xs <= mix2)   // all interpolation pixels inside image
                    cubint (original, xs, ys, Dx, Dy, &(transformed->r[y][x]), &(transformed->g[y][x]), &(transformed->b[y][x]), vignmul);
                else { // edge pixels
                    int y1 = (yc>0) ? yc : 0;
                    if (y1>miy) y1 = miy;
                    int y2 = (yc<miy) ? yc+1 : miy;
                    if (y2<0) y2 = 0;
                    int x1 = (xc>0) ? xc : 0;
                    if (x1>mix) x1 = mix;
                    int x2 = (xc<mix) ? xc+1 : mix;
                    if (x2<0) x2 = 0;
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

void ImProcFunctions::simpltransform_ (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to) {

  int oW = sizes.oW;
  int oH = sizes.oH;
  int cx = sizes.cx;
  int cy = sizes.cy;
  int sx = sizes.sx;
  int sy = sizes.sy;

  double  w2 = (double) oW  / 2.0 - 0.5;
  double  h2 = (double) oH  / 2.0 - 0.5;

  double cost = cos(params->rotate.degree * 3.14/180.0);
  double sint = sin(params->rotate.degree * 3.14/180.0);

  double  max_x = (double) (sx + original->width - 1);
  double  max_y = (double) (sy + original->height - 1);
  double  min_x = (double) sx;
  double  min_y = (double) sy;

  const int n2 = 2;
  const int n = 2;

  int mix  = original->width - 1; // maximum x-index src
  int miy  = original->height - 1;// maximum y-index src
  int mix2 = mix +1 - n;
  int miy2 = miy +1 - n;

  double scale = (oW>oH) ? (double)oW / 2.0 : (double)oH / 2.0 ;
  double radius = sqrt( (double)( oW*oW + oH*oH ) );
  radius /= (oW<oH) ? oW : oH;

  double a = params->distortion.amount;

  double d = 1.0 - a;


    // magnify image to keep size
    double rotmagn = 1.0;
    if (params->rotate.fill) {
        double beta = atan((double)MIN(oH,oW)/MAX(oW,oH));
        rotmagn = sin(beta) / sin(fabs(params->rotate.degree) * 3.14/180.0 + beta);
    }
    // 1. check upper and lower border
    double d1r = rotmagn - a*h2/scale - params->cacorrection.red;
    double d2r = rotmagn - a*w2/scale - params->cacorrection.red;
    double d3r = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.red;
    double dr = MIN(d,MIN(d1r,MIN(d2r,d3r)));
    double d1b = rotmagn - a*h2/scale - params->cacorrection.blue;
    double d2b = rotmagn - a*w2/scale - params->cacorrection.blue;
    double d3b = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.blue;
    double db = MIN(d,MIN(d1b,MIN(d2b,d3b)));
    double d1g = rotmagn - a*h2/scale;
    double d2g = rotmagn - a*w2/scale;
    double d3g = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;
    double dg = MIN(d,MIN(d1g,MIN(d2g,d3g)));

    d = MIN(dg,MIN(dr,db));

    // auxilary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2 / scale;

    double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;

    double mul = (1.0-v) / tanh(b);

    // main cycle
    double eps = 1e-10;
    bool calc_r=( (fabs(a)>eps) || (fabs(1.0-v)>eps) );
    bool do_vign = (fabs(1.0-v)>eps);

    for (int y=row_from; y<row_to; y++) {
        double y_d = (double) (y + cy) - h2 ;
        for (int x=0; x<transformed->width; x++) {
            double x_d = (double) (x + cx) - w2 ;

            double r=0.0;
            double s = d;//10000.0;
	    if (calc_r)
	    {
	            r=(sqrt(x_d*x_d + y_d*y_d)) / scale;
	            if (r<radius)
	            s += a * r ;
	    }

            double Dx = s*(x_d * cost - y_d * sint) + w2;
            double Dy = s*(x_d * sint + y_d * cost) + h2;

            if (fabs(Dx)<eps) Dx = 0;
            if (fabs(Dy)<eps) Dy = 0;
            if (fabs(Dx-max_x)<eps) Dx = nextafter(max_x,0);
            if (fabs(Dy-max_y)<eps) Dy = nextafter(max_y,0);

            bool valid = !((Dx >= max_x)   || (Dy >= max_y) || (Dx < min_x) || (Dy < min_y));

            // Convert only valid pixels
            if (valid) {
                // Extract integer and fractions of source screen coordinates
                int xc  =  (int) (Dx); Dx -= (double)xc;
                int yc  =  (int) (Dy); Dy -= (double)yc;
                int ys = yc +1 - n2 - sy; // smallest y-index used for interpolation
                int xs = xc +1 - n2 - sx; // smallest x-index used for interpolation

                double vignmul = 1.0;
		if (do_vign) vignmul /= (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                if (ys >= 0 && ys <= miy2 && xs >= 0 && xs <= mix2 && yc < miy-1) {   // all interpolation pixels inside image

                    int r = vignmul*(original->r[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->r[yc][xc+1]*Dx*(1.0-Dy) + original->r[yc+1][xc]*(1.0-Dx)*Dy + original->r[yc+1][xc+1]*Dx*Dy);
                    int g = vignmul*(original->g[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->g[yc][xc+1]*Dx*(1.0-Dy) + original->g[yc+1][xc]*(1.0-Dx)*Dy + original->g[yc+1][xc+1]*Dx*Dy);
                    int b = vignmul*(original->b[yc][xc]*(1.0-Dx)*(1.0-Dy) + original->b[yc][xc+1]*Dx*(1.0-Dy) + original->b[yc+1][xc]*(1.0-Dx)*Dy + original->b[yc+1][xc+1]*Dx*Dy);
                    transformed->r[y][x] = CLIP(r);
                    transformed->g[y][x] = CLIP(g);
                    transformed->b[y][x] = CLIP(b);
                }
                else { // edge pixels
                    int y1 = (yc>0) ? yc : 0;
                    if (y1>miy) y1 = miy;
                    int y2 = (yc<miy) ? yc+1 : miy;
                    if (y2<0) y2 = 0;
                    int x1 = (xc>0) ? xc : 0;
                    if (x1>mix) x1 = mix;
                    int x2 = (xc<mix) ? xc+1 : mix;
                    if (x2<0) x2 = 0;
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
void ImProcFunctions::transform_sep_ (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to) {

  int oW = sizes.oW;
  int oH = sizes.oH;
  int cx = sizes.cx;
  int cy = sizes.cy;
  int sx = sizes.sx;
  int sy = sizes.sy;

  double  w2 = (double) oW  / 2.0 - 0.5;
  double  h2 = (double) oH  / 2.0 - 0.5;

  double cost = cos(params->rotate.degree * 3.14/180.0);
  double sint = sin(params->rotate.degree * 3.14/180.0);

  double  max_x = (double) (sx + original->width - 1);
  double  max_y = (double) (sy + original->height - 1);
  double  min_x = (double) sx;
  double  min_y = (double) sy;

  const int n2 = 2;
  const int n = 4;

  int mix  = original->width - 1; // maximum x-index src
  int miy  = original->height - 1;// maximum y-index src
  int mix2 = mix +1 - n;
  int miy2 = miy +1 - n;

  double scale = (oW>oH) ? (double)oW / 2.0 : (double)oH / 2.0 ;
  double radius = sqrt( (double)( oW*oW + oH*oH ) );
  radius /= (oW<oH) ? oW : oH;

  double a = params->distortion.amount;
    double d = 1.0 - a;

    double cdist[3];
    cdist[0] = params->cacorrection.red;
    cdist[1] = 0.0;
    cdist[2] = params->cacorrection.blue;

    // magnify image to keep size
    double rotmagn = 1.0;
    if (params->rotate.fill) {
        double beta = atan((double)MIN(oH,oW)/MAX(oW,oH));
        rotmagn = sin(beta) / sin(fabs(params->rotate.degree) * 3.14/180.0 + beta);
    }
    // 1. check upper and lower border
    double d1r = rotmagn - a*h2/scale - params->cacorrection.red;
    double d2r = rotmagn - a*w2/scale - params->cacorrection.red;
    double d3r = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.red;
    double dr = MIN(d,MIN(d1r,MIN(d2r,d3r)));
    double d1b = rotmagn - a*h2/scale - params->cacorrection.blue;
    double d2b = rotmagn - a*w2/scale - params->cacorrection.blue;
    double d3b = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.blue;
    double db = MIN(d,MIN(d1b,MIN(d2b,d3b)));
    double d1g = rotmagn - a*h2/scale;
    double d2g = rotmagn - a*w2/scale;
    double d3g = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;
    double dg = MIN(d,MIN(d1g,MIN(d2g,d3g)));

    d = MIN(dg,MIN(dr,db));

    unsigned short** chorig[3];
    chorig[0] = original->r;
    chorig[1] = original->g;
    chorig[2] = original->b;

    unsigned short** chtrans[3];
    chtrans[0] = transformed->r;
    chtrans[1] = transformed->g;
    chtrans[2] = transformed->b;


    // auxilary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2 / scale;

    double v = 1.0 - params->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + params->vignetting.radius * 7.0 / 100.0;

    double mul = (1.0-v) / tanh(b);

    // main cycle
    double eps = 1e-10;
    for (int y=row_from; y<row_to; y++) {
        double y_d = (double) (y + cy) - h2 ;
        for (int x=0; x<transformed->width; x++) {
            double x_d = (double) (x + cx) - w2 ;

            double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
            double s = 10000.0;
            if (r<radius)
	            s = a * r + d;

            double vignmul = 1.0 / (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

            for (int c=0; c<3; c++) {

                double Dx = (s + cdist[c]) * (x_d * cost - y_d * sint) + w2;
                double Dy = (s + cdist[c]) * (x_d * sint + y_d * cost) + h2;

                if (fabs(Dx)<eps) Dx = 0;
                if (fabs(Dy)<eps) Dy = 0;
                if (fabs(Dx-max_x)<eps) Dx = nextafter(max_x,0);
                if (fabs(Dy-max_y)<eps) Dy = nextafter(max_y,0);

                bool valid = !((Dx >= max_x)   || (Dy >= max_y) || (Dx < min_x) || (Dy < min_y));

                // Convert only valid pixels
                if (valid) {
                    // Extract integer and fractions of source screen coordinates
                    int xc  =  (int) (Dx); Dx -= (double)xc;
                    int yc  =  (int) (Dy); Dy -= (double)yc;
                    int ys = yc +1 - n2 - sy; // smallest y-index used for interpolation
                    int xs = xc +1 - n2 - sx; // smallest x-index used for interpolation

                    if (ys >= 0 && ys <= miy2 && xs >= 0 && xs <= mix2)  // all interpolation pixels inside image
                        cubintch (chorig[c], xs, ys, Dx, Dy, &(chtrans[c][y][x]), vignmul);
                    else {// edge pixels, linear interpolation
                        int y1 = (yc>0) ? yc : 0;
                        if (y1>miy) y1 = miy;
                        int y2 = (yc<miy) ? yc+1 : miy;
                        if (y2<0) y2 = 0;
                        int x1 = (xc>0) ? xc : 0;
                        if (x1>mix) x1 = mix;
                        int x2 = (xc<mix) ? xc+1 : mix;
                        if (x2<0) x2 = 0;
                        int val = vignmul*(chorig[c][y1][x1]*(1.0-Dx)*(1.0-Dy) + chorig[c][y1][x2]*Dx*(1.0-Dy) + chorig[c][y2][x1]*(1.0-Dx)*Dy + chorig[c][y2][x2]*Dx*Dy);
                        chtrans[c][y][x] = CLIP(val);
                    }
                }
                else // not valid (source pixel x,y not inside source image, etc.)
                    chtrans[c][y][x] = 0;
            }
        }
    }
}

bool ImProcFunctions::transCoord (const ProcParams* params, int W, int H, std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue) {

    bool clipresize = true;
    bool clipped = false;

    red.clear ();
    green.clear ();
    blue.clear ();
    bool needstransform  = 0;// fabs(params->rotate.degree)>1e-15 || fabs(params->distortion.amount)>1e-15 || fabs(params->cacorrection.red)>1e-15 || fabs(params->cacorrection.blue)>1e-15;
    if (!needstransform) {
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
    double rW = W*params->resize.scale;
    double rH = H*params->resize.scale;
    double  w2 = (double) rW  / 2.0 - 0.5;
    double  h2 = (double) rH  / 2.0 - 0.5;
    double cost = cos(params->rotate.degree * 3.14/180.0);
    double sint = sin(params->rotate.degree * 3.14/180.0);

    double scale = (rW>rH) ? rW / 2.0 : rH / 2.0 ;
    double radius = sqrt ((double)(rW*rW + rH*rH ));
    radius /= (rW<rH) ? rW : rH;
    double a = params->distortion.amount;
    double d = 1.0 - a;

    // magnify image to keep size
    double rotmagn = 1.0;
    if (params->rotate.fill) {
        double beta = atan(MIN(rH,rW)/MAX(rW,rH));
        rotmagn = sin(beta) / sin(fabs(params->rotate.degree) * 3.14/180.0 + beta);
    }
    if (params->cacorrection.red==0 && params->cacorrection.blue==0) {
        // 1. check upper and lower border
        double d1 = rotmagn - a*h2/scale;
        double d2 = rotmagn - a*w2/scale;
        double d3 = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;
        d = MIN(d,MIN(d1,MIN(d2,d3)));

        for (int i=0; i<src.size(); i++) {
            double y_d = src[i].y - h2 ;
            double x_d = src[i].x - w2 ;
            double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
            double s = 10000.0;
            if (r<radius)
                s = a * r + d;
            red.push_back (Coord2D(s*(x_d * cost - y_d * sint) + w2, s*(x_d * sint + y_d * cost) + h2));
            green.push_back (Coord2D(s*(x_d * cost - y_d * sint) + w2, s*(x_d * sint + y_d * cost) + h2));
            blue.push_back (Coord2D(s*(x_d * cost - y_d * sint) + w2, s*(x_d * sint + y_d * cost) + h2));
        }
    }
    else {
        double cdist[3];
        cdist[0] = params->cacorrection.red;
        cdist[1] = 0.0;
        cdist[2] = params->cacorrection.blue;

        // 1. check upper and lower border
        double d1r = rotmagn - a*h2/scale - params->cacorrection.red;
        double d2r = rotmagn - a*w2/scale - params->cacorrection.red;
        double d3r = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.red;
        double dr = MIN(d,MIN(d1r,MIN(d2r,d3r)));
        double d1b = rotmagn - a*h2/scale - params->cacorrection.blue;
        double d2b = rotmagn - a*w2/scale - params->cacorrection.blue;
        double d3b = rotmagn - a*sqrt(h2*h2+w2*w2) / scale - params->cacorrection.blue;
        double db = MIN(d,MIN(d1b,MIN(d2b,d3b)));
        double d1g = rotmagn - a*h2/scale;
        double d2g = rotmagn - a*w2/scale;
        double d3g = rotmagn - a*sqrt(h2*h2+w2*w2) / scale;
        double dg = MIN(d,MIN(d1g,MIN(d2g,d3g)));

        d = MIN(dg,MIN(dr,db));

        for (int i=0; i<src.size(); i++) {
            double y_d = src[i].y - h2 ;
            double x_d = src[i].x - w2 ;
            double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
            double s = 10000.0;
            if (r<radius)
                s = a * r + d;
            src[i].x = s*(x_d * cost - y_d * sint) + w2;
            src[i].y  = s*(x_d * sint + y_d * cost) + h2;

            red.push_back (Coord2D((s+cdist[0])*(x_d * cost - y_d * sint) + w2, (s+cdist[0])*(x_d * sint + y_d * cost) + h2));
            green.push_back (Coord2D((s+cdist[1])*(x_d * cost - y_d * sint) + w2, (s+cdist[1])*(x_d * sint + y_d * cost) + h2));
            blue.push_back (Coord2D((s+cdist[2])*(x_d * cost - y_d * sint) + w2, (s+cdist[2])*(x_d * sint + y_d * cost) + h2));
        }
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

bool ImProcFunctions::transCoord (const ProcParams* params, int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv) {

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

    std::vector<Coord2D> r, g, b;

    bool result = transCoord (params, W, H, corners, r, g, b);

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

void ImProcFunctions::transform (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int sx, int sy, int oW, int oH) {

    STemp sizes;
    sizes.cx = 0;//cx;
    sizes.cy = 0;//cy;
    sizes.oW = oW;
    sizes.oH = oH;
    sizes.sx = 0;//sx;
    sizes.sy = 0;//sy;

    if (params->cacorrection.red==0 && params->cacorrection.blue==0) {
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::transform_), original, transformed, params, sizes, 0, transformed->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::transform_), original, transformed, params, sizes, transformed->height/2, transformed->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else
            transform_ (original, transformed, params, sizes, 0, transformed->height);
    }
    else {
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::transform_sep_), original, transformed, params, sizes, 0, transformed->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::transform_sep_), original, transformed, params, sizes, transformed->height/2, transformed->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else
            transform_sep_ (original, transformed, params, sizes, 0, transformed->height);
    }
}

void ImProcFunctions::simpltransform (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int sx, int sy, int oW, int oH) {

    STemp sizes;
    sizes.cx = 0;//cx;
    sizes.cy = 0;//cy;
    sizes.oW = oW;
    sizes.oH = oH;
    sizes.sx = 0;//sx;
    sizes.sy = 0;//sy;

    if (settings->dualThreadEnabled) {
        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::simpltransform_), original, transformed, params, sizes, 0, transformed->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &ImProcFunctions::simpltransform_), original, transformed, params, sizes, transformed->height/2, transformed->height), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else
        simpltransform_ (original, transformed, params, sizes, 0, transformed->height);
}
/*void ImProcFunctions::transform (Image16* original, Image16* transformed, const ProcParams* params, int ox, int oy) {

  if (!transformed)
    return;

  int oW = W, oH = H, tW = W, tH = H;

  double  w2 = (double) tW / 2.0 - 0.5;
  double  h2 = (double) tH / 2.0 - 0.5;
  double  sw2 = (double) oW  / 2.0 - 0.5;
  double  sh2 = (double) oH  / 2.0 - 0.5;

  double cost = cos(params->rotate_fine * 3.14/180.0);
  double sint = sin(params->rotate_fine * 3.14/180.0);

  double  max_x = (double) oW;
  double  max_y = (double) oH;
  double  min_x =  0.0;
  double  min_y =  0.0;

  const int n2 = 2;
  const int n = 4;

  int mix  = oW - 1; // maximum x-index src
  int miy  = oH - 1;// maximum y-index src
  int mix2 = mix +1 - n;
  int miy2 = miy +1 - n;

  double scale = (tW>tH) ? (double)tW / 2.0 : (double)tH / 2.0 ;
  double radius = sqrt( (double)( tW*tW + tH*tH ) );
  radius /= (tW<tH) ? tW : tH;

  double a = params->lens_distortion;

  for (int y=0; y<transformed->height; y++) {
    double y_d = (double) y + oy - h2 ;

    for (int x=0; x<transformed->width; x++) {
      double x_d = (double) x + ox - w2 ;

      double r = (sqrt(x_d*x_d + y_d*y_d)) / scale;
      double s = 10000.0;
      if (r<radius)
	        s = a * r + 1.0 - a;

      double Dx = s*(x_d * cost - y_d * sint) + sw2;
      double Dy = s*(x_d * sint + y_d * cost) + sh2;

      bool valid = !((Dx >= max_x)   || (Dy >= max_y) || (Dx < min_x) || (Dy < min_y));

      // Convert only valid pixels
      if (valid) {
        // Extract integer and fractions of source screen coordinates
        int xc  =  (int) floor (Dx) ; Dx -= (double)xc;
        int yc  =  (int) floor (Dy) ; Dy -= (double)yc;
        int ys = yc +1 - n2 ; // smallest y-index used for interpolation
        int xs = xc +1 - n2 ; // smallest x-index used for interpolation

        unsigned short sr[2][2], sg[2][2], sb[2][2];

        if (ys >= 0 && ys <= miy2 && xs >= 0 && xs <= mix2)   // all interpolation pixels inside image
          cubint (original, xs, ys, Dx, Dy, &(transformed->r[y][x]), &(transformed->g[y][x]), &(transformed->b[y][x]));
        else { // edge pixels
          transformed->r[y][x] = 0;
          transformed->g[y][x] = 0;
          transformed->b[y][x] = 0;
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
}*/

}

