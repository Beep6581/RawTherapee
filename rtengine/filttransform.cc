/*
 * filttransform.cc
 *
 *  Created on: Oct 4, 2010
 *      Author: gabor
 */

#include "filttransform.h"
#include "rtengine.h"
#include "macros.h"

namespace rtengine {

TransformFilterDescriptor transformFilterDescriptor;

TransformFilterDescriptor::TransformFilterDescriptor ()
	: FilterDescriptor ("Transform", MultiImage::RGB, MultiImage::RGB, true) {

	addTriggerEvent (EvROTDegree);
    addTriggerEvent (EvTransAutoFill);
    addTriggerEvent (EvDISTAmount);
    addTriggerEvent (EvVignetting);
    addTriggerEvent (EvPerspCorr);
    addTriggerEvent (EvCACorr);
}

void TransformFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new TransformFilter ());
}

TransformFilter::TransformFilter ()
	: Filter (&transformFilterDescriptor) {
}

ImageView TransformFilter::calculateSourceImageView (const ImageView& requestedImView) {

    ImageView src;
    transCoord (getFullImageSize (), requestedImView, src);
    return src;
}

void TransformFilter::reverseTransPoint (int x, int y, int& xv, int& yv) {

    std::vector<Coord2D> corners (1);
    corners[0].set (x, y);
    std::vector<Coord2D> r, g, b;
    transCoord (getFullImageSize(), corners, r, g, b);
    xv = g[0].x;
    yv = g[0].y;
}

void TransformFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    if (!needsTransform()) {
        if (targetImage != sourceImage)
            targetImage->copyFrom (sourceImage);
        return;
    }

    if (!(needsCA() || needsDistortion() || needsRotation() || needsPerspective()) && needsVignetting())
        vignetting (sourceImage, targetImage);
    else if (!needsCA()) {
        if (getScale() < 1.0)
            simpltransform (sourceImage, targetImage);
        else
            transformNonSep (sourceImage, targetImage);
    }
    else
        transformSep (sourceImage, targetImage);
}

void TransformFilter::transformSep (MultiImage* sourceImage, MultiImage* targetImage) {

    Dim fsize = getFullImageSize ();
    double scale = getScale ();

    int oW = (int)round (fsize.width * scale);
    int oH = (int)round (fsize.height * scale);

    int cx = getScaledTargetImageView().x;
    int cy = getScaledTargetImageView().y;

    int sx = getScaledSourceImageView().x;
    int sy = getScaledSourceImageView().y;

    double  w2 = oW / 2.0 - 0.5;
    double  h2 = oH  / 2.0 - 0.5;

    // auxiliary variables for c/a correction
    double cdist[3];
    cdist[0] = procParams->cacorrection.red;
    cdist[1] = 0.0;
    cdist[2] = procParams->cacorrection.blue;
    unsigned short** chorig[3];
    chorig[0] = sourceImage->r;
    chorig[1] = sourceImage->g;
    chorig[2] = sourceImage->b;
    unsigned short** chtrans[3];
    chtrans[0] = targetImage->r;
    chtrans[1] = targetImage->g;
    chtrans[2] = targetImage->b;

    // auxiliary variables for distortion correction
    double a = procParams->distortion.amount;

    // auxiliary variables for rotation
    double cost = cos(procParams->rotate.degree * 3.14/180.0);
    double sint = sin(procParams->rotate.degree * 3.14/180.0);

    // auxiliary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;
    double v = 1.0 - procParams->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + procParams->vignetting.radius * 7.0 / 100.0;
    double mul = (1.0-v) / tanh(b);
    bool dovign = procParams->vignetting.amount != 0;

    // auxiliary variables for vertical perspective correction
    double vpdeg = procParams->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * 3.14;
    double vpteta  = fabs(vpalpha-3.14/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

    // auxiliary variables for horizontal perspective correction
    double hpdeg = procParams->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * 3.14;
    double hpteta  = fabs(hpalpha-3.14/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

    double ascale = procParams->commonTrans.autofill ? getTransformAutoFill () : 1.0;

    // main cycle
    #pragma omp parallel for if (multiThread)
    for (int y=0; y<targetImage->height; y++) {
        for (int x=0; x<targetImage->width; x++) {
            double x_d = ascale * (x + cx - w2);        // centering x coord & scale
            double y_d = ascale * (y + cy - h2);        // centering y coord & scale

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
            double r = sqrt(Dxc*Dxc + Dyc*Dyc) / maxRadius;
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
                if (yc>=0 && yc<sourceImage->height && xc>=0 && xc<sourceImage->width) {

                    // multiplier for vignetting correction
                    double vignmul = 1.0;
                    if (dovign)
                        vignmul /= (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                    if (yc > 0 && yc < sourceImage->height-2 && xc > 0 && xc < sourceImage->width-2)   // all interpolation pixels inside image
                        cubintch (chorig[c], xc-1, yc-1, Dx, Dy, &(chtrans[c][y][x]), vignmul);
                    else { // edge pixels
                        int y1 = CLIPTO(yc,   0, sourceImage->height-1);
                        int y2 = CLIPTO(yc+1, 0, sourceImage->height-1);
                        int x1 = CLIPTO(xc,   0, sourceImage->width-1);
                        int x2 = CLIPTO(xc+1, 0, sourceImage->width-1);
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

void TransformFilter::transformNonSep (MultiImage* sourceImage, MultiImage* targetImage) {

    Dim fsize = getFullImageSize ();
    double scale = getScale ();

    int oW = (int)round (fsize.width * scale);
    int oH = (int)round (fsize.height * scale);

    int cx = getScaledTargetImageView().x;
    int cy = getScaledTargetImageView().y;

    int sx = getScaledSourceImageView().x;
    int sy = getScaledSourceImageView().y;

    double  w2 = (double) oW  / 2.0 - 0.5;
    double  h2 = (double) oH  / 2.0 - 0.5;

    // auxiliary variables for distortion correction
    double a = procParams->distortion.amount;

    // auxiliary variables for rotation
    double cost = cos(procParams->rotate.degree * 3.14/180.0);
    double sint = sin(procParams->rotate.degree * 3.14/180.0);

    // auxiliary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;
    double v = 1.0 - procParams->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + procParams->vignetting.radius * 7.0 / 100.0;
    double mul = (1.0-v) / tanh(b);
    bool dovign = procParams->vignetting.amount != 0;

    // auxiliary variables for vertical perspective correction
    double vpdeg = procParams->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * 3.14;
    double vpteta  = fabs(vpalpha-3.14/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

    // auxiliary variables for horizontal perspective correction
    double hpdeg = procParams->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * 3.14;
    double hpteta  = fabs(hpalpha-3.14/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

    double ascale = procParams->commonTrans.autofill ? getTransformAutoFill () : 1.0;

    // main cycle
    #pragma omp parallel for if (multiThread)
    for (int y=0; y<targetImage->height; y++) {
        for (int x=0; x<targetImage->width; x++) {
            double x_d = ascale * (x + cx - w2);        // centering x coord & scale
            double y_d = ascale * (y + cy - h2);        // centering y coord & scale

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
            double r = sqrt(Dx*Dx + Dy*Dy) / maxRadius;
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
            if (yc>=0 && yc<sourceImage->height && xc>=0 && xc<sourceImage->width) {

                // multiplier for vignetting correction
                double vignmul = 1.0;
                if (dovign)
                    vignmul /= (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                if (yc > 0 && yc < sourceImage->height-2 && xc > 0 && xc < sourceImage->width-2)   // all interpolation pixels inside image
                    cubint (sourceImage, xc-1, yc-1, Dx, Dy, &(targetImage->r[y][x]), &(targetImage->g[y][x]), &(targetImage->b[y][x]), vignmul);
                else { // edge pixels
                    int y1 = CLIPTO(yc,   0, sourceImage->height-1);
                    int y2 = CLIPTO(yc+1, 0, sourceImage->height-1);
                    int x1 = CLIPTO(xc,   0, sourceImage->width-1);
                    int x2 = CLIPTO(xc+1, 0, sourceImage->width-1);
                    int r = vignmul*(sourceImage->r[y1][x1]*(1.0-Dx)*(1.0-Dy) + sourceImage->r[y1][x2]*Dx*(1.0-Dy) + sourceImage->r[y2][x1]*(1.0-Dx)*Dy + sourceImage->r[y2][x2]*Dx*Dy);
                    int g = vignmul*(sourceImage->g[y1][x1]*(1.0-Dx)*(1.0-Dy) + sourceImage->g[y1][x2]*Dx*(1.0-Dy) + sourceImage->g[y2][x1]*(1.0-Dx)*Dy + sourceImage->g[y2][x2]*Dx*Dy);
                    int b = vignmul*(sourceImage->b[y1][x1]*(1.0-Dx)*(1.0-Dy) + sourceImage->b[y1][x2]*Dx*(1.0-Dy) + sourceImage->b[y2][x1]*(1.0-Dx)*Dy + sourceImage->b[y2][x2]*Dx*Dy);
                    targetImage->r[y][x] = CLIP(r);
                    targetImage->g[y][x] = CLIP(g);
                    targetImage->b[y][x] = CLIP(b);
                }
            }
            else {
                // not valid (source pixel x,y not inside source image, etc.)
                targetImage->r[y][x] = 0;
                targetImage->g[y][x] = 0;
                targetImage->b[y][x] = 0;
            }
        }
    }
}

void TransformFilter::simpltransform (MultiImage* sourceImage, MultiImage* targetImage) {

    Dim fsize = getFullImageSize ();
    double scale = getScale ();

    int oW = (int)round (fsize.width * scale);
    int oH = (int)round (fsize.height * scale);

    int cx = getScaledTargetImageView().x;
    int cy = getScaledTargetImageView().y;

    int sx = getScaledSourceImageView().x;
    int sy = getScaledSourceImageView().y;

    double  w2 = (double) oW  / 2.0 - 0.5;
    double  h2 = (double) oH  / 2.0 - 0.5;

    // auxiliary variables for distortion correction
    double a = procParams->distortion.amount;

    // auxiliary variables for rotation
    double cost = cos(procParams->rotate.degree * 3.14/180.0);
    double sint = sin(procParams->rotate.degree * 3.14/180.0);

    // auxiliary variables for vignetting
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;
    double v = 1.0 - procParams->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + procParams->vignetting.radius * 7.0 / 100.0;
    double mul = (1.0-v) / tanh(b);
    bool dovign = procParams->vignetting.amount != 0;

    // auxiliary variables for vertical perspective correction
    double vpdeg = procParams->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90 - vpdeg) / 180.0 * 3.14;
    double vpteta  = fabs(vpalpha-3.14/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);

    // auxiliary variables for horizontal perspective correction
    double hpdeg = procParams->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90 - hpdeg) / 180.0 * 3.14;
    double hpteta  = fabs(hpalpha-3.14/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

    double ascale = procParams->commonTrans.autofill ? getTransformAutoFill () : 1.0;

    // main cycle
    #pragma omp parallel for if (multiThread)
    for (int y=0; y<targetImage->height; y++) {
        for (int x=0; x<targetImage->width; x++) {
            double y_d = ascale * (y + cy - h2);        // centering y coord & scale
            double x_d = ascale * (x + cx - w2);        // centering x coord & scale

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
            double r = sqrt(Dx*Dx + Dy*Dy) / maxRadius;
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
            if (yc>=0 && yc<sourceImage->height && xc>=0 && xc<sourceImage->width) {

                // multiplier for vignetting correction
                double vignmul = 1.0;
                if (dovign)
                    vignmul /= (v + mul * tanh (b*(maxRadius-s*r) / maxRadius));

                if (yc < sourceImage->height-1 && xc < sourceImage->width-1) {  // all interpolation pixels inside image
                    int r = vignmul*(sourceImage->r[yc][xc]*(1.0-Dx)*(1.0-Dy) + sourceImage->r[yc][xc+1]*Dx*(1.0-Dy) + sourceImage->r[yc+1][xc]*(1.0-Dx)*Dy + sourceImage->r[yc+1][xc+1]*Dx*Dy);
                    int g = vignmul*(sourceImage->g[yc][xc]*(1.0-Dx)*(1.0-Dy) + sourceImage->g[yc][xc+1]*Dx*(1.0-Dy) + sourceImage->g[yc+1][xc]*(1.0-Dx)*Dy + sourceImage->g[yc+1][xc+1]*Dx*Dy);
                    int b = vignmul*(sourceImage->b[yc][xc]*(1.0-Dx)*(1.0-Dy) + sourceImage->b[yc][xc+1]*Dx*(1.0-Dy) + sourceImage->b[yc+1][xc]*(1.0-Dx)*Dy + sourceImage->b[yc+1][xc+1]*Dx*Dy);
                    targetImage->r[y][x] = CLIP(r);
                    targetImage->g[y][x] = CLIP(g);
                    targetImage->b[y][x] = CLIP(b);
                }
                else { // edge pixels
                    int y1 = CLIPTO(yc,   0, sourceImage->height-1);
                    int y2 = CLIPTO(yc+1, 0, sourceImage->height-1);
                    int x1 = CLIPTO(xc,   0, sourceImage->width-1);
                    int x2 = CLIPTO(xc+1, 0, sourceImage->width-1);
                    int r = vignmul*(sourceImage->r[y1][x1]*(1.0-Dx)*(1.0-Dy) + sourceImage->r[y1][x2]*Dx*(1.0-Dy) + sourceImage->r[y2][x1]*(1.0-Dx)*Dy + sourceImage->r[y2][x2]*Dx*Dy);
                    int g = vignmul*(sourceImage->g[y1][x1]*(1.0-Dx)*(1.0-Dy) + sourceImage->g[y1][x2]*Dx*(1.0-Dy) + sourceImage->g[y2][x1]*(1.0-Dx)*Dy + sourceImage->g[y2][x2]*Dx*Dy);
                    int b = vignmul*(sourceImage->b[y1][x1]*(1.0-Dx)*(1.0-Dy) + sourceImage->b[y1][x2]*Dx*(1.0-Dy) + sourceImage->b[y2][x1]*(1.0-Dx)*Dy + sourceImage->b[y2][x2]*Dx*Dy);
                    targetImage->r[y][x] = CLIP(r);
                    targetImage->g[y][x] = CLIP(g);
                    targetImage->b[y][x] = CLIP(b);
                }
            }
            else {
                // not valid (source pixel x,y not inside source image, etc.)
                targetImage->r[y][x] = 0;
                targetImage->g[y][x] = 0;
                targetImage->b[y][x] = 0;
            }
        }
    }
}

void TransformFilter::vignetting (MultiImage* sourceImage, MultiImage* targetImage) {

    Dim fsize = getFullImageSize ();
    double scale = getScale ();

    int oW = (int)round (fsize.width * scale);
    int oH = (int)round (fsize.height * scale);

    int cx = getScaledTargetImageView().x;
    int cy = getScaledTargetImageView().y;

    double  w2 = (double) oW  / 2.0 - 0.5;
    double  h2 = (double) oH  / 2.0 - 0.5;

    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;

    double v = 1.0 - procParams->vignetting.amount * 3.0 / 400.0;
    double b = 1.0 + procParams->vignetting.radius * 7.0 / 100.0;

    double mul = (1.0-v) / tanh(b);

    #pragma omp parallel for if (multiThread)
    for (int y=0; y<targetImage->height; y++) {
        double y_d = (double) (y + cy) - h2 ;
        int val;
        for (int x=0; x<targetImage->width; x++) {
            double x_d = (double) (x + cx) - w2 ;
            double r = sqrt(x_d*x_d + y_d*y_d);
            double vign = v + mul * tanh (b*(maxRadius-r) / maxRadius);
            val = sourceImage->r[y][x] / vign;
            targetImage->r[y][x] = CLIP(val);
            val =  sourceImage->g[y][x] / vign;
            targetImage->g[y][x] = CLIP(val);
            val = sourceImage->b[y][x] / vign;
            targetImage->b[y][x] = CLIP(val);
        }
    }
}

double TransformFilter::getTransformAutoFill () {

    Dim fsize = getFullImageSize ();
    ImageView dst (0, 0, fsize.width, fsize.height);
    ImageView src;
    double scaleU = 1.0;
    double scaleL = 0.001;
    while (scaleU - scaleL > 0.001) {
        double scale = (scaleU + scaleL) / 2.0;
        bool clipped = transCoord (fsize, dst, src, scale);
        if (clipped)
            scaleU = scale;
        else
            scaleL = scale;
    }
    return scaleL;
}

bool TransformFilter::transCoord (Dim fullSize, ImageView target, ImageView& source, double ascaleDef) {

    int x1 = target.x, y1 = target.y;
    int x2 = x1 + target.w - 1;
    int y2 = y1 + target.h - 1;

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

    bool result = transCoord (fullSize, corners, r, g, b, ascaleDef);

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

    source = ImageView (x1v, y1v, x2v - x1v + 1, y2v - y1v + 1);

    return result;
}


bool TransformFilter::transCoord (Dim fullSize, std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef) {

    bool clipresize = true;
    bool clipped = false;

    red.clear ();
    green.clear ();
    blue.clear ();

    if (!needsCA() && !needsDistortion() && !needsRotation() && !needsPerspective()) {
        for (int i=0; i<src.size(); i++) {
            red.push_back   (Coord2D (src[i].x, src[i].y));
            green.push_back (Coord2D (src[i].x, src[i].y));
            blue.push_back  (Coord2D (src[i].x, src[i].y));
        }
        return false;
    }

    int oW = fullSize.width;
    int oH = fullSize.height;
    double w2 = (double) oW  / 2.0 - 0.5;
    double h2 = (double) oH  / 2.0 - 0.5;
    double a = procParams->distortion.amount;
    double cost = cos(procParams->rotate.degree * 3.14/180.0);
    double sint = sin(procParams->rotate.degree * 3.14/180.0);
    double maxRadius = sqrt( (double)( oW*oW + oH*oH ) ) / 2;
    double vpdeg = procParams->perspective.vertical / 100.0 * 45.0;
    double vpalpha = (90.0 - vpdeg) / 180.0 * 3.14;
    double vpteta  = fabs(vpalpha-3.14/2)<1e-3 ? 0.0 : acos ((vpdeg>0 ? 1.0 : -1.0) * sqrt((-oW*oW*tan(vpalpha)*tan(vpalpha) + (vpdeg>0 ? 1.0 : -1.0) * oW*tan(vpalpha)*sqrt(16*maxRadius*maxRadius+oW*oW*tan(vpalpha)*tan(vpalpha)))/(maxRadius*maxRadius*8)));
    double vpcospt = (vpdeg>=0 ? 1.0 : -1.0) * cos (vpteta), vptanpt = tan (vpteta);
    double hpdeg = procParams->perspective.horizontal / 100.0 * 45.0;
    double hpalpha = (90.0 - hpdeg) / 180.0 * 3.14;
    double hpteta  = fabs(hpalpha-3.14/2)<1e-3 ? 0.0 : acos ((hpdeg>0 ? 1.0 : -1.0) * sqrt((-oH*oH*tan(hpalpha)*tan(hpalpha) + (hpdeg>0 ? 1.0 : -1.0) * oH*tan(hpalpha)*sqrt(16*maxRadius*maxRadius+oH*oH*tan(hpalpha)*tan(hpalpha)))/(maxRadius*maxRadius*8)));
    double hpcospt = (hpdeg>=0 ? 1.0 : -1.0) * cos (hpteta), hptanpt = tan (hpteta);

    double ascale = ascaleDef>0 ? ascaleDef : (procParams->commonTrans.autofill ? getTransformAutoFill () : 1.0);

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

        red.push_back   (Coord2D(Dx*(s+procParams->cacorrection.red)+w2, Dy*(s+procParams->cacorrection.red)+h2));
        green.push_back (Coord2D(Dx*s+w2, Dy*s+h2));
        blue.push_back  (Coord2D(Dx*(s+procParams->cacorrection.blue)+w2, Dy*(s+procParams->cacorrection.blue)+h2));
    }

    for (int i=0; i<src.size(); i++) {
        red[i].x   = CLIPTOC(red[i].x,0,oW-1,clipped);
        red[i].y   = CLIPTOC(red[i].y,0,oH-1,clipped);
        green[i].x = CLIPTOC(green[i].x,0,oW-1,clipped);
        green[i].y = CLIPTOC(green[i].y,0,oH-1,clipped);
        blue[i].x  = CLIPTOC(blue[i].x,0,oW-1,clipped);
        blue[i].y  = CLIPTOC(blue[i].y,0,oH-1,clipped);
    }

    return clipped;
}

bool TransformFilter::needsCA () {

    return fabs (procParams->cacorrection.red) > 1e-15 || fabs (procParams->cacorrection.blue) > 1e-15;
}
bool TransformFilter::needsDistortion () {

    return fabs (procParams->distortion.amount) > 1e-15;
}
bool TransformFilter::needsRotation () {

    return fabs (procParams->rotate.degree) > 1e-15;
}
bool TransformFilter::needsPerspective () {

    return procParams->perspective.horizontal || procParams->perspective.vertical;
}
bool TransformFilter::needsVignetting () {

    return procParams->vignetting.amount;
}
bool TransformFilter::needsTransform () {

    return needsCA () || needsDistortion () || needsRotation () || needsPerspective () || needsVignetting ();
}

#define A   (-0.85)
void TransformFilter::cubintch (float** src, int xs, int ys, double Dx, double Dy, float *r, double mul) {

	float w[4];

  { float t1, t2;
  t1 = -A*(Dx-1.0)*Dx;
  t2 = (3.0-2.0*Dx)*Dx*Dx;
  w[3] = t1*Dx;
  w[2] = t1*(Dx-1.0) + t2;
  w[1] = -t1*Dx + 1.0 - t2;
  w[0] = -t1*(Dx-1.0);
  }

  float rd;
  float yr[4];

  for (int k=ys, kx=0; k<ys+4; k++, kx++) {
    rd = 0.0;
    for (int i=xs, ix=0; i<xs+4; i++, ix++) {
      rd += src[k][i] * w[ix];
    }
    yr[kx] = rd;
  }


  { float t1, t2;
  t1 = -A*(Dy-1.0)*Dy;
  t2 = (3.0-2.0*Dy)*Dy*Dy;
  w[3] = t1*Dy;
  w[2] = t1*(Dy-1.0) + t2;
  w[1] = -t1*Dy + 1.0 - t2;
  w[0] = -t1*(Dy-1.0);
  }

  rd = 0.0;
  for (int i=0; i<4; i++)
    rd += yr[i] * w[i];

  rd*=mul;

  *r = rd;
}

void TransformFilter::cubint (MultiImage* src, int xs, int ys, double Dx, double Dy, float *r, float *g, float *b, double mul) {

  float w[4];

  { float t1, t2;
  t1 = -A*(Dx-1.0)*Dx;
  t2 = (3.0-2.0*Dx)*Dx*Dx;
  w[3] = t1*Dx;
  w[2] = t1*(Dx-1.0) + t2;
  w[1] = -t1*Dx + 1.0 - t2;
  w[0] = -t1*(Dx-1.0);
  }

  float rd, gd, bd;
  float yr[4], yg[4], yb[4];

  for (int k=ys, kx=0; k<ys+4; k++, kx++) {
    rd = gd = bd = 0.0;
    for (int i=xs, ix=0; i<xs+4; i++, ix++) {
      rd += src->r[k][i] * w[ix];
      gd += src->g[k][i] * w[ix];
      bd += src->b[k][i] * w[ix];
    }
    yr[kx] = rd; yg[kx] = gd; yb[kx] = bd;
  }


  { float t1, t2;
  t1 = -A*(Dy-1.0)*Dy;
  t2 = (3.0-2.0*Dy)*Dy*Dy;
  w[3] = t1*Dy;
  w[2] = t1*(Dy-1.0) + t2;
  w[1] = -t1*Dy + 1.0 - t2;
  w[0] = -t1*(Dy-1.0);
  }

  rd = gd = bd = 0.0;
  for (int i=0; i<4; i++) {
    rd += yr[i] * w[i];
    gd += yg[i] * w[i];
    bd += yb[i] * w[i];
  }

  rd*=mul;
  gd*=mul;
  bd*=mul;

  *r = rd;
  *g = gd;
  *b = bd;
}


}
