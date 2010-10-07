/*
 * filtcoarse.cc
 *
 *  Created on: Sep 17, 2010
 *      Author: gabor
 */

#include "filtresize.h"
#include "rtengine.h"
#include "macros.h"
#include "filterchain.h"

namespace rtengine {

ResizeFilterDescriptor resizeFilterDescriptor;

ResizeFilterDescriptor::ResizeFilterDescriptor ()
	: FilterDescriptor ("Resize", MultiImage::RGB, MultiImage::RGB, true) {

	addTriggerEvent (EvResizeScale);
    addTriggerEvent (EvResizeMethod);
    addTriggerEvent (EvResizeSpec);
    addTriggerEvent (EvResizeWidth);
    addTriggerEvent (EvResizeHeight);
    addTriggerEvent (EvResizeEnabled);
}

void ResizeFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ResizeFilter ());
}

ResizeFilter::ResizeFilter ()
	: Filter (&resizeFilterDescriptor) {
}

ImageView ResizeFilter::calculateSourceImageView (const ImageView& requestedImView) {

    int x1, y1, x2, y2;
    reverseTransPoint (requestedImView.x, requestedImView.y, x1, y1);
    reverseTransPoint (requestedImView.x + requestedImView.w - 1, requestedImView.y + requestedImView.h - 1, x2, y2);

    return ImageView (std::min(x1,x2), std::min(y1,y2), ABS(x2-x1)+1, ABS(y2-y1)+1, 1);
}

double ResizeFilter::getScale () {

    double s = getParentFilter()->getScale();

    // if we are processing a thumbnail, do not apply resize filter, just update the "Scale", that is the
    // ratio of the image obtained compared to the image requested
    if (getFilterChain()->getImageSource()->isThumbnail())
        return s / getResizeScale ();
    else
        return s;
}

double ResizeFilter::getTargetScale (int skip) {

    double s = Filter::getTargetScale (skip);

    if (getFilterChain()->getImageSource()->isThumbnail())
        return s / getResizeScale ();
    else
        return s;
}

double ResizeFilter::getResizeScale () {

    if (procParams->resize.enabled) {
        Dim pdim = getPreviousFilter()->getFullImageSize ();
        if (procParams->resize.dataspec==1)
            return procParams->resize.width / pdim.width;
        else if (procParams->resize.dataspec==2)
            return procParams->resize.height / pdim.height;
        else if (procParams->resize.dataspec==0)
            return procParams->resize.scale;
        else
           return 1.0;
    }
    else
        return 1.0;
}

Dim ResizeFilter::getFullImageSize () {

    Dim pdim = getPreviousFilter()->getFullImageSize ();
    if (procParams->resize.enabled) {
        if (procParams->resize.dataspec==1)
            return Dim (procParams->resize.width, pdim.height * procParams->resize.width / pdim.width);
        else if (procParams->resize.dataspec==2)
            return Dim (pdim.width * procParams->resize.height / pdim.height);
        else if (procParams->resize.dataspec==0)
            return Dim (pdim.width * procParams->resize.scale, pdim.height * procParams->resize.scale);
        else
            return pdim;
    }
    else
        return pdim;
}

void ResizeFilter::reverseTransPoint (int x, int y, int& xv, int& yv) {

    xv = x / getResizeScale ();
    yv = y / getResizeScale ();
}

void ResizeFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    if (!getFilterChain()->getImageSource()->isThumbnail() && fabs(getResizeScale()-1.0) > 1e-12) {
        if (procParams->resize.method.substr(0,7)=="Bicubic")
            bicubic (sourceImage, targetImage);
        else if (procParams->resize.method.substr(0,7)=="Bilinear")
            bilinear (sourceImage, targetImage);
        else if (procParams->resize.method.substr(0,7)=="Average")
            average (sourceImage, targetImage);
        else
            nearest (sourceImage, targetImage);
    }
    else if (targetImage!=sourceImage)
        targetImage->copyFrom (sourceImage);
}

void ResizeFilter::bicubic (MultiImage* sourceImage, MultiImage* targetImage) {

    double scale = getResizeScale ();
    double Av = -0.5;
    if (procParams->resize.method=="Bicubic (Sharper)")
        Av = -0.75;
    else if (procParams->resize.method=="Bicubic (Softer)")
        Av = -0.25;
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<targetImage->height; i++) {
        double wx[4], wy[4];
        double Dy = i / scale;
        int yc  =  (int) Dy; Dy -= (double)yc;
        int ys = yc - 1; // smallest y-index used for interpolation
        // compute vertical weights
        double t1y = -Av*(Dy-1.0)*Dy;
        double t2y = (3.0-2.0*Dy)*Dy*Dy;
        wy[3] = t1y*Dy;
        wy[2] = t1y*(Dy-1.0) + t2y;
        wy[1] = -t1y*Dy + 1.0 - t2y;
        wy[0] = -t1y*(Dy-1.0);
        for (int j=0; j<targetImage->width; j++) {
            double Dx = j / scale;
            int xc  =  (int) Dx; Dx -= (double)xc;
            int xs = xc - 1; // smallest x-index used for interpolation
            if (ys >= 0 && ys <sourceImage->height-3 && xs >= 0 && xs <= sourceImage->width-3) {
                // compute horizontal weights
                double t1 = -Av*(Dx-1.0)*Dx;
                double t2 = (3.0-2.0*Dx)*Dx*Dx;
                wx[3] = t1*Dx;
                wx[2] = t1*(Dx-1.0) + t2;
                wx[1] = -t1*Dx + 1.0 - t2;
                wx[0] = -t1*(Dx-1.0);
                // compute weighted sum
                int r = 0;
                int g = 0;
                int b = 0;
                for (int x=0; x<4; x++)
                    for (int y=0; y<4; y++) {
                        double w = wx[x]*wy[y];
                        r += w*sourceImage->r[ys+y][xs+x];
                        g += w*sourceImage->g[ys+y][xs+x];
                        b += w*sourceImage->b[ys+y][xs+x];
                    }
                targetImage->r[i][j] = CLIP(r);
                targetImage->g[i][j] = CLIP(g);
                targetImage->b[i][j] = CLIP(b);
            }
            else {
                xc = CLIPTO(xc, 0, sourceImage->width-1);
                yc = CLIPTO(yc, 0, sourceImage->height-1);
                int nx = xc + 1;
                if (nx>=sourceImage->width)
                    nx = xc;
                int ny = yc + 1;
                if (ny>=sourceImage->height)
                    ny = yc;
                targetImage->r[i][j] = (1-Dx)*(1-Dy)*sourceImage->r[yc][xc] + (1-Dx)*Dy*sourceImage->r[ny][xc] + Dx*(1-Dy)*sourceImage->r[yc][nx] + Dx*Dy*sourceImage->r[ny][nx];
                targetImage->g[i][j] = (1-Dx)*(1-Dy)*sourceImage->g[yc][xc] + (1-Dx)*Dy*sourceImage->g[ny][xc] + Dx*(1-Dy)*sourceImage->g[yc][nx] + Dx*Dy*sourceImage->g[ny][nx];
                targetImage->b[i][j] = (1-Dx)*(1-Dy)*sourceImage->b[yc][xc] + (1-Dx)*Dy*sourceImage->b[ny][xc] + Dx*(1-Dy)*sourceImage->b[yc][nx] + Dx*Dy*sourceImage->b[ny][nx];
            }
        }
    }
}

void ResizeFilter::bilinear (MultiImage* sourceImage, MultiImage* targetImage) {

    double scale = getResizeScale ();
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<targetImage->height; i++) {
        int sy = i / scale;
        sy = CLIPTO(sy, 0, sourceImage->height-1);
        double dy = i / scale - sy;
        int ny = sy+1;
        if (ny>=sourceImage->height)
            ny = sy;
        for (int j=0; j<targetImage->width; j++) {
            int sx = j / scale;
            sx = CLIPTO(sx, 0, sourceImage->width-1);
            double dx = j / scale - sx;
            int nx = sx+1;
            if (nx>=sourceImage->width)
                nx = sx;
            targetImage->r[i][j] = (1-dx)*(1-dy)*sourceImage->r[sy][sx] + (1-dx)*dy*sourceImage->r[ny][sx] + dx*(1-dy)*sourceImage->r[sy][nx] + dx*dy*sourceImage->r[ny][nx];
            targetImage->g[i][j] = (1-dx)*(1-dy)*sourceImage->g[sy][sx] + (1-dx)*dy*sourceImage->g[ny][sx] + dx*(1-dy)*sourceImage->g[sy][nx] + dx*dy*sourceImage->g[ny][nx];
            targetImage->b[i][j] = (1-dx)*(1-dy)*sourceImage->b[sy][sx] + (1-dx)*dy*sourceImage->b[ny][sx] + dx*(1-dy)*sourceImage->b[sy][nx] + dx*dy*sourceImage->b[ny][nx];
        }
    }
}

void ResizeFilter::average (MultiImage* sourceImage, MultiImage* targetImage) {

    // small-scale algorithm by Ilia
    // provides much better quality on small scales
    // calculates mean value over source pixels which current destination pixel covers
    // works only for scales < 1
    // for scales ~1 it is analogous to bilinear
    // possibly, for even less scale factors (< 0.2 possibly) boundary pixels are not needed, omitting them can give a speedup
    // this algorithm is much slower on small factors than others, because it uses all pixels of the SOURCE image
    // Ilia Popov ilia_popov@rambler.ru 2010

    double scale = getResizeScale ();
    double delta = 1.0 / scale;
    double k = scale * scale;

    #pragma omp parallel for if (multiThread)
    for(int i = 0; i < targetImage->height; i++) {
        // top and bottom boundary coordinates
        double y0 = i * delta;
        double y1 = (i + 1) * delta;

        int m0 = y0;
        m0 = CLIPTO(m0, 0, sourceImage->height-1);

        int m1 = y1;
        m1 = CLIPTO(m1, 0, sourceImage->height-1);

        // weights of boundary pixels
        double wy0 = 1.0 - (y0 - m0);
        double wy1 = y1 - m1;

        for(int j = 0; j < targetImage->width; j++) {
            // left and right boundary coordinates
            double x0 = j * delta;
            double x1 = (j + 1) * delta;

            int n0 = x0;
            n0 = CLIPTO(n0, 0, sourceImage->width-1);
            int n1 = x1;
            n1 = CLIPTO(n1, 0, sourceImage->width-1);

            double wx0 = 1.0 - (x0 - n0);
            double wx1 = x1 - n1;

            double r = 0;
            double g = 0;
            double b = 0;

            // integration
            // corners
            r += wy0 * wx0 * sourceImage->r[m0][n0] + wy0 * wx1 * sourceImage->r[m0][n1] + wy1 * wx0 * sourceImage->r[m1][n0] + wy1 * wx1 * sourceImage->r[m1][n1];
            g += wy0 * wx0 * sourceImage->g[m0][n0] + wy0 * wx1 * sourceImage->g[m0][n1] + wy1 * wx0 * sourceImage->g[m1][n0] + wy1 * wx1 * sourceImage->g[m1][n1];
            b += wy0 * wx0 * sourceImage->b[m0][n0] + wy0 * wx1 * sourceImage->b[m0][n1] + wy1 * wx0 * sourceImage->b[m1][n0] + wy1 * wx1 * sourceImage->b[m1][n1];

            // top and bottom boundaries
            for(int n = n0 + 1; n < n1; n++) {
                r += wy0 * sourceImage->r[m0][n] + wy1 * sourceImage->r[m1][n];
                g += wy0 * sourceImage->g[m0][n] + wy1 * sourceImage->g[m1][n];
                b += wy0 * sourceImage->b[m0][n] + wy1 * sourceImage->b[m1][n];
            }

            // inner rows
            for(int m = m0 + 1; m < m1; m++) {
                // left and right boundaries
                r += wx0 * sourceImage->r[m][n0] + wx1 * sourceImage->r[m][n1];
                g += wx0 * sourceImage->g[m][n0] + wx1 * sourceImage->g[m][n1];
                b += wx0 * sourceImage->b[m][n0] + wx1 * sourceImage->b[m][n1];
                // inner pixels
                for(int n = n0 + 1; n < n1; n++) {
                    r += sourceImage->r[m][n];
                    g += sourceImage->g[m][n];
                    b += sourceImage->b[m][n];
                }
            }
            // overall weight is equal to the DST pixel area in SRC coordinates
            r *= k;
            g *= k;
            b *= k;

            targetImage->r[i][j] = CLIP((int)r);
            targetImage->g[i][j] = CLIP((int)g);
            targetImage->b[i][j] = CLIP((int)b);
        }
    }
}

void ResizeFilter::nearest (MultiImage* sourceImage, MultiImage* targetImage) {

    double scale = getResizeScale ();
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<targetImage->height; i++) {
        int sy = i / scale;
        sy = CLIPTO(sy, 0, sourceImage->height-1);
        for (int j=0; j<targetImage->width; j++) {
            int sx = j / scale;
            sx = CLIPTO(sx, 0, sourceImage->width-1);
            targetImage->r[i][j] = sourceImage->r[sy][sx];
            targetImage->g[i][j] = sourceImage->g[sy][sx];
            targetImage->b[i][j] = sourceImage->b[sy][sx];
        }
    }
}


}
