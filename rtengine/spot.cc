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

#include "improcfun.h"
#include "alpha.h"
#include "procparams.h"
#include "imagesource.h"

namespace rtengine
{

/* Code taken from Gimp 2.8.10 and converted for RawTherapee by Jean-Christophe FRISCH (aka Hombre) on 02.19.2014
 *
 * ORIGINAL NOTES
 *
 *     The method used here is similar to the lighting invariant correction
 *     method but slightly different: we do not divide the RGB components,
 *     but substract them I2 = I0 - I1, where I0 is the sample image to be
 *     corrected, I1 is the reference pattern. Then we solve DeltaI=0
 *     (Laplace) with I2 Dirichlet conditions at the borders of the
 *     mask. The solver is a unoptimized red/black checker Gauss-Siedel
 *     with an over-relaxation factor of 1.8. It can benefit from a
 *     multi-grid evaluation of an initial solution before the main
 *     iteration loop.
 *
 *     I reduced the convergence criteria to 0.1% (0.001) as we are
 *     dealing here with RGB integer components, more is overkill.
 *
 *     Jean-Yves Couleaud cjyves@free.fr
 */

/* Original Algorithm Design:
 *
 * T. Georgiev, "Photoshop Healing Brush: a Tool for Seamless Cloning
 * http://www.tgeorgiev.net/Photoshop_Healing.pdf
 */

#if 1

class SpotBox {

public:
    enum class Type {
        SOURCE,
        TARGET
    };

private:
    Type type;

public:
    int topLeftX;
    int topLeftY;
    int bottomRightX;
    int bottomRightY;
    Imagefloat* img;

    SpotBox (int tl_x, int tl_y, int br_x, int br_y, Type type) :
       type(type),
       topLeftX(tl_x),
       topLeftY(tl_y),
       bottomRightX(br_x),
       bottomRightY(br_y),
       img(nullptr)
    {}

    SpotBox (int tl_x, int tl_y, Imagefloat* image, Type type) :
       type(type),
       topLeftX(tl_x),
       topLeftY(tl_y),
       bottomRightX(image ? tl_x + image->getWidth() - 1 : 0),
       bottomRightY(image ? tl_y + image->getHeight() - 1 : 0),
       img(image)
    {}

    SpotBox (SpotEntry &spot, Type type) :
        type(type),
        img(nullptr)
    {
        float featherRadius = spot.radius * (1.f + spot.feather);
        topLeftX = int ((type == Type::SOURCE ? spot.sourcePos.x : spot.targetPos.x) - featherRadius);
        bottomRightX = int ((type == Type::SOURCE ? spot.sourcePos.x : spot.targetPos.x) + featherRadius);
        topLeftY = int ((type == Type::SOURCE ? spot.sourcePos.y : spot.targetPos.y) - featherRadius);
        bottomRightY = int ((type == Type::SOURCE ? spot.sourcePos.y : spot.targetPos.y) + featherRadius);
    }

    void translate(int dx, int dy) {
        topLeftX += dx;
        topLeftY += dy;
        bottomRightX += dx;
        bottomRightY += dy;
    }

    void operator /(float v) {
        topLeftX = int(topLeftX / v + 0.5f);
        topLeftY = int(topLeftY / v + 0.5f);
        bottomRightX = int(bottomRightX / v + 0.5f);
        bottomRightY = int(bottomRightY / v + 0.5f);
    }

    void operator *(float v) {
        topLeftX *= v;
        topLeftY *= v;
        bottomRightX *= v;
        bottomRightY *= v;
    }

    bool intersects(const SpotBox &other) const {
        return (other.topLeftX <= bottomRightX && other.bottomRightX >= topLeftX)
            && (other.topLeftY <= bottomRightY && other.bottomRightY >= topLeftY);
    }

    int getWidth() {
        return bottomRightX - topLeftX + 1;
    }

    int getHeight() {
        return bottomRightY - topLeftY + 1;
    }
};


void ImProcFunctions::removeSpots (Imagefloat* img, ImageSource* imgsrc, const std::vector<SpotEntry> &entries, const PreviewProps &pp, const ColorTemp &currWB, int tr)
{
    // ---------- Get the image areas (src & dst) from the source image

    printf("\n=======================================================================\n\n");

    std::vector< std::shared_ptr<SpotBox> > srcSpotBoxs;
    std::vector< std::shared_ptr<SpotBox> > dstSpotBoxs;
    for (auto entry : params->spot.entries) {
        Coord origin;
        int size = int(entry.getFeatherRadius() * 2.f + 0.5f);
        int scaledSize = int(entry.getFeatherRadius() * 2.f / float(pp.getSkip()) + 0.5f);
        //printf("size: %d - skip: %d  -> scaledSize: %d", size, pp.getSkip(), scaledSize);

        // ------ Source area
        Imagefloat *currSrcSpot = new Imagefloat(scaledSize, scaledSize);
        for (int y = 0; y < currSrcSpot->getHeight(); ++y) {
            for (int x = 0; x < currSrcSpot->getWidth(); ++x) {
                currSrcSpot->r(y,x) = 0.f;
                currSrcSpot->g(y,x) = 0.f;
                currSrcSpot->b(y,x) = 0.f;
            }
        }
        entry.sourcePos.get(origin.x, origin.y);
        origin.x -= entry.getFeatherRadius();
        origin.y -= entry.getFeatherRadius();
        PreviewProps spp(origin.x, origin.y, size, size, pp.getSkip());
        imgsrc->getImage(currWB, tr, currSrcSpot, spp, params->toneCurve, params->raw);
        //printf("  /  src size: %d,%d", currSrcSpot->getWidth(), currSrcSpot->getHeight());

        std::shared_ptr<SpotBox> srcSpotBox(new SpotBox(origin.x / pp.getSkip(), origin.y / pp.getSkip(), currSrcSpot,  SpotBox::Type::SOURCE));
        srcSpotBoxs.push_back(srcSpotBox);

        // ------ Destination area
        Imagefloat *currDstSpot = new Imagefloat(scaledSize, scaledSize);
        for (int y = 0; y < currDstSpot->getHeight(); ++y) {
            for (int x = 0; x < currDstSpot->getWidth(); ++x) {
                currDstSpot->r(y,x) = 0.f;
                currDstSpot->g(y,x) = 0.f;
                currDstSpot->b(y,x) = 0.f;
            }
        }
        entry.targetPos.get(origin.x, origin.y);
        origin.x -= entry.getFeatherRadius();
        origin.y -= entry.getFeatherRadius();
        spp.set(origin.x, origin.y, size, size, pp.getSkip());
        imgsrc->getImage(currWB, tr, currDstSpot, spp, params->toneCurve, params->raw);
        //printf("  /  dst size: %d,%d\n", currDstSpot->getWidth(), currDstSpot->getHeight());

        std::shared_ptr<SpotBox> dstSpotBox(new SpotBox(origin.x / pp.getSkip(), origin.y / pp.getSkip(), currDstSpot,  SpotBox::Type::TARGET));

        dstSpotBoxs.push_back(dstSpotBox);
    }

    // Filter out out of preview Spots

    /*
    for (size_t i = entries.size(); i >= 0; ++i) {
        float featherRadius = entries.at(i).radius * (1.f + entries.at(i).feather);

        SpotBox srcBox(entries.at(i), SpotBox::Type::SOURCE);
        srcBox.translate(-pp.getX(), -pp.getY());
        srcBox /= float (pp.getSkip());

        SpotBox dstBox(entries.at(i), SpotBox::Type::TARGET);
        dstBox.translate(-pp.getX(), -pp.getY());
        dstBox /= float (pp.getSkip());

    }
    */



    // ---------- Copy spots from src to dst

    for (int i = entries.size() - 1; i >= 0; --i) {
        // 1. copy src to dst
        std::shared_ptr<SpotBox> srcSpotBox = srcSpotBoxs.at(i);
        std::shared_ptr<SpotBox> dstSpotBox = dstSpotBoxs.at(i);
        float scaledRadius = float(entries.at(i).radius) / float(pp.getSkip());
        float scaledFeatherRadius = entries.at(i).getFeatherRadius() / float(pp.getSkip());
        Imagefloat *srcImg = srcSpotBox->img;
        Imagefloat *dstImg = dstSpotBox->img;

        //printf("#%d: srcSpotBox @ %p - img @ %p  /  dstSpotBox @ %p - img @ %p\n", i,
        //        srcSpotBox.get(), srcSpotBox->img, dstSpotBox.get(), dstSpotBox->img );

        //printf("#%d: srcSpotBox(%d,%d) srcImg(%d,%d)  /  dstSpotBox(%d,%d) dstImg(%d,%d)\n", i,
        //        srcSpotBox->getWidth(), srcSpotBox->getHeight(), srcImg->getWidth(), srcImg->getHeight(),
        //        dstSpotBox->getWidth(), dstSpotBox->getHeight(), dstImg->getWidth(), dstImg->getHeight()
        //        );

        for (int y = 0; y < srcSpotBox->getHeight(); ++y) {
            float  dy = float(y - float(srcSpotBox->getHeight()) / 2.f);
            for (int x = 0; x < srcSpotBox->getWidth(); ++x) {
                float dx = float(x - float(srcSpotBox->getWidth()) / 2.f);
                float r = sqrt(dx * dx + dy * dy);
                if (r >= scaledFeatherRadius) {
                    continue;
                }
                if (r <= scaledRadius) {
                    dstImg->r(y, x) = srcImg->r(y, x);
                    dstImg->g(y, x) = srcImg->g(y, x);
                    dstImg->b(y, x) = srcImg->b(y, x);
                } else {
                    float opacity = (scaledFeatherRadius - r) / (scaledFeatherRadius - scaledRadius);
                    dstImg->r(y, x) = (srcImg->r(y, x) - dstImg->r(y, x)) * opacity + dstImg->r(y,x);
                    dstImg->g(y, x) = (srcImg->g(y, x) - dstImg->g(y, x)) * opacity + dstImg->g(y,x);
                    dstImg->b(y, x) = (srcImg->b(y, x) - dstImg->b(y, x)) * opacity + dstImg->b(y,x);
                }
            }
        }
        //printf("\n\n");

        // 2. copy dst to later src and dst

    }

    // 3. copy all dst to the finale image

    // Putting the dest image in a SpotBox
    SpotBox imgSpotBox(pp.getX() / pp.getSkip(), pp.getY() / pp.getSkip(), img, SpotBox::Type::TARGET);
    /*
    printf("#--: spotBox(X1:%d, Y1:%d, X2:%d, Y2:%d, W:%d, H:%d)  img(W:%d, H:%d)\n\n",
            imgSpotBox.topLeftX, imgSpotBox.topLeftY, imgSpotBox.bottomRightX, imgSpotBox.bottomRightY,
            imgSpotBox.getWidth(), imgSpotBox.getHeight(),
            imgSpotBox.img->getWidth(), imgSpotBox.img->getHeight()
            );
    */

    for (size_t i = 0; i < entries.size(); ++i) {
        // 1. copy src to dst
        std::shared_ptr<SpotBox> dstSpotBox = dstSpotBoxs.at(i);
        Imagefloat *dstImg = dstSpotBox->img;

        /*
        printf("#%llu: spotBox(X1:%d, Y1:%d, X2:%d, Y2:%d, W:%d, H:%d)  img(W:%d, H:%d)\n", i,
                dstSpotBox->topLeftX, dstSpotBox->topLeftY, dstSpotBox->bottomRightX, dstSpotBox->bottomRightY,
                dstSpotBox->getWidth(), dstSpotBox->getHeight(),
                dstImg->getWidth(), dstImg->getHeight()
                );
        */

        if (dstSpotBox->intersects(imgSpotBox)) {
            int beginX = rtengine::max(dstSpotBox->topLeftX, imgSpotBox.topLeftX);
            int endX = rtengine::min(dstSpotBox->bottomRightX, imgSpotBox.bottomRightX);
            int beginY = rtengine::max(dstSpotBox->topLeftY, imgSpotBox.topLeftY);
            int endY = rtengine::min(dstSpotBox->bottomRightY, imgSpotBox.bottomRightY);

            //printf("--- Intersection:  X1:%d, Y1:%d -> X2:%d, Y2:%d\n", beginX, beginY, endX, endY);

            int dstSpotOffsetY = beginY - dstSpotBox->topLeftY;
            int imgOffsetY = beginY - imgSpotBox.topLeftY;

            for (int y = beginY; y <= endY; ++y) {
                int dstSpotOffsetX = beginX - dstSpotBox->topLeftX;
                int imgOffsetX = beginX - imgSpotBox.topLeftX;

                for (int x = beginX; x <= endX; ++x) {
                    /*
                    if (y == beginY && x == beginX) {
                        printf("--- dstSpotOffsetX = beginX - dstSpotBox->topLeftX = %d - %d = %d\n", beginX, dstSpotBox->topLeftX, dstSpotOffsetX);
                        printf("--- dstSpotOffsetY = beginY - dstSpotBox->topLeftY = %d - %d = %d\n", beginY, dstSpotBox->topLeftY, dstSpotOffsetY);
                        printf("--- imgOffsetX = beginX - imgSpotBox.topLeftX = %d - %d = %d\n", beginX, imgSpotBox.topLeftX, imgOffsetX);
                        printf("--- imgOffsetX = beginY - imgSpotBox.topLeftY = %d - %d = %d\n", beginY, imgSpotBox.topLeftY, imgOffsetY);
                    }
                    */
                    img->r(imgOffsetY, imgOffsetX) = dstImg->r(dstSpotOffsetY, dstSpotOffsetX);
                    img->g(imgOffsetY, imgOffsetX) = dstImg->g(dstSpotOffsetY, dstSpotOffsetX);
                    img->b(imgOffsetY, imgOffsetX) = dstImg->b(dstSpotOffsetY, dstSpotOffsetX);
                    ++imgOffsetX;
                    ++dstSpotOffsetX;
                }
                ++imgOffsetY;
                ++dstSpotOffsetY;
            }
        //} else {
        //    printf("#%llu: No intersection !\n", i);
        }
    }

    for (auto srcSpotBox : srcSpotBoxs) {
        delete srcSpotBox->img;
    }
    for (auto dstSpotBox : dstSpotBoxs) {
        delete dstSpotBox->img;
    }
}

#endif





#if 0
void ImProcFunctions::removeSpots (Imagefloat* img, ImageSource* imgsrc, const std::vector<SpotEntry> &entries, const PreviewProps &pp, const ColorTemp &currWB, int tr)
{

    Alpha mask;


    for (const auto entry : entries) {
        float featherRadius = entry.getFeatherRadius();
        int scaledFeatherRadius = featherRadius / pp.getSkip ();

        SpotBox srcBox(entry, SpotBox::Type::SOURCE);
        srcBox.translate(-pp.getX(), -pp.getY());
        srcBox /= float (pp.getSkip());

        SpotBox dstBox(entry, SpotBox::Type::TARGET);
        dstBox.translate(-pp.getX(), -pp.getY());
        dstBox /= float (pp.getSkip());

        //printf("  -> X: %04d > %04d\n  -> Y: %04d > %04d\n", dst_XMin, dst_XMax, dst_YMin, dst_YMax);

        // scaled spot is too small, we do not preview it
        if (scaledFeatherRadius < 2 && pp.getSkip() != 1) {
#ifndef NDEBUG

            if (options.rtSettings.verbose) {
                printf ("Skipping spot located at %d x %d, too small for the preview zoom rate\n", entry.sourcePos.x, entry.sourcePos.y);
            }

#endif
            continue;
        }

        // skipping entries totally transparent
        if (entry.opacity == 0.) {
#ifndef NDEBUG

            if (options.rtSettings.verbose) {
                printf ("Skipping spot located at %d x %d: opacity=%.3f\n", entry.sourcePos.x, entry.sourcePos.y, entry.opacity);
            }

            continue;
#endif
        }

        // skipping entries where the source circle isn't inside the image bounds, even partially
        if (src_XMin < 0 || src_XMax >= img->getWidth() || src_YMin < 0 || src_YMax >= img->getHeight()) {
#ifndef NDEBUG

            if (options.rtSettings.verbose) {
                printf ("Skipping spot located at %d x %d, from the data at %d x %d, radius=%d, feather=%.3f, opacity=%.3f: source out of bounds\n", entry.sourcePos.x, entry.sourcePos.y, entry.targetPos.x, entry.targetPos.y, entry.radius, entry.feather, entry.opacity);
                printf ("%d < 0 || %d >= %d || %d < 0 || %d >= %d\n",
                        src_XMin, src_XMax, img->getWidth(), src_YMin, src_YMax, img->getHeight());
            }

#endif
            continue;
        }

        // skipping entries where the dest circle is completely outside the image bounds
        /*
        if (dst_XMin >= img->getWidth() || dst_XMax <= 0 || dst_YMin >= img->getHeight() || dst_YMax <= 0) {
#ifndef NDEBUG

            if (options.rtSettings.verbose) {
                printf ("Skipping spot located at %d x %d, from the data at %d x %d, radius=%d, feather=%.3f, opacity=%.3f: source out of bounds\n", entry.sourcePos.x, entry.sourcePos.y, entry.targetPos.x, entry.targetPos.y, entry.radius, entry.feather, entry.opacity);
                printf ("%d >= %d || %d <= 0 || %d >= %d || %d <= 0\n",
                        dst_XMin, img->getWidth(), dst_XMax, dst_YMin, img->getHeight(), dst_YMax);
            }

#endif
            continue;
        }
        */

        // ----------------- Core function -----------------

#if 0
        int scaledPPX = pp.getX() / pp.skip;
        int scaledPPY = pp.getY() / pp.skip;
        int scaledPPW = pp.getWidth() / pp.skip + (pp.getWidth() % pp.getSkip() > 0);
        int scaledPPH = pp.getHeight() / pp.skip + (pp.getHeight() % pp.skip > 0);

        int sizeX = dst_XMax - dst_XMin + 1;
        int sizeY = dst_YMax - dst_YMin + 1;

        Imagefloat matrix (sizeX, sizeY);
        Imagefloat solution (sizeX, sizeY);

        // allocate the mask and draw it
        mask.setSize (sizeX, sizeY);
        {
            Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create (mask.getSurface());

            // clear the bitmap
            cr->set_source_rgba (0., 0., 0., 0.);
            cr->rectangle (0., 0., sizeX, sizeY);
            cr->set_line_width (0.);
            cr->fill();

            // draw the mask
            cr->set_antialias (Cairo::ANTIALIAS_GRAY);
            cr->set_line_width (featherRadius);
            double gradientCenterX = double (sizeX) / 2.;
            double gradientCenterY = double (sizeY) / 2.;
            {
                Cairo::RefPtr<Cairo::RadialGradient> radialGradient = Cairo::RadialGradient::create (
                            gradientCenterX, gradientCenterY, radius,
                            gradientCenterX, gradientCenterY, featherRadius
                        );
                radialGradient->add_color_stop_rgb (0., 0., 0., 1.);
                radialGradient->add_color_stop_rgb (1., 0., 0., 0.);
                cr->set_source_rgba (0., 0., 0., 1.);
                cr->mask (radialGradient);
                cr->rectangle (0., 0., sizeX, sizeY);
                cr->fill();
            }
        }

        // copy the src part to a temporary buffer to avoid possible self modified source
        Imagefloat *srcBuff = img->copySubRegion (srcX, srcY, sizeX, sizeY);


        // subtract pattern to image and store the result as a double in matrix
        for (int i = 0, i2 = dst_YMin;     i2 < sizeY - 1;     ++i, ++i2) {
            for (int j = 0, j2 = dst_XMin;     i2 < sizeX - 1;     ++j, ++j2) {
                matrix.r (i, j) = img->r (i2, j2) - srcBuff->r (i, j);
                matrix.g (i, j) = img->g (i2, j2) - srcBuff->g (i, j);
                matrix.b (i, j) = img->b (i2, j2) - srcBuff->b (i, j);
            }
        }


        // FIXME: is a faster implementation needed?
#define EPSILON   0.001
#define MAX_ITER  500

        // repeat until convergence or max iterations
        for (int n = 0; n < MAX_ITER; ++n) {

            printf ("<<< n=#%d\n", n);
            // ----------------------------------------------------------------

            /* Perform one iteration of the Laplace solver for matrix.  Store the
             * result in solution and get the square of the cumulative error
             * of the solution.
             */
            int i, j;
            double tmp, diff;
            double sqr_err_r = 0.0;
            double sqr_err_g = 0.0;
            double sqr_err_b = 0.0;
            const double w = 1.80 * 0.25; /* Over-relaxation = 1.8 */

            // we use a red/black checker model of the discretization grid

            // do reds
            for (i = 0; i < matrix.getHeight(); ++i) {
                for (j = i % 2; j < matrix.getWidth(); j += 2) {
                    printf ("/%d,%d", j, i);

                    if ((0 == mask (i, j)) || (i == 0) || (i == (matrix.getHeight() - 1)) || (j == 0) || (j == (matrix.getWidth() - 1))) {
                        // do nothing at the boundary or outside mask
                        solution.r (i, j) = matrix.r (i, j);
                        solution.g (i, j) = matrix.g (i, j);
                        solution.b (i, j) = matrix.b (i, j);
                    } else {
                        // Use Gauss Siedel to get the correction factor then over-relax it
                        tmp = solution.r (i, j);
                        solution.r (i, j) = (matrix.r (i, j) + w *
                                             (
                                                 matrix.r (i,  j - 1) +  // west
                                                 matrix.r (i,  j + 1) +  // east
                                                 matrix.r (i - 1,  j) +  // north
                                                 matrix.r (i + 1,  j) - 4.0 * matrix.r (i, j) // south
                                             )
                                            );

                        diff = solution.r (i, j) - tmp;
                        sqr_err_r += diff * diff;


                        tmp = solution.g (i, j);
                        solution.g (i, j) = (matrix.g (i, j) + w *
                                             (
                                                 matrix.g (i,  j - 1) +  // west
                                                 matrix.g (i,  j + 1) +  // east
                                                 matrix.g (i - 1,  j) +  // north
                                                 matrix.g (i + 1,  j) - 4.0 * matrix.g (i, j) // south
                                             )
                                            );

                        diff = solution.g (i, j) - tmp;
                        sqr_err_g += diff * diff;



                        tmp = solution.b (i, j);
                        solution.b (i, j) = (matrix.b (i, j) + w *
                                             (
                                                 matrix.b (i,  j - 1) +  // west
                                                 matrix.b (i,  j + 1) +  // east
                                                 matrix.b (i - 1,  j) +  // north
                                                 matrix.b (i + 1,  j) - 4.0 * matrix.b (i, j) // south
                                             )
                                            );

                        diff = solution.b (i, j) - tmp;
                        sqr_err_b += diff * diff;

                    }
                }
            }


            /* Do blacks
            *
            * As we've done the reds earlier, we can use them right now to
            * accelerate the convergence. So we have "solution" in the solver
            * instead of "matrix" above
            */
            for (i = 0; i < matrix.getHeight(); i++) {
                for (j = (i % 2) ? 0 : 1; j < matrix.getWidth(); j += 2) {
                    printf (":%d,%d", j, i);

                    if ((0 == mask (i, j)) || (i == 0) || (i == (matrix.getHeight() - 1)) || (j == 0) || (j == (matrix.getWidth() - 1))) {
                        // do nothing at the boundary or outside mask
                        solution.r (i, j) = matrix.r (i, j);
                        solution.g (i, j) = matrix.g (i, j);
                        solution.b (i, j) = matrix.b (i, j);
                    } else {
                        // Use Gauss Siedel to get the correction factor then over-relax it
                        tmp = solution.r (i, j);
                        solution.r (i, j) = (matrix.r (i, j) + w *
                                             (
                                                 matrix.r (i,  j - 1) +  // west
                                                 matrix.r (i,  j + 1) +  // east
                                                 matrix.r (i - 1,  j) +  // north
                                                 matrix.r (i + 1,  j) - 4.0 * matrix.r (i, j) // south
                                             )
                                            );

                        diff = solution.r (i, j) - tmp;
                        sqr_err_r += diff * diff;



                        tmp = solution.g (i, j);
                        solution.g (i, j) = (matrix.g (i, j) + w *
                                             (
                                                 matrix.g (i,  j - 1) +  // west
                                                 matrix.g (i,  j + 1) +  // east
                                                 matrix.g (i - 1,  j) +  // north
                                                 matrix.g (i + 1,  j) - 4.0 * matrix.g (i, j) // south
                                             )
                                            );

                        diff = solution.g (i, j) - tmp;
                        sqr_err_g += diff * diff;



                        tmp = solution.b (i, j);
                        solution.b (i, j) = (matrix.b (i, j) + w *
                                             (
                                                 matrix.b (i,  j - 1) +  // west
                                                 matrix.b (i,  j + 1) +  // east
                                                 matrix.b (i - 1,  j) +  // north
                                                 matrix.b (i + 1,  j) - 4.0 * matrix.b (i, j) // south
                                             )
                                            );

                        diff = solution.b (i, j) - tmp;
                        sqr_err_b += diff * diff;
                    }
                }
            }

            // ----------------------------------------------------------------

            // copy solution to matrix
            solution.copyData (&matrix);

            if (sqr_err_r < EPSILON && sqr_err_g < EPSILON && sqr_err_b < EPSILON) {
                break;
            }

            printf ("\n>>> n=#%d\n", n);
        }

        printf ("\n");
#endif




        // add solution to original image and store in tempPR
        for     (int i = 0, i2 = dst_YMin; i2 < dst_YMax - 1; ++i, ++i2) {
            if (i2 < 0 || i2 >= img->getHeight()) {
                continue;
            }

            for (int j = 0, j2 = dst_XMin; j2 < dst_XMax - 1; ++j, ++j2) {
                if (j2 < 0 || j2 >= img->getWidth()) {
                    continue;
                }

                //float c2 = float (mask (i, j)) / 255.f;
                //float c1 = 1.f - c2;
                //resultPR->r(i,j) = (unsigned char) CLAMP0255 ( ROUND( double(first->r(i,j)) + double(secondPR->r(i,j)) ) );


                img->r (i2, j2) = 65535.0f; //img->r(i2,j2)*c1 + srcBuff->r(i,j)*c2;
                img->g (i2, j2) = 0.0f; //img->g(i2,j2)*c1 + srcBuff->g(i,j)*c2;
                img->b (i2, j2) = 0.0f; //img->b(i2,j2)*c1 + srcBuff->b(i,j)*c2;
                /*
                img->r(i2,j2) = img->r(i2,j2)*c1 + (solution.r(i,j) + srcBuff->r(i,j))*c2;
                img->g(i2,j2) = img->g(i2,j2)*c1 + (solution.g(i,j) + srcBuff->g(i,j))*c2;
                img->b(i2,j2) = img->b(i2,j2)*c1 + (solution.b(i,j) + srcBuff->b(i,j))*c2;
                */
            }
        }

    }
}
#endif

}

