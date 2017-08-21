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
void ImProcFunctions::removeSpots (Imagefloat* img, const std::vector<SpotEntry> &entries, const PreviewProps &pp)
{
    Alpha mask;

    //printf("img(%04d, %04d)\n", img->width, img->height);

    for (const auto entry : entries) {
        float srcX = float (entry.sourcePos.x);
        float srcY = float (entry.sourcePos.y);
        float dstX = float (entry.targetPos.x);
        float dstY = float (entry.targetPos.y);
        //float radius = float (entry.radius) + 0.5f;

        float featherRadius = entry.radius * (1.f + entry.feather);
        int scaledFeatherRadius = featherRadius / pp.getSkip ();

        int src_XMin = int ((srcX - featherRadius - pp.getX()) / float (pp.getSkip()) + 0.5f);
        int src_XMax = int ((srcX + featherRadius - pp.getX()) / float (pp.getSkip()) + 0.5f);
        int src_YMin = int ((srcY - featherRadius - pp.getY()) / float (pp.getSkip()) + 0.5f);
        int src_YMax = int ((srcY + featherRadius - pp.getY()) / float (pp.getSkip()) + 0.5f);

        int dst_XMin = int ((dstX - featherRadius - pp.getX()) / float (pp.getSkip()) + 0.5f);
        int dst_XMax = int ((dstX + featherRadius - pp.getX()) / float (pp.getSkip()) + 0.5f);
        int dst_YMin = int ((dstY - featherRadius - pp.getY()) / float (pp.getSkip()) + 0.5f);
        int dst_YMax = int ((dstY + featherRadius - pp.getY()) / float (pp.getSkip()) + 0.5f);

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

        // skipping entries where the source circle isn't completely inside the image bounds
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

}

