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
#include <cmath>
#include <iostream>

#include "rtengine.h"
#include "rawimagesource.h"
#include "rawimagesource_i.h"
#include "jaggedarray.h"
#include "median.h"
#include "rawimage.h"
#include "mytime.h"
#include "iccstore.h"
#include "curves.h"
#include "dfmanager.h"
#include "ffmanager.h"
#include "dcp.h"
#include "rt_math.h"
#include "improcfun.h"
#include "rtlensfun.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "opthelper.h"
#define BENCHMARK
#include "StopWatch.h"
#define clipretinex( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )
#undef CLIPD
#define CLIPD(a) ((a)>0.0f?((a)<1.0f?(a):1.0f):0.0f)

namespace
{

void rotateLine (const float* const line, rtengine::PlanarPtr<float> &channel, const int tran, const int i, const int w, const int h)
{
    switch(tran & TR_ROT) {
        case TR_R180:
            for (int j = 0; j < w; j++) {
                channel(h - 1 - i, w - 1 - j) = line[j];
            }

            break;

        case TR_R90:
            for (int j = 0; j < w; j++) {
                channel(j, h - 1 - i) = line[j];
            }

            break;

        case TR_R270:
            for (int j = 0; j < w; j++) {
                channel(w - 1 - j, i) = line[j];
            }

            break;

        case TR_NONE:
        default:
            for (int j = 0; j < w; j++) {
                channel(i, j) = line[j];
            }
    }
}

void transLineStandard (const float* const red, const float* const green, const float* const blue, const int i, rtengine::Imagefloat* const image, const int tran, const int imwidth, const int imheight)
{
    // conventional CCD coarse rotation
    rotateLine (red, image->r, tran, i, imwidth, imheight);
    rotateLine (green, image->g, tran, i, imwidth, imheight);
    rotateLine (blue, image->b, tran, i, imwidth, imheight);
}

void transLineFuji (const float* const red, const float* const green, const float* const blue, const int i, rtengine::Imagefloat* const image, const int tran, const int imheight, const int fw)
{

    // Fuji SuperCCD rotation + coarse rotation
    int start = ABS(fw - i);
    int w = fw * 2 + 1;
    int h = (imheight - fw) * 2 + 1;
    int end = min(h + fw - i, w - fw + i);

    switch(tran & TR_ROT) {
        case TR_R180:
            for (int j = start; j < end; j++) {
                int y = i + j - fw;
                int x = fw - i + j;

                if (x >= 0 && y < image->getHeight() && y >= 0 && x < image->getWidth()) {
                    image->r(image->getHeight() - 1 - y, image->getWidth() - 1 - x) = red[j];
                    image->g(image->getHeight() - 1 - y, image->getWidth() - 1 - x) = green[j];
                    image->b(image->getHeight() - 1 - y, image->getWidth() - 1 - x) = blue[j];
                }
            }

            break;

        case TR_R270:
            for (int j = start; j < end; j++) {
                int y = i + j - fw;
                int x = fw - i + j;

                if (x >= 0 && x < image->getHeight() && y >= 0 && y < image->getWidth()) {
                    image->r(image->getHeight() - 1 - x, y) = red[j];
                    image->g(image->getHeight() - 1 - x, y) = green[j];
                    image->b(image->getHeight() - 1 - x, y) = blue[j];
                }
            }

            break;

        case TR_R90:
            for (int j = start; j < end; j++) {
                int y = i + j - fw;
                int x = fw - i + j;

                if (x >= 0 && y < image->getWidth() && y >= 0 && x < image->getHeight()) {
                    image->r(x, image->getWidth() - 1 - y) = red[j];
                    image->g(x, image->getWidth() - 1 - y) = green[j];
                    image->b(x, image->getWidth() - 1 - y) = blue[j];
                }
            }

            break;

        case TR_NONE:
        default:
            for (int j = start; j < end; j++) {
                int y = i + j - fw;
                int x = fw - i + j;

                if (x >= 0 && y < image->getHeight() && y >= 0 && x < image->getWidth()) {
                    image->r(y, x) = red[j];
                    image->g(y, x) = green[j];
                    image->b(y, x) = blue[j];
                }
            }
    }
}

void transLineD1x (const float* const red, const float* const green, const float* const blue, const int i, rtengine::Imagefloat* const image, const int tran, const int imwidth, const int imheight, const bool oddHeight, const bool clip)
{
    // Nikon D1X has an uncommon sensor with 4028 x 1324 sensels.
    // Vertical sensel size is 2x horizontal sensel size
    // We have to do vertical interpolation for the 'missing' rows
    // We do that in combination with coarse rotation

    switch(tran & TR_ROT) {
        case TR_R180: // rotate 180 degree
            for (int j = 0; j < imwidth; j++) {
                image->r(2 * (imheight - 1 - i), imwidth - 1 - j) = red[j];
                image->g(2 * (imheight - 1 - i), imwidth - 1 - j) = green[j];
                image->b(2 * (imheight - 1 - i), imwidth - 1 - j) = blue[j];
            }

            if (i == 0) {
                for (int j = 0; j < imwidth; j++) {
                    image->r(2 * imheight - 1, imwidth - 1 - j) = red[j];
                    image->g(2 * imheight - 1, imwidth - 1 - j) = green[j];
                    image->b(2 * imheight - 1, imwidth - 1 - j) = blue[j];
                }
            }

            if (i == 1 || i == 2) { // linear interpolation
                int row = 2 * imheight - 1 - 2 * i;

                for (int j = 0; j < imwidth; j++) {
                    int col = imwidth - 1 - j;
                    image->r(row, col) = (red[j] + image->r(row + 1, col)) / 2;
                    image->g(row, col) = (green[j] + image->g(row + 1, col)) / 2;
                    image->b(row, col) = (blue[j] + image->b(row + 1, col)) / 2;
                }

                if(i == 2 && oddHeight) {
                    int row = 2 * imheight;

                    for (int j = 0; j < imwidth; j++) {
                        int col = imwidth - 1 - j;
                        image->r(row, col) = (red[j] + image->r(row - 2, col)) / 2;
                        image->g(row, col) = (green[j] + image->g(row - 2, col)) / 2;
                        image->b(row, col) = (blue[j] + image->b(row - 2, col)) / 2;
                    }
                }
            } else if (i == imheight - 1 || i == imheight - 2) {
                int row = 2 * imheight - 1 - 2 * i;

                for (int j = 0; j < imwidth; j++) {
                    int col = imwidth - 1 - j;
                    image->r(row, col) = (red[j] + image->r(row + 1, col)) / 2;
                    image->g(row, col) = (green[j] + image->g(row + 1, col)) / 2;
                    image->b(row, col) = (blue[j] + image->b(row + 1, col)) / 2;
                }

                row = 2 * imheight - 1 - 2 * i + 2;

                for (int j = 0; j < imwidth; j++) {
                    int col = imwidth - 1 - j;
                    image->r(row, col) = (red[j] + image->r(row + 1, col)) / 2;
                    image->g(row, col) = (green[j] + image->g(row + 1, col)) / 2;
                    image->b(row, col) = (blue[j] + image->b(row + 1, col)) / 2;
                }
            } else if (i > 2 && i < imheight - 1) { // vertical bicubic interpolation
                int row = 2 * imheight - 1 - 2 * i + 2;

                for (int j = 0; j < imwidth; j++) {
                    int col = imwidth - 1 - j;
                    image->r(row, col) = MAX(0.f, -0.0625f * (red[j] + image->r(row + 3, col)) + 0.5625f * (image->r(row - 1, col) + image->r(row + 1, col)));
                    image->g(row, col) = MAX(0.f, -0.0625f * (green[j] + image->g(row + 3, col)) + 0.5625f * (image->g(row - 1, col) + image->g(row + 1, col)));
                    image->b(row, col) = MAX(0.f, -0.0625f * (blue[j] + image->b(row + 3, col)) + 0.5625f * (image->b(row - 1, col) + image->b(row + 1, col)));

                    if(clip) {
                        image->r(row, col) = MIN(image->r(row, col), rtengine::MAXVALF);
                        image->g(row, col) = MIN(image->g(row, col), rtengine::MAXVALF);
                        image->b(row, col) = MIN(image->b(row, col), rtengine::MAXVALF);
                    }
                }
            }

            break;

        case TR_R90: // rotate right
            if( i == 0) {
                for (int j = 0; j < imwidth; j++) {
                    image->r(j, 2 * imheight - 1) = red[j];
                    image->g(j, 2 * imheight - 1) = green[j];
                    image->b(j, 2 * imheight - 1) = blue[j];
                }
            }

            for (int j = 0; j < imwidth; j++) {
                image->r(j, 2 * (imheight - 1 - i)) = red[j];
                image->g(j, 2 * (imheight - 1 - i)) = green[j];
                image->b(j, 2 * (imheight - 1 - i)) = blue[j];
            }

            if (i == 1 || i == 2) { // linear interpolation
                int col = 2 * imheight - 1 - 2 * i;

                for (int j = 0; j < imwidth; j++) {
                    image->r(j, col) = (red[j] + image->r(j, col + 1)) / 2;
                    image->g(j, col) = (green[j] + image->g(j, col + 1)) / 2;
                    image->b(j, col) = (blue[j] + image->b(j, col + 1)) / 2;

                    if(oddHeight && i == 2) {
                        image->r(j, 2 * imheight) = (red[j] + image->r(j, 2 * imheight - 2)) / 2;
                        image->g(j, 2 * imheight) = (green[j] + image->g(j, 2 * imheight - 2)) / 2;
                        image->b(j, 2 * imheight) = (blue[j] + image->b(j, 2 * imheight - 2)) / 2;
                    }
                }
            } else if (i == imheight - 1) {
                int col = 2 * imheight - 1 - 2 * i;

                for (int j = 0; j < imwidth; j++) {
                    image->r(j, col) = (red[j] + image->r(j, col + 1)) / 2;
                    image->g(j, col) = (green[j] + image->g(j, col + 1)) / 2;
                    image->b(j, col) = (blue[j] + image->b(j, col + 1)) / 2;
                }

                col = 2 * imheight - 1 - 2 * i + 2;

                for (int j = 0; j < imwidth; j++) {
                    image->r(j, col) = (red[j] + image->r(j, col + 1)) / 2;
                    image->g(j, col) = (green[j] + image->g(j, col + 1)) / 2;
                    image->b(j, col) = (blue[j] + image->b(j, col + 1)) / 2;
                }
            } else if (i > 2 && i < imheight - 1) { // vertical bicubic interpolation
                int col = 2 * imheight - 1 - 2 * i + 2;

                for (int j = 0; j < imwidth; j++) {
                    image->r(j, col) = MAX(0.f, -0.0625f * (red[j] + image->r(j, col + 3)) + 0.5625f * (image->r(j, col - 1) + image->r(j, col + 1)));
                    image->g(j, col) = MAX(0.f, -0.0625f * (green[j] + image->g(j, col + 3)) + 0.5625f * (image->g(j, col - 1) + image->g(j, col + 1)));
                    image->b(j, col) = MAX(0.f, -0.0625f * (blue[j] + image->b(j, col + 3)) + 0.5625f * (image->b(j, col - 1) + image->b(j, col + 1)));

                    if(clip) {
                        image->r(j, col) = MIN(image->r(j, col), rtengine::MAXVALF);
                        image->g(j, col) = MIN(image->g(j, col), rtengine::MAXVALF);
                        image->b(j, col) = MIN(image->b(j, col), rtengine::MAXVALF);
                    }
                }
            }

            break;

        case TR_R270: // rotate left
            if (i == 0) {
                for (int j = imwidth - 1, row = 0; j >= 0; j--, row++) {
                    image->r(row, 2 * i) = red[j];
                    image->g(row, 2 * i) = green[j];
                    image->b(row, 2 * i) = blue[j];
                }
            } else if (i == 1 || i == 2) { // linear interpolation
                for (int j = imwidth - 1, row = 0; j >= 0; j--, row++) {
                    image->r(row, 2 * i) = red[j];
                    image->g(row, 2 * i) = green[j];
                    image->b(row, 2 * i) = blue[j];
                    image->r(row, 2 * i - 1) = (red[j] + image->r(row, 2 * i - 2)) * 0.5f;
                    image->g(row, 2 * i - 1) = (green[j] + image->g(row, 2 * i - 2)) * 0.5f;
                    image->b(row, 2 * i - 1) = (blue[j] + image->b(row, 2 * i - 2)) * 0.5f;
                }
            } else if (i > 0 && i < imheight) { // vertical bicubic interpolation
                for (int j = imwidth - 1, row = 0; j >= 0; j--, row++) {
                    image->r(row, 2 * i - 3) = MAX(0.f, -0.0625f * (red[j] + image->r(row, 2 * i - 6)) + 0.5625f * (image->r(row, 2 * i - 2) + image->r(row, 2 * i - 4)));
                    image->g(row, 2 * i - 3) = MAX(0.f, -0.0625f * (green[j] + image->g(row, 2 * i - 6)) + 0.5625f * (image->g(row, 2 * i - 2) + image->g(row, 2 * i - 4)));
                    image->b(row, 2 * i - 3) = MAX(0.f, -0.0625f * (blue[j] + image->b(row, 2 * i - 6)) + 0.5625f * (image->b(row, 2 * i - 2) + image->b(row, 2 * i - 4)));

                    if(clip) {
                        image->r(row, 2 * i - 3) = MIN(image->r(row, 2 * i - 3), rtengine::MAXVALF);
                        image->g(row, 2 * i - 3) = MIN(image->g(row, 2 * i - 3), rtengine::MAXVALF);
                        image->b(row, 2 * i - 3) = MIN(image->b(row, 2 * i - 3), rtengine::MAXVALF);
                    }

                    image->r(row, 2 * i) = red[j];
                    image->g(row, 2 * i) = green[j];
                    image->b(row, 2 * i) = blue[j];
                }
            }

            if (i == imheight - 1) {
                for (int j = imwidth - 1, row = 0; j >= 0; j--, row++) {
                    image->r(row, 2 * i - 1) = MAX(0.f, -0.0625f * (red[j] + image->r(row, 2 * i - 4)) + 0.5625f * (image->r(row, 2 * i) + image->r(row, 2 * i - 2)));
                    image->g(row, 2 * i - 1) = MAX(0.f, -0.0625f * (green[j] + image->g(row, 2 * i - 4)) + 0.5625f * (image->g(row, 2 * i) + image->g(row, 2 * i - 2)));
                    image->b(row, 2 * i - 1) = MAX(0.f, -0.0625f * (blue[j] + image->b(row, 2 * i - 4)) + 0.5625f * (image->b(row, 2 * i) + image->b(row, 2 * i - 2)));

                    if(clip) {
                        image->r(j, 2 * i - 1) = MIN(image->r(j, 2 * i - 1), rtengine::MAXVALF);
                        image->g(j, 2 * i - 1) = MIN(image->g(j, 2 * i - 1), rtengine::MAXVALF);
                        image->b(j, 2 * i - 1) = MIN(image->b(j, 2 * i - 1), rtengine::MAXVALF);
                    }

                    image->r(row, 2 * i + 1) = (red[j] + image->r(row, 2 * i - 1)) / 2;
                    image->g(row, 2 * i + 1) = (green[j] + image->g(row, 2 * i - 1)) / 2;
                    image->b(row, 2 * i + 1) = (blue[j] + image->b(row, 2 * i - 1)) / 2;

                    if (oddHeight) {
                        image->r(row, 2 * i + 2) = (red[j] + image->r(row, 2 * i - 2)) / 2;
                        image->g(row, 2 * i + 2) = (green[j] + image->g(row, 2 * i - 2)) / 2;
                        image->b(row, 2 * i + 2) = (blue[j] + image->b(row, 2 * i - 2)) / 2;
                    }
                }
            }

            break;

        case TR_NONE: // no coarse rotation
        default:
            rotateLine (red, image->r, tran, 2 * i, imwidth, imheight);
            rotateLine (green, image->g, tran, 2 * i, imwidth, imheight);
            rotateLine (blue, image->b, tran, 2 * i, imwidth, imheight);

            if (i == 1 || i == 2) { // linear interpolation
                for (int j = 0; j < imwidth; j++) {
                    image->r(2 * i - 1, j) = (red[j] + image->r(2 * i - 2, j)) / 2;
                    image->g(2 * i - 1, j) = (green[j] + image->g(2 * i - 2, j)) / 2;
                    image->b(2 * i - 1, j) = (blue[j] + image->b(2 * i - 2, j)) / 2;
                }
            } else if (i > 2 && i < imheight) { // vertical bicubic interpolation
                for (int j = 0; j < imwidth; j++) {
                    image->r(2 * i - 3, j) = MAX(0.f, -0.0625f * (red[j] + image->r(2 * i - 6, j)) + 0.5625f * (image->r(2 * i - 2, j) + image->r(2 * i - 4, j)));
                    image->g(2 * i - 3, j) = MAX(0.f, -0.0625f * (green[j] + image->g(2 * i - 6, j)) + 0.5625f * (image->g(2 * i - 2, j) + image->g(2 * i - 4, j)));
                    image->b(2 * i - 3, j) = MAX(0.f, -0.0625f * (blue[j] + image->b(2 * i - 6, j)) + 0.5625f * (image->b(2 * i - 2, j) + image->b(2 * i - 4, j)));

                    if(clip) {
                        image->r(2 * i - 3, j) = MIN(image->r(2 * i - 3, j), rtengine::MAXVALF);
                        image->g(2 * i - 3, j) = MIN(image->g(2 * i - 3, j), rtengine::MAXVALF);
                        image->b(2 * i - 3, j) = MIN(image->b(2 * i - 3, j), rtengine::MAXVALF);
                    }
                }
            }

            if (i == imheight - 1) {
                for (int j = 0; j < imwidth; j++) {
                    image->r(2 * i - 1, j) = MAX(0.f, -0.0625f * (red[j] + image->r(2 * i - 4, j)) + 0.5625f * (image->r(2 * i, j) + image->r(2 * i - 2, j)));
                    image->g(2 * i - 1, j) = MAX(0.f, -0.0625f * (green[j] + image->g(2 * i - 4, j)) + 0.5625f * (image->g(2 * i, j) + image->g(2 * i - 2, j)));
                    image->b(2 * i - 1, j) = MAX(0.f, -0.0625f * (blue[j] + image->b(2 * i - 4, j)) + 0.5625f * (image->b(2 * i, j) + image->b(2 * i - 2, j)));

                    if(clip) {
                        image->r(2 * i - 1, j) = MIN(image->r(2 * i - 1, j), rtengine::MAXVALF);
                        image->g(2 * i - 1, j) = MIN(image->g(2 * i - 1, j), rtengine::MAXVALF);
                        image->b(2 * i - 1, j) = MIN(image->b(2 * i - 1, j), rtengine::MAXVALF);
                    }

                    image->r(2 * i + 1, j) = (red[j] + image->r(2 * i - 1, j)) / 2;
                    image->g(2 * i + 1, j) = (green[j] + image->g(2 * i - 1, j)) / 2;
                    image->b(2 * i + 1, j) = (blue[j] + image->b(2 * i - 1, j)) / 2;

                    if (oddHeight) {
                        image->r(2 * i + 2, j) = (red[j] + image->r(2 * i - 2, j)) / 2;
                        image->g(2 * i + 2, j) = (green[j] + image->g(2 * i - 2, j)) / 2;
                        image->b(2 * i + 2, j) = (blue[j] + image->b(2 * i - 2, j)) / 2;
                    }
                }
            }
    }
}

}


namespace rtengine
{

extern const Settings* settings;
#undef ABS
#undef DIST

#define ABS(a) ((a)<0?-(a):(a))
#define DIST(a,b) (ABS(a-b))

RawImageSource::RawImageSource ()
    : ImageSource()
    , W(0), H(0)
    , plistener(nullptr)
    , scale_mul{}
    , c_black{}
    , c_white{}
    , cblacksom{}
    , ref_pre_mul{}
    , refwb_red(0.0)
    , refwb_green(0.0)
    , refwb_blue(0.0)
    , rgb_cam{}
    , cam_rgb{}
    , xyz_cam{}
    , cam_xyz{}
    , fuji(false)
    , d1x(false)
    , border(4)
    , chmax{}
    , hlmax{}
    , clmax{}
    , initialGain(0.0)
    , camInitialGain(0.0)
    , defGain(0.0)
    , ri(nullptr)
    , lc00(0.0)
    , lc01(0.0)
    , lc02(0.0)
    , lc10(0.0)
    , lc11(0.0)
    , lc12(0.0)
    , lc20(0.0)
    , lc21(0.0)
    , lc22(0.0)
    , cache(nullptr)
    , threshold(0)
    , rawData(0, 0)
    , green(0, 0)
    , red(0, 0)
    , blue(0, 0)
    , rawDirty(true)
{
    camProfile = nullptr;
    embProfile = nullptr;
    rgbSourceModified = false;
    for(int i = 0; i < 4; ++i) {
        psRedBrightness[i] = psGreenBrightness[i] = psBlueBrightness[i] = 1.f;
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RawImageSource::~RawImageSource ()
{

    delete idata;

    for(size_t i = 0; i < numFrames; ++i) {
        delete riFrames[i];
    }

    for(size_t i = 0; i < numFrames - 1; ++i) {
        delete rawDataBuffer[i];
    }

    flushRGB();
    flushRawData();

    if( cache ) {
        delete [] cache;
    }

    if (camProfile) {
        cmsCloseProfile (camProfile);
    }

    if (embProfile) {
        cmsCloseProfile (embProfile);
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::transformRect (const PreviewProps &pp, int tran, int &ssx1, int &ssy1, int &width, int &height, int &fw)
{
    int pp_x = pp.getX() + border;
    int pp_y = pp.getY() + border;
    int pp_width = pp.getWidth();
    int pp_height = pp.getHeight();

    if (d1x) {
        if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
            pp_x /= 2;
            pp_width = pp_width / 2 + 1;
        } else {
            pp_y /= 2;
            pp_height = pp_height / 2 + 1;
        }
    }

    int w = W, h = H;

    if (fuji) {
        w = ri->get_FujiWidth() * 2 + 1;
        h = (H - ri->get_FujiWidth()) * 2 + 1;
    }

    int sw = w, sh = h;

    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = h;
        sh = w;
    }

    if( pp_width > sw - 2 * border) {
        pp_width = sw - 2 * border;
    }

    if( pp_height > sh - 2 * border) {
        pp_height = sh - 2 * border;
    }

    int ppx = pp_x, ppy = pp_y;

    if (tran & TR_HFLIP) {
        ppx = max(sw - pp_x - pp_width, 0);
    }

    if (tran & TR_VFLIP) {
        ppy = max(sh - pp_y - pp_height, 0);
    }

    int sx1 = ppx;        // assuming it's >=0
    int sy1 = ppy;        // assuming it's >=0
    int sx2 = min(ppx + pp_width, w - 1);
    int sy2 = min(ppy + pp_height, h - 1);

    if ((tran & TR_ROT) == TR_R180) {
        sx1 = max(w - ppx - pp_width, 0);
        sy1 = max(h - ppy - pp_height, 0);
        sx2 = min(sx1 + pp_width, w - 1);
        sy2 = min(sy1 + pp_height, h - 1);
    } else if ((tran & TR_ROT) == TR_R90) {
        sx1 = ppy;
        sy1 = max(h - ppx - pp_width, 0);
        sx2 = min(sx1 + pp_height, w - 1);
        sy2 = min(sy1 + pp_width, h - 1);
    } else if ((tran & TR_ROT) == TR_R270) {
        sx1 = max(w - ppy - pp_height, 0);
        sy1 = ppx;
        sx2 = min(sx1 + pp_height, w - 1);
        sy2 = min(sy1 + pp_width, h - 1);
    }

    if (fuji) {
        // atszamoljuk a koordinatakat fuji-ra:
        // recalculate the coordinates fuji-ra:
        ssx1 = (sx1 + sy1) / 2;
        ssy1 = (sy1 - sx2 ) / 2 + ri->get_FujiWidth();
        int ssx2 = (sx2 + sy2) / 2 + 1;
        int ssy2 = (sy2 - sx1) / 2 + ri->get_FujiWidth();
        fw   = (sx2 - sx1) / 2 / pp.getSkip();
        width  = (ssx2 - ssx1) / pp.getSkip() + ((ssx2 - ssx1) % pp.getSkip() > 0);
        height = (ssy2 - ssy1) / pp.getSkip() + ((ssy2 - ssy1) % pp.getSkip() > 0);
    } else {
        ssx1 = sx1;
        ssy1 = sy1;
        width  = (sx2 - sx1) / pp.getSkip() + ((sx2 - sx1) % pp.getSkip() > 0);
        height = (sy2 - sy1) / pp.getSkip() + ((sy2 - sy1) % pp.getSkip() > 0);
    }
}

float calculate_scale_mul(float scale_mul[4], const float pre_mul_[4], const float c_white[4], const float c_black[4], bool isMono, int colors)
{
    if (isMono || colors == 1) {
        for (int c = 0; c < 4; c++) {
            scale_mul[c] = 65535.0 / (c_white[c] - c_black[c]);
        }
    } else {
        float pre_mul[4];

        for (int c = 0; c < 4; c++) {
            pre_mul[c] = pre_mul_[c];
        }

        if (pre_mul[3] == 0) {
            pre_mul[3] = pre_mul[1]; // G2 == G1
        }

        float maxpremul = max(pre_mul[0], pre_mul[1], pre_mul[2], pre_mul[3]);

        for (int c = 0; c < 4; c++) {
            scale_mul[c] = (pre_mul[c] / maxpremul) * 65535.0 / (c_white[c] - c_black[c]);
        }
    }

    float gain = max(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]) / min(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);
    return gain;
}

void RawImageSource::getImage (const ColorTemp &ctemp, int tran, Imagefloat* image, const PreviewProps &pp, const ToneCurveParams &hrp, const RAWParams &raw )
{
    MyMutex::MyLock lock(getImageMutex);

    tran = defTransform (tran);

    // compute channel multipliers
    double r, g, b;
    float rm, gm, bm;

    if (ctemp.getTemp() < 0) {
        // no white balance, ie revert the pre-process white balance to restore original unbalanced raw camera color
        rm = ri->get_pre_mul(0);
        gm = ri->get_pre_mul(1);
        bm = ri->get_pre_mul(2);
    } else {
        ctemp.getMultipliers (r, g, b);
        rm = imatrices.cam_rgb[0][0] * r + imatrices.cam_rgb[0][1] * g + imatrices.cam_rgb[0][2] * b;
        gm = imatrices.cam_rgb[1][0] * r + imatrices.cam_rgb[1][1] * g + imatrices.cam_rgb[1][2] * b;
        bm = imatrices.cam_rgb[2][0] * r + imatrices.cam_rgb[2][1] * g + imatrices.cam_rgb[2][2] * b;
    }

    if (true) {
        // adjust gain so the maximum raw value of the least scaled channel just hits max
        const float new_pre_mul[4] = { ri->get_pre_mul(0) / rm, ri->get_pre_mul(1) / gm, ri->get_pre_mul(2) / bm, ri->get_pre_mul(3) / gm };
        float new_scale_mul[4];

        bool isMono = (ri->getSensorType() == ST_FUJI_XTRANS && raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::mono])
                      || (ri->getSensorType() == ST_BAYER && raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::mono]);
        float gain = calculate_scale_mul(new_scale_mul, new_pre_mul, c_white, cblacksom, isMono, ri->get_colors());
        rm = new_scale_mul[0] / scale_mul[0] * gain;
        gm = new_scale_mul[1] / scale_mul[1] * gain;
        bm = new_scale_mul[2] / scale_mul[2] * gain;
        //fprintf(stderr, "camera gain: %f, current wb gain: %f, diff in stops %f\n", camInitialGain, gain, log2(camInitialGain) - log2(gain));
    } else {
        // old scaling: used a fixed reference gain based on camera (as-shot) white balance

        // how much we need to scale each channel to get our new white balance
        rm = refwb_red / rm;
        gm = refwb_green / gm;
        bm = refwb_blue / bm;
        // normalize so larger multiplier becomes 1.0
        float minval = min(rm, gm, bm);
        rm /= minval;
        gm /= minval;
        bm /= minval;
        // multiply with reference gain, ie as-shot WB
        rm *= camInitialGain;
        gm *= camInitialGain;
        bm *= camInitialGain;
    }

    defGain = 0.0;
    // compute image area to render in order to provide the requested part of the image
    int sx1, sy1, imwidth, imheight, fw, d1xHeightOdd = 0;
    transformRect (pp, tran, sx1, sy1, imwidth, imheight, fw);

    // check possible overflows
    int maximwidth, maximheight;

    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        maximwidth = image->getHeight();
        maximheight = image->getWidth();
    } else {
        maximwidth = image->getWidth();
        maximheight = image->getHeight();
    }

    if (d1x) {
        // D1X has only half of the required rows
        // we interpolate the missing ones later to get correct aspect ratio
        // if the height is odd we also have to add an additional row to avoid a black line
        d1xHeightOdd = maximheight & 1;
        maximheight /= 2;
        imheight = maximheight;
    }

    // correct if overflow (very rare), but not fuji because it is corrected in transline
    if (!fuji && imwidth > maximwidth) {
        imwidth = maximwidth;
    }

    if (!fuji && imheight > maximheight) {
        imheight = maximheight;
    }

    if (fuji) { // zero image to avoid access to uninitialized values in further processing because fuji super-ccd processing is not clean...
        for (int i = 0; i < image->getHeight(); ++i) {
            for (int j = 0; j < image->getWidth(); ++j) {
                image->r(i, j) = image->g(i, j) = image->b(i, j) = 0;
            }
        }
    }

    int maxx = this->W, maxy = this->H, skip = pp.getSkip();

    // raw clip levels after white balance
    hlmax[0] = clmax[0] * rm;
    hlmax[1] = clmax[1] * gm;
    hlmax[2] = clmax[2] * bm;

    const bool doClip = (chmax[0] >= clmax[0] || chmax[1] >= clmax[1] || chmax[2] >= clmax[2]) && !hrp.hrenabled;

    float area = skip * skip;
    rm /= area;
    gm /= area;
    bm /= area;
    bool doHr = (hrp.hrenabled && hrp.method != "Color");
#ifdef _OPENMP
    #pragma omp parallel if(!d1x)       // omp disabled for D1x to avoid race conditions (see Issue 1088 http://code.google.com/p/rawtherapee/issues/detail?id=1088)
    {
#endif
        // render the requested image part
        float line_red[imwidth] ALIGNED16;
        float line_grn[imwidth] ALIGNED16;
        float line_blue[imwidth] ALIGNED16;

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int ix = 0; ix < imheight; ix++) {
            int i = sy1 + skip * ix;

            if (i >= maxy - skip) {
                i = maxy - skip - 1;    // avoid trouble
            }

            if (ri->getSensorType() == ST_BAYER || ri->getSensorType() == ST_FUJI_XTRANS || ri->get_colors() == 1) {
                for (int j = 0, jx = sx1; j < imwidth; j++, jx += skip) {
                    jx = std::min(jx, maxx - skip - 1); // avoid trouble

                    float rtot = 0.f, gtot = 0.f, btot = 0.f;

                    for (int m = 0; m < skip; m++)
                        for (int n = 0; n < skip; n++) {
                            rtot += red[i + m][jx + n];
                            gtot += green[i + m][jx + n];
                            btot += blue[i + m][jx + n];
                        }

                    rtot *= rm;
                    gtot *= gm;
                    btot *= bm;

                    if (doClip) {
                        // note: as hlmax[] can be larger than CLIP and we can later apply negative
                        // exposure this means that we can clip away local highlights which actually
                        // are not clipped. We have to do that though as we only check pixel by pixel
                        // and don't know if this will transition into a clipped area, if so we need
                        // to clip also surrounding to make a good colour transition
                        rtot = CLIP(rtot);
                        gtot = CLIP(gtot);
                        btot = CLIP(btot);
                    }

                    line_red[j] = rtot;
                    line_grn[j] = gtot;
                    line_blue[j] = btot;
                }
            } else {
                for (int j = 0, jx = sx1; j < imwidth; j++, jx += skip) {
                    if (jx > maxx - skip) {
                        jx = maxx - skip - 1;
                    }

                    float rtot, gtot, btot;
                    rtot = gtot = btot = 0;

                    for (int m = 0; m < skip; m++)
                        for (int n = 0; n < skip; n++) {
                            rtot += rawData[i + m][(jx + n) * 3 + 0];
                            gtot += rawData[i + m][(jx + n) * 3 + 1];
                            btot += rawData[i + m][(jx + n) * 3 + 2];
                        }

                    rtot *= rm;
                    gtot *= gm;
                    btot *= bm;

                    if (doClip) {
                        rtot = CLIP(rtot);
                        gtot = CLIP(gtot);
                        btot = CLIP(btot);
                    }

                    line_red[j] = rtot;
                    line_grn[j] = gtot;
                    line_blue[j] = btot;

                }
            }

            //process all highlight recovery other than "Color"
            if (doHr) {
                hlRecovery (hrp.method, line_red, line_grn, line_blue, imwidth, hlmax);
            }

            if(d1x) {
                transLineD1x (line_red, line_grn, line_blue, ix, image, tran, imwidth, imheight, d1xHeightOdd, doClip);
            } else if(fuji) {
                transLineFuji (line_red, line_grn, line_blue, ix, image, tran, imheight, fw);
            } else {
                transLineStandard (line_red, line_grn, line_blue, ix, image, tran, imwidth, imheight);
            }

        }

#ifdef _OPENMP
    }
#endif

    if (fuji) {
        int a = ((tran & TR_ROT) == TR_R90 && image->getWidth() % 2 == 0) || ((tran & TR_ROT) == TR_R180 && image->getHeight() % 2 + image->getWidth() % 2 == 1) || ((tran & TR_ROT) == TR_R270 && image->getHeight() % 2 == 0);

        // first row
        for (int j = 1 + a; j < image->getWidth() - 1; j += 2) {
            image->r(0, j) = (image->r(1, j) + image->r(0, j + 1) + image->r(0, j - 1)) / 3;
            image->g(0, j) = (image->g(1, j) + image->g(0, j + 1) + image->g(0, j - 1)) / 3;
            image->b(0, j) = (image->b(1, j) + image->b(0, j + 1) + image->b(0, j - 1)) / 3;
        }

        // other rows
        for (int i = 1; i < image->getHeight() - 1; i++) {
            for (int j = 2 - (a + i + 1) % 2; j < image->getWidth() - 1; j += 2) {
                // edge-adaptive interpolation
                double dh = (ABS(image->r(i, j + 1) - image->r(i, j - 1)) + ABS(image->g(i, j + 1) - image->g(i, j - 1)) + ABS(image->b(i, j + 1) - image->b(i, j - 1))) / 1.0;
                double dv = (ABS(image->r(i + 1, j) - image->r(i - 1, j)) + ABS(image->g(i + 1, j) - image->g(i - 1, j)) + ABS(image->b(i + 1, j) - image->b(i - 1, j))) / 1.0;
                double eh = 1.0 / (1.0 + dh);
                double ev = 1.0 / (1.0 + dv);
                image->r(i, j) = (eh * (image->r(i, j + 1) + image->r(i, j - 1)) + ev * (image->r(i + 1, j) + image->r(i - 1, j))) / (2.0 * (eh + ev));
                image->g(i, j) = (eh * (image->g(i, j + 1) + image->g(i, j - 1)) + ev * (image->g(i + 1, j) + image->g(i - 1, j))) / (2.0 * (eh + ev));
                image->b(i, j) = (eh * (image->b(i, j + 1) + image->b(i, j - 1)) + ev * (image->b(i + 1, j) + image->b(i - 1, j))) / (2.0 * (eh + ev));
            }

            // first pixel
            if (2 - (a + i + 1) % 2 == 2) {
                image->r(i, 0) = (image->r(i + 1, 0) + image->r(i - 1, 0) + image->r(i, 1)) / 3;
                image->g(i, 0) = (image->g(i + 1, 0) + image->g(i - 1, 0) + image->g(i, 1)) / 3;
                image->b(i, 0) = (image->b(i + 1, 0) + image->b(i - 1, 0) + image->b(i, 1)) / 3;
            }

            // last pixel
            if (2 - (a + i + image->getWidth()) % 2 == 2) {
                image->r(i, image->getWidth() - 1) = (image->r(i + 1, image->getWidth() - 1) + image->r(i - 1, image->getWidth() - 1) + image->r(i, image->getWidth() - 2)) / 3;
                image->g(i, image->getWidth() - 1) = (image->g(i + 1, image->getWidth() - 1) + image->g(i - 1, image->getWidth() - 1) + image->g(i, image->getWidth() - 2)) / 3;
                image->b(i, image->getWidth() - 1) = (image->b(i + 1, image->getWidth() - 1) + image->b(i - 1, image->getWidth() - 1) + image->b(i, image->getWidth() - 2)) / 3;
            }
        }

        // last row
        int b = (a == 1 && image->getHeight() % 2) || (a == 0 && image->getHeight() % 2 == 0);

        for (int j = 1 + b; j < image->getWidth() - 1; j += 2) {
            image->r(image->getHeight() - 1, j) = (image->r(image->getHeight() - 2, j) + image->r(image->getHeight() - 1, j + 1) + image->r(image->getHeight() - 1, j - 1)) / 3;
            image->g(image->getHeight() - 1, j) = (image->g(image->getHeight() - 2, j) + image->g(image->getHeight() - 1, j + 1) + image->g(image->getHeight() - 1, j - 1)) / 3;
            image->b(image->getHeight() - 1, j) = (image->b(image->getHeight() - 2, j) + image->b(image->getHeight() - 1, j + 1) + image->b(image->getHeight() - 1, j - 1)) / 3;
        }
    }

    // Flip if needed
    if (tran & TR_HFLIP) {
        hflip (image);
    }

    if (tran & TR_VFLIP) {
        vflip (image);
    }

    // Colour correction (only when running on full resolution)
    if(pp.getSkip() == 1) {
        switch(ri->getSensorType()) {
            case ST_BAYER:
                processFalseColorCorrection (image, raw.bayersensor.ccSteps);
                break;

            case ST_FUJI_XTRANS:
                processFalseColorCorrection (image, raw.xtranssensor.ccSteps);
                break;

            case ST_FOVEON:
            case ST_NONE:
                break;
        }
    }
}

DCPProfile *RawImageSource::getDCP(const ColorManagementParams &cmp, DCPProfile::ApplyState &as)
{
    DCPProfile *dcpProf = nullptr;
    cmsHPROFILE dummy;
    findInputProfile(cmp.input, nullptr, (static_cast<const FramesData*>(getMetaData()))->getCamera(), &dcpProf, dummy);

    if (dcpProf == nullptr) {
        if (settings->verbose) {
            printf("Can't load DCP profile '%s'!\n", cmp.input.c_str());
        }
        return nullptr;
    }

    dcpProf->setStep2ApplyState(cmp.working, cmp.toneCurve, cmp.applyLookTable, cmp.applyBaselineExposureOffset, as);
    return dcpProf;
}

void RawImageSource::convertColorSpace(Imagefloat* image, const ColorManagementParams &cmp, const ColorTemp &wb)
{
    double pre_mul[3] = { ri->get_pre_mul(0), ri->get_pre_mul(1), ri->get_pre_mul(2) };
    colorSpaceConversion (image, cmp, wb, pre_mul, embProfile, camProfile, imatrices.xyz_cam, (static_cast<const FramesData*>(getMetaData()))->getCamera());
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/* interpolateBadPixelsBayer: correct raw pixels looking at the bitmap
 * takes into consideration if there are multiple bad pixels in the neighbourhood
 */
int RawImageSource::interpolateBadPixelsBayer( PixelsMap &bitmapBads, array2D<float> &rawData )
{
    static const float eps = 1.f;
    int counter = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:counter) schedule(dynamic,16)
#endif

    for( int row = 2; row < H - 2; row++ ) {
        for(int col = 2; col < W - 2; col++ ) {
            int sk = bitmapBads.skipIfZero(col, row); //optimization for a stripe all zero

            if( sk ) {
                col += sk - 1; //-1 is because of col++ in cycle
                continue;
            }

            if(!bitmapBads.get(col, row)) {
                continue;
            }

            float wtdsum = 0.f, norm = 0.f;

            // diagonal interpolation
            if(FC(row, col) == 1) {
                // green channel. We can use closer pixels than for red or blue channel. Distance to centre pixel is sqrt(2) => weighting is 0.70710678
                // For green channel following pixels will be used for interpolation. Pixel to be interpolated is in centre.
                // 1 means that pixel is used in this step, if itself and his counterpart are not marked bad
                // 0 0 0 0 0
                // 0 1 0 1 0
                // 0 0 0 0 0
                // 0 1 0 1 0
                // 0 0 0 0 0
                for( int dx = -1; dx <= 1; dx += 2) {
                    if( bitmapBads.get(col + dx, row - 1) || bitmapBads.get(col - dx, row + 1)) {
                        continue;
                    }

                    float dirwt = 0.70710678f / ( fabsf( rawData[row - 1][col + dx] - rawData[row + 1][col - dx]) + eps);
                    wtdsum += dirwt * (rawData[row - 1][col + dx] + rawData[row + 1][col - dx]);
                    norm += dirwt;
                }
            } else {
                // red and blue channel. Distance to centre pixel is sqrt(8) => weighting is 0.35355339
                // For red and blue channel following pixels will be used for interpolation. Pixel to be interpolated is in centre.
                // 1 means that pixel is used in this step, if itself and his counterpart are not marked bad
                // 1 0 0 0 1
                // 0 0 0 0 0
                // 0 0 0 0 0
                // 0 0 0 0 0
                // 1 0 0 0 1
                for( int dx = -2; dx <= 2; dx += 4) {
                    if( bitmapBads.get(col + dx, row - 2) || bitmapBads.get(col - dx, row + 2)) {
                        continue;
                    }

                    float dirwt = 0.35355339f / ( fabsf( rawData[row - 2][col + dx] - rawData[row + 2][col - dx]) + eps);
                    wtdsum += dirwt * (rawData[row - 2][col + dx] + rawData[row + 2][col - dx]);
                    norm += dirwt;
                }
            }

            // channel independent. Distance to centre pixel is 2 => weighting is 0.5
            // Additionally for all channel following pixels will be used for interpolation. Pixel to be interpolated is in centre.
            // 1 means that pixel is used in this step, if itself and his counterpart are not marked bad
            // 0 0 1 0 0
            // 0 0 0 0 0
            // 1 0 0 0 1
            // 0 0 0 0 0
            // 0 0 1 0 0

            // horizontal interpolation
            if(!(bitmapBads.get(col - 2, row) || bitmapBads.get(col + 2, row))) {
                float dirwt = 0.5f / ( fabsf( rawData[row][col - 2] - rawData[row][col + 2]) + eps);
                wtdsum += dirwt * (rawData[row][col - 2] + rawData[row][col + 2]);
                norm += dirwt;
            }

            // vertical interpolation
            if(!(bitmapBads.get(col, row - 2) || bitmapBads.get(col, row + 2))) {
                float dirwt = 0.5f / ( fabsf( rawData[row - 2][col] - rawData[row + 2][col]) + eps);
                wtdsum += dirwt * (rawData[row - 2][col] + rawData[row + 2][col]);
                norm += dirwt;
            }

            if (LIKELY(norm > 0.f)) { // This means, we found at least one pair of valid pixels in the steps above, likelihood of this case is about 99.999%
                rawData[row][col] = wtdsum / (2.f * norm); //gradient weighted average, Factor of 2.f is an optimization to avoid multiplications in former steps
                counter++;
            } else { //backup plan -- simple average. Same method for all channels. We could improve this, but it's really unlikely that this case happens
                int tot = 0;
                float sum = 0;

                for( int dy = -2; dy <= 2; dy += 2) {
                    for( int dx = -2; dx <= 2; dx += 2) {
                        if(bitmapBads.get(col + dx, row + dy)) {
                            continue;
                        }

                        sum += rawData[row + dy][col + dx];
                        tot++;
                    }
                }

                if (tot > 0) {
                    rawData[row][col] = sum / tot;
                    counter ++;
                }
            }
        }
    }

    return counter; // Number of interpolated pixels.
}

/* interpolateBadPixels3Colours: correct raw pixels looking at the bitmap
 * takes into consideration if there are multiple bad pixels in the neighbourhood
 */
int RawImageSource::interpolateBadPixelsNColours( PixelsMap &bitmapBads, const int colours )
{
    static const float eps = 1.f;
    int counter = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:counter) schedule(dynamic,16)
#endif

    for( int row = 2; row < H - 2; row++ ) {
        for(int col = 2; col < W - 2; col++ ) {
            int sk = bitmapBads.skipIfZero(col, row); //optimization for a stripe all zero

            if( sk ) {
                col += sk - 1; //-1 is because of col++ in cycle
                continue;
            }

            if(!bitmapBads.get(col, row)) {
                continue;
            }

            float wtdsum[colours], norm[colours];

            for (int i = 0; i < colours; ++i) {
                wtdsum[i] = norm[i] = 0.f;
            }

            // diagonal interpolation
            for( int dx = -1; dx <= 1; dx += 2) {
                if( bitmapBads.get(col + dx, row - 1) || bitmapBads.get(col - dx, row + 1)) {
                    continue;
                }

                for(int c = 0; c < colours; c++) {
                    float dirwt = 0.70710678f / ( fabsf( rawData[row - 1][(col + dx) * colours + c] - rawData[row + 1][(col - dx) * colours + c]) + eps);
                    wtdsum[c] += dirwt * (rawData[row - 1][(col + dx) * colours + c] + rawData[row + 1][(col - dx) * colours + c]);
                    norm[c] += dirwt;
                }
            }

            // horizontal interpolation
            if(!(bitmapBads.get(col - 1, row) || bitmapBads.get(col + 1, row))) {
                for(int c = 0; c < colours; c++) {
                    float dirwt = 1.f / ( fabsf( rawData[row][(col - 1) * colours + c] - rawData[row][(col + 1) * colours + c]) + eps);
                    wtdsum[c] += dirwt * (rawData[row][(col - 1) * colours + c] + rawData[row][(col + 1) * colours + c]);
                    norm[c] += dirwt;
                }
            }

            // vertical interpolation
            if(!(bitmapBads.get(col, row - 1) || bitmapBads.get(col, row + 1))) {
                for(int c = 0; c < colours; c++) {
                    float dirwt = 1.f / ( fabsf( rawData[row - 1][col * colours + c] - rawData[row + 1][col * colours + c]) + eps);
                    wtdsum[c] += dirwt * (rawData[row - 1][col * colours + c] + rawData[row + 1][col * colours + c]);
                    norm[c] += dirwt;
                }
            }

            if (LIKELY(norm[0] > 0.f)) { // This means, we found at least one pair of valid pixels in the steps above, likelihood of this case is about 99.999%
                for(int c = 0; c < colours; c++) {
                    rawData[row][col * colours + c] = wtdsum[c] / (2.f * norm[c]); //gradient weighted average, Factor of 2.f is an optimization to avoid multiplications in former steps
                }

                counter++;
            } else { //backup plan -- simple average. Same method for all channels. We could improve this, but it's really unlikely that this case happens
                int tot = 0;
                float sum[colours];

                for (int i = 0; i < colours; ++i) {
                    sum[i] = 0.f;
                }

                for( int dy = -2; dy <= 2; dy += 2) {
                    for( int dx = -2; dx <= 2; dx += 2) {
                        if(bitmapBads.get(col + dx, row + dy)) {
                            continue;
                        }

                        for(int c = 0; c < colours; c++) {
                            sum[c] += rawData[row + dy][(col + dx) * colours + c];
                        }

                        tot++;
                    }
                }

                if (tot > 0) {
                    for(int c = 0; c < colours; c++) {
                        rawData[row][col * colours + c] = sum[c] / tot;
                    }

                    counter ++;
                }
            }
        }
    }

    return counter; // Number of interpolated pixels.
}
/* interpolateBadPixelsXtrans: correct raw pixels looking at the bitmap
 * takes into consideration if there are multiple bad pixels in the neighbourhood
 */
int RawImageSource::interpolateBadPixelsXtrans( PixelsMap &bitmapBads )
{
    static const float eps = 1.f;
    int counter = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:counter) schedule(dynamic,16)
#endif

    for( int row = 2; row < H - 2; row++ ) {
        for(int col = 2; col < W - 2; col++ ) {
            int skip = bitmapBads.skipIfZero(col, row); //optimization for a stripe all zero

            if( skip ) {
                col += skip - 1; //-1 is because of col++ in cycle
                continue;
            }

            if(!bitmapBads.get(col, row)) {
                continue;
            }

            float wtdsum = 0.f, norm = 0.f;
            unsigned int pixelColor = ri->XTRANSFC(row, col);

            if(pixelColor == 1) {
                // green channel. A green pixel can either be a solitary green pixel or a member of a 2x2 square of green pixels
                if(ri->XTRANSFC(row, col - 1) == ri->XTRANSFC(row, col + 1)) {
                    // If left and right neighbour have same colour, then this is a solitary green pixel
                    // For these the following pixels will be used for interpolation. Pixel to be interpolated is in centre and marked with a P.
                    // Pairs of pixels used in this step are numbered. A pair will be used if none of the pixels of the pair is marked bad
                    // 0 means, the pixel has a different colour and will not be used
                    // 0 1 0 2 0
                    // 3 5 0 6 4
                    // 0 0 P 0 0
                    // 4 6 0 5 3
                    // 0 2 0 1 0
                    for( int dx = -1; dx <= 1; dx += 2) { // pixels marked 5 or 6 in above example. Distance to P is sqrt(2) => weighting is 0.70710678f
                        if( bitmapBads.get(col + dx, row - 1) || bitmapBads.get(col - dx, row + 1)) {
                            continue;
                        }

                        float dirwt = 0.70710678f / ( fabsf( rawData[row - 1][col + dx] - rawData[row + 1][col - dx]) + eps);
                        wtdsum += dirwt * (rawData[row - 1][col + dx] + rawData[row + 1][col - dx]);
                        norm += dirwt;
                    }

                    for( int dx = -1; dx <= 1; dx += 2) { // pixels marked 1 or 2 on above example. Distance to P is sqrt(5) => weighting is 0.44721359f
                        if( bitmapBads.get(col + dx, row - 2) || bitmapBads.get(col - dx, row + 2)) {
                            continue;
                        }

                        float dirwt = 0.44721359f / ( fabsf( rawData[row - 2][col + dx] - rawData[row + 2][col - dx]) + eps);
                        wtdsum += dirwt * (rawData[row - 2][col + dx] + rawData[row + 2][col - dx]);
                        norm += dirwt;
                    }

                    for( int dx = -2; dx <= 2; dx += 4) { // pixels marked 3 or 4 on above example. Distance to P is sqrt(5) => weighting is 0.44721359f
                        if( bitmapBads.get(col + dx, row - 1) || bitmapBads.get(col - dx, row + 1)) {
                            continue;
                        }

                        float dirwt = 0.44721359f / ( fabsf( rawData[row - 1][col + dx] - rawData[row + 1][col - dx]) + eps);
                        wtdsum += dirwt * (rawData[row - 1][col + dx] + rawData[row + 1][col - dx]);
                        norm += dirwt;
                    }
                } else {
                    // this is a member of a 2x2 square of green pixels
                    // For these the following pixels will be used for interpolation. Pixel to be interpolated is at position P in the example.
                    // Pairs of pixels used in this step are numbered. A pair will be used if none of the pixels of the pair is marked bad
                    // 0 means, the pixel has a different colour and will not be used
                    // 1 0 0 3
                    // 0 P 2 0
                    // 0 2 1 0
                    // 3 0 0 0

                    // pixels marked 1 in above example. Distance to P is sqrt(2) => weighting is 0.70710678f
                    int offset1 = ri->XTRANSFC(row - 1, col - 1) == ri->XTRANSFC(row + 1, col + 1) ? 1 : -1;

                    if( !(bitmapBads.get(col - offset1, row - 1) || bitmapBads.get(col + offset1, row + 1))) {
                        float dirwt = 0.70710678f / ( fabsf( rawData[row - 1][col - offset1] - rawData[row + 1][col + offset1]) + eps);
                        wtdsum += dirwt * (rawData[row - 1][col - offset1] + rawData[row + 1][col + offset1]);
                        norm += dirwt;
                    }

                    // pixels marked 2 in above example. Distance to P is 1 => weighting is 1.f
                    int offsety = (ri->XTRANSFC(row - 1, col) != 1 ? 1 : -1);
                    int offsetx = offset1 * offsety;

                    if( !(bitmapBads.get(col + offsetx, row) || bitmapBads.get(col, row + offsety))) {
                        float dirwt = 1.f / ( fabsf( rawData[row][col + offsetx] - rawData[row + offsety][col]) + eps);
                        wtdsum += dirwt * (rawData[row][col + offsetx] + rawData[row + offsety][col]);
                        norm += dirwt;
                    }

                    int offsety2 = -offsety;
                    int offsetx2 = -offsetx;
                    offsetx *= 2;
                    offsety *= 2;

                    // pixels marked 3 in above example. Distance to P is sqrt(5) => weighting is 0.44721359f
                    if( !(bitmapBads.get(col + offsetx, row + offsety2) || bitmapBads.get(col + offsetx2, row + offsety))) {
                        float dirwt = 0.44721359f / ( fabsf( rawData[row + offsety2][col + offsetx] - rawData[row + offsety][col + offsetx2]) + eps);
                        wtdsum += dirwt * (rawData[row + offsety2][col + offsetx] + rawData[row + offsety][col + offsetx2]);
                        norm += dirwt;
                    }
                }
            } else {
                // red and blue channel.
                // Each red or blue pixel has exactly one neighbour of same colour in distance 2 and four neighbours of same colour which can be reached by a move of a knight in chess.
                // For the distance 2 pixel (marked with an X) we generate a virtual counterpart (marked with a V)
                // For red and blue channel following pixels will be used for interpolation. Pixel to be interpolated is in centre and marked with a P.
                // Pairs of pixels used in this step are numbered except for distance 2 pixels which are marked X and V. A pair will be used if none of the pixels of the pair is marked bad
                // 0 1 0 0 0    0 0 X 0 0   remaining cases are symmetric
                // 0 0 0 0 2    1 0 0 0 2
                // X 0 P 0 V    0 0 P 0 0
                // 0 0 0 0 1    0 0 0 0 0
                // 0 2 0 0 0    0 2 V 1 0

                // Find two knight moves landing on a pixel of same colour as the pixel to be interpolated.
                // If we look at first and last row of 5x5 square, we will find exactly two knight pixels.
                // Additionally we know that the column of this pixel has 1 or -1 horizontal distance to the centre pixel
                // When we find a knight pixel, we get its counterpart, which has distance (+-3,+-3), where the signs of distance depend on the corner of the found knight pixel.
                // These pixels are marked 1 or 2 in above examples. Distance to P is sqrt(5) => weighting is 0.44721359f
                // The following loop simply scans the four possible places. To keep things simple, it does not stop after finding two knight pixels, because it will not find more than two
                for(int d1 = -2, offsety = 3; d1 <= 2; d1 += 4, offsety -= 6) {
                    for(int d2 = -1, offsetx = 3; d2 < 1; d2 += 2, offsetx -= 6) {
                        if(ri->XTRANSFC(row + d1, col + d2) == pixelColor) {
                            if( !(bitmapBads.get(col + d2, row + d1) || bitmapBads.get(col + d2 + offsetx, row + d1 + offsety))) {
                                float dirwt = 0.44721359f / ( fabsf( rawData[row + d1][col + d2] - rawData[row + d1 + offsety][col + d2 + offsetx]) + eps);
                                wtdsum += dirwt * (rawData[row + d1][col + d2] + rawData[row + d1 + offsety][col + d2 + offsetx]);
                                norm += dirwt;
                            }
                        }
                    }
                }

                // now scan for the pixel of same colour in distance 2 in each direction (marked with an X in above examples).
                bool distance2PixelFound = false;
                int dx, dy;

                // check horizontal
                for(dx = -2, dy = 0; dx <= 2 && !distance2PixelFound; dx += 4)
                    if(ri->XTRANSFC(row, col + dx) == pixelColor) {
                        distance2PixelFound = true;
                    }

                if(!distance2PixelFound)

                    // no distance 2 pixel on horizontal, check vertical
                    for(dx = 0, dy = -2; dy <= 2 && !distance2PixelFound; dy += 4)
                        if(ri->XTRANSFC(row + dy, col) == pixelColor) {
                            distance2PixelFound = true;
                        }

                // calculate the value of its virtual counterpart (marked with a V in above examples)
                float virtualPixel;

                if(dy == 0) {
                    virtualPixel = 0.5f * (rawData[row - 1][col - dx] + rawData[row + 1][col - dx]);
                } else {
                    virtualPixel = 0.5f * (rawData[row - dy][col - 1] + rawData[row - dy][col + 1]);
                }

                // and weight as usual. Distance to P is 2 => weighting is 0.5f
                float dirwt = 0.5f / ( fabsf( virtualPixel - rawData[row + dy][col + dx]) + eps);
                wtdsum += dirwt * (virtualPixel + rawData[row + dy][col + dx]);
                norm += dirwt;
            }

            if (LIKELY(norm > 0.f)) { // This means, we found at least one pair of valid pixels in the steps above, likelihood of this case is about 99.999%
                rawData[row][col] = wtdsum / (2.f * norm); //gradient weighted average, Factor of 2.f is an optimization to avoid multiplications in former steps
                counter++;
            }
        }
    }

    return counter; // Number of interpolated pixels.
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*  Search for hot or dead pixels in the image and update the map
 *  For each pixel compare its value to the average of similar colour surrounding
 *  (Taken from Emil Martinec idea)
 *  (Optimized by Ingo Weyrich 2013 and 2015)
 */
SSEFUNCTION int RawImageSource::findHotDeadPixels( PixelsMap &bpMap, float thresh, bool findHotPixels, bool findDeadPixels )
{
    float varthresh = (20.0 * (thresh / 100.0) + 1.0 ) / 24.f;

    // allocate temporary buffer
    float (*cfablur);
    cfablur = (float (*)) malloc (H * W * sizeof * cfablur);

    // counter for dead or hot pixels
    int counter = 0;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16) nowait
#endif

        for (int i = 2; i < H - 2; i++) {
            for (int j = 2; j < W - 2; j++) {
                const float& temp = median(rawData[i - 2][j - 2], rawData[i - 2][j], rawData[i - 2][j + 2],
                                           rawData[i][j - 2], rawData[i][j], rawData[i][j + 2],
                                           rawData[i + 2][j - 2], rawData[i + 2][j], rawData[i + 2][j + 2]);
                cfablur[i * W + j] = rawData[i][j] - temp;
            }
        }

        // process borders. Former version calculated the median using mirrored border which does not make sense because the original pixel loses weight
        // Setting the difference between pixel and median for border pixels to zero should do the job not worse then former version
#ifdef _OPENMP
        #pragma omp single
#endif
        {
            for(int i = 0; i < 2; i++) {
                for(int j = 0; j < W; j++) {
                    cfablur[i * W + j] = 0.f;
                }
            }

            for(int i = 2; i < H - 2; i++) {
                for(int j = 0; j < 2; j++) {
                    cfablur[i * W + j] = 0.f;
                }

                for(int j = W - 2; j < W; j++) {
                    cfablur[i * W + j] = 0.f;
                }
            }

            for(int i = H - 2; i < H; i++) {
                for(int j = 0; j < W; j++) {
                    cfablur[i * W + j] = 0.f;
                }
            }
        }
#ifdef _OPENMP
        #pragma omp barrier // barrier because of nowait clause above

        #pragma omp for reduction(+:counter) schedule(dynamic,16)
#endif

        //cfa pixel heat/death evaluation
        for (int rr = 2; rr < H - 2; rr++) {
            int rrmWpcc = rr * W + 2;

            for (int cc = 2; cc < W - 2; cc++, rrmWpcc++) {
                //evaluate pixel for heat/death
                float pixdev = cfablur[rrmWpcc];

                if(pixdev == 0.f) {
                    continue;
                }

                if((!findDeadPixels) && pixdev < 0) {
                    continue;
                }

                if((!findHotPixels) && pixdev > 0) {
                    continue;
                }

                pixdev = fabsf(pixdev);
                float hfnbrave = -pixdev;

#ifdef __SSE2__
                // sum up 5*4 = 20 values using SSE
                // 10 fabs function calls and float 10 additions with SSE
                vfloat sum = vabsf(LVFU(cfablur[(rr - 2) * W + cc - 2])) + vabsf(LVFU(cfablur[(rr - 1) * W + cc - 2]));
                sum += vabsf(LVFU(cfablur[(rr) * W + cc - 2]));
                sum += vabsf(LVFU(cfablur[(rr + 1) * W + cc - 2]));
                sum += vabsf(LVFU(cfablur[(rr + 2) * W + cc - 2]));
                // horizontally add the values and add the result to hfnbrave
                hfnbrave += vhadd(sum);

                // add remaining 5 values of last column
                for (int mm = rr - 2; mm <= rr + 2; mm++) {
                    hfnbrave += fabsf(cfablur[mm * W + cc + 2]);
                }

#else

                //  25 fabs function calls and 25 float additions without SSE
                for (int mm = rr - 2; mm <= rr + 2; mm++) {
                    for (int nn = cc - 2; nn <= cc + 2; nn++) {
                        hfnbrave += fabsf(cfablur[mm * W + nn]);
                    }
                }

#endif

                if (pixdev > varthresh * hfnbrave) {
                    // mark the pixel as "bad"
                    bpMap.set(cc, rr);
                    counter++;
                }
            }//end of pixel evaluation
        }
    }//end of parallel processing
    free (cfablur);
    return counter;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::getFullSize (int& w, int& h, int tr)
{

    tr = defTransform (tr);

    if (fuji) {
        w = ri->get_FujiWidth() * 2 + 1;
        h = (H - ri->get_FujiWidth()) * 2 + 1;
    } else if (d1x) {
        w = W;
        h = 2 * H;
    } else {
        w = W;
        h = H;
    }

    if ((tr & TR_ROT) == TR_R90 || (tr & TR_ROT) == TR_R270) {
        int tmp = w;
        w = h;
        h = tmp;
    }

    w -= 2 * border;
    h -= 2 * border;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::getSize (const PreviewProps &pp, int& w, int& h)
{
    w = pp.getWidth() / pp.getSkip() + (pp.getWidth() % pp.getSkip() > 0);
    h = pp.getHeight() / pp.getSkip() + (pp.getHeight() % pp.getSkip() > 0);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hflip (Imagefloat* image)
{
    image->hflip();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::vflip (Imagefloat* image)
{
    image->vflip();
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int RawImageSource::load (const Glib::ustring &fname)
{

    MyTime t1, t2;
    t1.set();
    fileName = fname;

    if (plistener) {
        plistener->setProgressStr ("Decoding...");
        plistener->setProgress (0.0);
    }
    ri = new RawImage(fname);
    int errCode = ri->loadRaw (false, 0, false);

    if (errCode) {
        return errCode;
    }
    numFrames = ri->getFrameCount();

    errCode = 0;

    if(numFrames > 1) {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            int errCodeThr = 0;
#ifdef _OPENMP
            #pragma omp for nowait
#endif
            for(unsigned int i = 0; i < numFrames; ++i) {
                if(i == 0) {
                    riFrames[i] = ri;
                    errCodeThr = riFrames[i]->loadRaw (true, i, true, plistener, 0.8);
                } else {
                    riFrames[i] = new RawImage(fname);
                    errCodeThr = riFrames[i]->loadRaw (true, i);
                }
            }
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                errCode = errCodeThr ? errCodeThr : errCode;
            }
        }
    } else {
        riFrames[0] = ri;
        errCode = riFrames[0]->loadRaw (true, 0, true, plistener, 0.8);
    }

    if(!errCode) {
        for(unsigned int i = 0; i < numFrames; ++i) {
            riFrames[i]->compress_image(i);
        }
    } else {
        return errCode;
    }

    if(numFrames > 1 ) { // this disables multi frame support for Fuji S5 until I found a solution to handle different dimensions
        if(riFrames[0]->get_width() != riFrames[1]->get_width() || riFrames[0]->get_height() != riFrames[1]->get_height()) {
            numFrames = 1;
        }
    }

    if (plistener) {
        plistener->setProgress (0.9);
    }

    /***** Copy once constant data extracted from raw *******/
    W = ri->get_width();
    H = ri->get_height();
    fuji = ri->get_FujiWidth() != 0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            imatrices.rgb_cam[i][j] = ri->get_rgb_cam(i, j);
        }

    // compute inverse of the color transformation matrix
    // first arg is matrix, second arg is inverse
    inverse33 (imatrices.rgb_cam, imatrices.cam_rgb);

    d1x  = ! ri->get_model().compare("D1X");

    if(ri->getSensorType() == ST_FUJI_XTRANS) {
        border = 7;
    } else if(ri->getSensorType() == ST_FOVEON) {
        border = 0;
    }

    if ( ri->get_profile() ) {
        embProfile = cmsOpenProfileFromMem (ri->get_profile(), ri->get_profileLen());
    }

    // create profile
    memset (imatrices.xyz_cam, 0, sizeof(imatrices.xyz_cam));

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) {
                imatrices.xyz_cam[i][j] += xyz_sRGB[i][k] * imatrices.rgb_cam[k][j];
            }

    camProfile = ICCStore::getInstance()->createFromMatrix (imatrices.xyz_cam, false, "Camera");
    inverse33 (imatrices.xyz_cam, imatrices.cam_xyz);

    for (int c = 0; c < 4; c++) {
        c_white[c] = ri->get_white(c);
    }

    // First we get the "as shot" ("Camera") white balance and store it
    float pre_mul[4];
    // FIXME: get_colorsCoeff not so much used nowadays, when we have calculate_scale_mul() function here
    ri->get_colorsCoeff( pre_mul, scale_mul, c_black, false);//modify  for black level
    camInitialGain = max(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]) / min(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);

    double camwb_red = ri->get_pre_mul(0) / pre_mul[0];
    double camwb_green = ri->get_pre_mul(1) / pre_mul[1];
    double camwb_blue = ri->get_pre_mul(2) / pre_mul[2];
    double cam_r = imatrices.rgb_cam[0][0] * camwb_red + imatrices.rgb_cam[0][1] * camwb_green + imatrices.rgb_cam[0][2] * camwb_blue;
    double cam_g = imatrices.rgb_cam[1][0] * camwb_red + imatrices.rgb_cam[1][1] * camwb_green + imatrices.rgb_cam[1][2] * camwb_blue;
    double cam_b = imatrices.rgb_cam[2][0] * camwb_red + imatrices.rgb_cam[2][1] * camwb_green + imatrices.rgb_cam[2][2] * camwb_blue;
    camera_wb = ColorTemp (cam_r, cam_g, cam_b, 1.); // as shot WB

    ColorTemp ReferenceWB;
    double ref_r, ref_g, ref_b;
    {
        // ...then we re-get the constants but now with auto which gives us better demosaicing and CA auto-correct
        // performance for strange white balance settings (such as UniWB)
        ri->get_colorsCoeff( ref_pre_mul, scale_mul, c_black, true);
        refwb_red = ri->get_pre_mul(0) / ref_pre_mul[0];
        refwb_green = ri->get_pre_mul(1) / ref_pre_mul[1];
        refwb_blue = ri->get_pre_mul(2) / ref_pre_mul[2];
        initialGain = max(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]) / min(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);
        ref_r = imatrices.rgb_cam[0][0] * refwb_red + imatrices.rgb_cam[0][1] * refwb_green + imatrices.rgb_cam[0][2] * refwb_blue;
        ref_g = imatrices.rgb_cam[1][0] * refwb_red + imatrices.rgb_cam[1][1] * refwb_green + imatrices.rgb_cam[1][2] * refwb_blue;
        ref_b = imatrices.rgb_cam[2][0] * refwb_red + imatrices.rgb_cam[2][1] * refwb_green + imatrices.rgb_cam[2][2] * refwb_blue;
        ReferenceWB = ColorTemp (ref_r, ref_g, ref_b, 1.);
    }

    if (settings->verbose) {
        printf("Raw As Shot White balance: temp %f, tint %f\n", camera_wb.getTemp(), camera_wb.getGreen());
        printf("Raw Reference (auto) white balance: temp %f, tint %f, multipliers [%f %f %f | %f %f %f]\n", ReferenceWB.getTemp(), ReferenceWB.getGreen(), ref_r, ref_g, ref_b, refwb_red, refwb_blue, refwb_green);
    }

    /*{
            // Test code: if you want to test a specific white balance
        ColorTemp d50wb = ColorTemp(5000.0, 1.0, 1.0, "Custom");
        double rm,gm,bm,r,g,b;
        d50wb.getMultipliers(r, g, b);
        camwb_red   = imatrices.cam_rgb[0][0]*r + imatrices.cam_rgb[0][1]*g + imatrices.cam_rgb[0][2]*b;
        camwb_green = imatrices.cam_rgb[1][0]*r + imatrices.cam_rgb[1][1]*g + imatrices.cam_rgb[1][2]*b;
        camwb_blue  = imatrices.cam_rgb[2][0]*r + imatrices.cam_rgb[2][1]*g + imatrices.cam_rgb[2][2]*b;
        double pre_mul[3], dmax = 0;
        pre_mul[0] = ri->get_pre_mul(0) / camwb_red;
        pre_mul[1] = ri->get_pre_mul(1) / camwb_green;
        pre_mul[2] = ri->get_pre_mul(2) / camwb_blue;
        for (int c = 0; c < 3; c++) {
            if (dmax < pre_mul[c])
                dmax = pre_mul[c];
                }
                for (int c = 0; c < 3; c++) {
            pre_mul[c] /= dmax;
                }
                camwb_red *= dmax;
                camwb_green *= dmax;
                camwb_blue *= dmax;
                for (int c = 0; c < 3; c++) {
            int sat = ri->get_white(c) - ri->get_cblack(c);
            scale_mul[c] = pre_mul[c] * 65535.0 / sat;
                }
                scale_mul[3] = pre_mul[1] * 65535.0 / (ri->get_white(3) - ri->get_cblack(3));
                initialGain = 1.0 / min(pre_mul[0], pre_mul[1], pre_mul[2]);
    }*/

    for(unsigned int i = 0;i < numFrames; ++i) {
        riFrames[i]->set_prefilters();
    }


    // Load complete Exif informations
    std::unique_ptr<RawMetaDataLocation> rml(new RawMetaDataLocation (ri->get_exifBase(), ri->get_ciffBase(), ri->get_ciffLen()));
    idata = new FramesData (fname, std::move(rml));
    idata->setDCRawFrameCount (numFrames);

    green(W, H);
    red(W, H);
    blue(W, H);
    //hpmap = allocArray<char>(W, H);

    if (plistener) {
        plistener->setProgress (1.0);
    }

    plistener = nullptr; // This must be reset, because only load() is called through progressConnector
    t2.set();

    if( settings->verbose ) {
        printf("Load %s: %d usec\n", fname.c_str(), t2.etime(t1));
    }

    return 0; // OK!
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::preprocess  (const RAWParams &raw, const LensProfParams &lensProf, const CoarseTransformParams& coarse, bool prepareDenoise)
{
//    BENCHFUN
    MyTime t1, t2;
    t1.set();

    Glib::ustring newDF = raw.dark_frame;
    RawImage *rid = nullptr;

    if (!raw.df_autoselect) {
        if( !raw.dark_frame.empty()) {
            rid = dfm.searchDarkFrame( raw.dark_frame );
        }
    } else {
        rid = dfm.searchDarkFrame(idata->getMake(), idata->getModel(), idata->getISOSpeed(), idata->getShutterSpeed(), idata->getDateTimeAsTS());
    }

    if( rid && settings->verbose) {
        printf( "Subtracting Darkframe:%s\n", rid->get_filename().c_str());
    }

    PixelsMap *bitmapBads = nullptr;

    int totBP = 0; // Hold count of bad pixels to correct

    if(ri->zeroIsBad()) { // mark all pixels with value zero as bad, has to be called before FF and DF. dcraw sets this flag only for some cameras (mainly Panasonic and Leica)
        bitmapBads = new PixelsMap(W, H);
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:totBP) schedule(dynamic,16)
#endif

        for(int i = 0; i < H; i++)
            for(int j = 0; j < W; j++) {
                if(ri->data[i][j] == 0.f) {
                    bitmapBads->set(j, i);
                    totBP++;
                }
            }

        if( settings->verbose) {
            printf( "%d pixels with value zero marked as bad pixels\n", totBP);
        }
    }

    //FLATFIELD start
    RawImage *rif = nullptr;

    if (!raw.ff_AutoSelect) {
        if( !raw.ff_file.empty()) {
            rif = ffm.searchFlatField( raw.ff_file );
        }
    } else {
        rif = ffm.searchFlatField( idata->getMake(), idata->getModel(), idata->getLens(), idata->getFocalLen(), idata->getFNumber(), idata->getDateTimeAsTS());
    }


    bool hasFlatField = (rif != nullptr);

    if( hasFlatField && settings->verbose) {
        printf( "Flat Field Correction:%s\n", rif->get_filename().c_str());
    }

    if(numFrames == 4) {
        int bufferNumber = 0;
        for(unsigned int i=0; i<4; ++i) {
            if(i==currFrame) {
                copyOriginalPixels(raw, ri, rid, rif, rawData);
                rawDataFrames[i] = &rawData;
            } else {
                if(!rawDataBuffer[bufferNumber]) {
                    rawDataBuffer[bufferNumber] = new array2D<float>;
                }
                rawDataFrames[i] = rawDataBuffer[bufferNumber];
                ++bufferNumber;
                copyOriginalPixels(raw, riFrames[i], rid, rif, *rawDataFrames[i]);
            }
        }
    } else {
        copyOriginalPixels(raw, ri, rid, rif, rawData);
    }
    //FLATFIELD end


    // Always correct camera badpixels from .badpixels file
    std::vector<badPix> *bp = dfm.getBadPixels( ri->get_maker(), ri->get_model(), idata->getSerialNumber() );

    if( bp ) {
        if(!bitmapBads) {
            bitmapBads = new PixelsMap(W, H);
        }

        totBP += bitmapBads->set( *bp );

        if( settings->verbose ) {
            std::cout << "Correcting " << bp->size() << " pixels from .badpixels" << std::endl;
        }
    }

    // If darkframe selected, correct hotpixels found on darkframe
    bp = nullptr;

    if( raw.df_autoselect ) {
        bp = dfm.getHotPixels(idata->getMake(), idata->getModel(), idata->getISOSpeed(), idata->getShutterSpeed(), idata->getDateTimeAsTS());
    } else if( !raw.dark_frame.empty() ) {
        bp = dfm.getHotPixels( raw.dark_frame );
    }

    if(bp) {
        if(!bitmapBads) {
            bitmapBads = new PixelsMap(W, H);
        }

        totBP += bitmapBads->set( *bp );

        if( settings->verbose && !bp->empty()) {
            std::cout << "Correcting " << bp->size() << " hotpixels from darkframe" << std::endl;
        }
    }

    if(numFrames == 4) {
        for(int i=0; i<4; ++i) {
            scaleColors( 0, 0, W, H, raw, *rawDataFrames[i]);
        }
    } else {
        scaleColors( 0, 0, W, H, raw, rawData); //+ + raw parameters for black level(raw.blackxx)
    }

    // Correct vignetting of lens profile
    if (!hasFlatField && lensProf.useVign) {
        std::unique_ptr<LensCorrection> pmap;
        if (lensProf.useLensfun()) {
            pmap = LFDatabase::findModifier(lensProf, idata, W, H, coarse, -1);
        } else {
            const std::shared_ptr<LCPProfile> pLCPProf = LCPStore::getInstance()->getProfile(lensProf.lcpFile);

            if (pLCPProf) { // don't check focal length to allow distortion correction for lenses without chip, also pass dummy focal length 1 in case of 0
                pmap.reset(new LCPMapper(pLCPProf, max(idata->getFocalLen(), 1.0), idata->getFocalLen35mm(), idata->getFocusDist(), idata->getFNumber(), true, false, W, H, coarse, -1));
            }
        }

        if (pmap) {
            LensCorrection &map = *pmap;
            if (ri->getSensorType() == ST_BAYER || ri->getSensorType() == ST_FUJI_XTRANS || ri->get_colors() == 1) {
                if(numFrames == 4) {
                    for(int i = 0; i < 4; ++i) {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                        for (int y = 0; y < H; y++) {
                            map.processVignetteLine(W, y, (*rawDataFrames[i])[y]);
                        }
                    }
                } else {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int y = 0; y < H; y++) {
                        map.processVignetteLine(W, y, rawData[y]);
                    }
                }
            } else if(ri->get_colors() == 3) {
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < H; y++) {
                    map.processVignetteLine3Channels(W, y, rawData[y]);
                }
            }
        }
    }

    defGain = 0.0;//log(initialGain) / log(2.0);

    if ( ri->getSensorType() == ST_BAYER && (raw.hotPixelFilter > 0 || raw.deadPixelFilter > 0)) {
        if (plistener) {
            plistener->setProgressStr ("Hot/Dead Pixel Filter...");
            plistener->setProgress (0.0);
        }

        if(!bitmapBads) {
            bitmapBads = new PixelsMap(W, H);
        }

        int nFound = findHotDeadPixels( *bitmapBads, raw.hotdeadpix_thresh, raw.hotPixelFilter, raw.deadPixelFilter );
        totBP += nFound;

        if( settings->verbose && nFound > 0) {
            printf( "Correcting %d hot/dead pixels found inside image\n", nFound );
        }
    }

    // check if it is an olympus E camera or green equilibration is enabled. If yes, compute G channel pre-compensation factors
    if ( ri->getSensorType() == ST_BAYER && (raw.bayersensor.greenthresh || (((idata->getMake().size() >= 7 && idata->getMake().substr(0, 7) == "OLYMPUS" && idata->getModel()[0] == 'E') || (idata->getMake().size() >= 9 && idata->getMake().substr(0, 9) == "Panasonic")) && raw.bayersensor.method != RAWParams::BayerSensor::methodstring[ RAWParams::BayerSensor::vng4])) ) {
        // global correction
        if(numFrames == 4) {
            for(int i = 0; i < 4; ++i) {
                green_equilibrate_global(*rawDataFrames[i]);
            }
        } else {
            green_equilibrate_global(rawData);
        }
    }

    if ( ri->getSensorType() == ST_BAYER && raw.bayersensor.greenthresh > 0) {
        if (plistener) {
            plistener->setProgressStr ("Green equilibrate...");
            plistener->setProgress (0.0);
        }

        if(numFrames == 4) {
            for(int i = 0; i < 4; ++i) {
                green_equilibrate(0.01 * raw.bayersensor.greenthresh, *rawDataFrames[i]);
            }
        } else {
            green_equilibrate(0.01 * raw.bayersensor.greenthresh, rawData);
        }
    }


    if( totBP ) {
        if ( ri->getSensorType() == ST_BAYER ) {
            if(numFrames == 4) {
                for(int i = 0; i < 4; ++i) {
                    interpolateBadPixelsBayer( *bitmapBads, *rawDataFrames[i] );
                }
            } else {
                interpolateBadPixelsBayer( *bitmapBads, rawData );
            }
        } else if ( ri->getSensorType() == ST_FUJI_XTRANS ) {
            interpolateBadPixelsXtrans( *bitmapBads );
        } else {
            interpolateBadPixelsNColours( *bitmapBads, ri->get_colors() );
        }
    }

    if ( ri->getSensorType() == ST_BAYER && raw.bayersensor.linenoise > 0 ) {
        if (plistener) {
            plistener->setProgressStr ("Line Denoise...");
            plistener->setProgress (0.0);
        }

        cfa_linedn(0.00002 * (raw.bayersensor.linenoise));
    }

    if ( (raw.ca_autocorrect || fabs(raw.cared) > 0.001 || fabs(raw.cablue) > 0.001) && ri->getSensorType() == ST_BAYER ) { // Auto CA correction disabled for X-Trans, for now...
        if (plistener) {
            plistener->setProgressStr ("CA Auto Correction...");
            plistener->setProgress (0.0);
        }
        if(numFrames == 4) {
            for(int i=0; i<4; ++i) {
                CA_correct_RT(raw.ca_autocorrect, raw.cared, raw.cablue, 8.0, *rawDataFrames[i]);
            }
        } else {
            CA_correct_RT(raw.ca_autocorrect, raw.cared, raw.cablue, 8.0, rawData);
        }
    }

    if ( raw.expos != 1 ) {
        if(numFrames == 4) {
            for(int i = 0; i < 4; ++i) {
                processRawWhitepoint(raw.expos, raw.preser, *rawDataFrames[i]);
            }
        } else {
            processRawWhitepoint(raw.expos, raw.preser, rawData);
        }
    }

    if(prepareDenoise && dirpyrdenoiseExpComp == INFINITY) {
        LUTu aehist;
        int aehistcompr;
        double clip = 0;
        int brightness, contrast, black, hlcompr, hlcomprthresh;
        getAutoExpHistogram (aehist, aehistcompr);
        ImProcFunctions::getAutoExp (aehist, aehistcompr, clip, dirpyrdenoiseExpComp, brightness, contrast, black, hlcompr, hlcomprthresh);
    }

    t2.set();

    if( settings->verbose ) {
        printf("Preprocessing: %d usec\n", t2.etime(t1));
    }

    if(bitmapBads) {
        delete bitmapBads;
    }

    rawDirty = true;
    return;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::demosaic(const RAWParams &raw)
{
    MyTime t1, t2;
    t1.set();

    if (ri->getSensorType() == ST_BAYER) {
        if ( raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::hphd] ) {
            hphd_demosaic ();
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::vng4] ) {
            vng4_demosaic ();
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::ahd] ) {
            ahd_demosaic ();
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::amaze] ) {
            amaze_demosaic_RT (0, 0, W, H, rawData, red, green, blue);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::pixelshift] ) {
            pixelshift(0, 0, W, H, raw.bayersensor, currFrame, ri->get_model(), raw.expos);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::dcb] ) {
            dcb_demosaic(raw.bayersensor.dcb_iterations, raw.bayersensor.dcb_enhance);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::eahd]) {
            eahd_demosaic ();
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::igv]) {
            igv_interpolate(W, H);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::lmmse]) {
            lmmse_interpolate_omp(W, H, rawData, red, green, blue, raw.bayersensor.lmmse_iterations);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::fast] ) {
            fast_demosaic (0, 0, W, H);
        } else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::mono] ) {
            nodemosaic(true);
        } else {
            nodemosaic(false);
        }

        //if (raw.all_enhance) refinement_lassus();

    } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
        if (raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::fast] ) {
            fast_xtrans_interpolate();
        } else if (raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::onePass]) {
            xtrans_interpolate(1, false);
        } else if (raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::threePass] ) {
            xtrans_interpolate(3, true);
        } else if(raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::mono] ) {
            nodemosaic(true);
        } else {
            nodemosaic(false);
        }
    } else if (ri->get_colors() == 1) {
        // Monochrome
        nodemosaic(true);
    }

    t2.set();


    rgbSourceModified = false;


    if( settings->verbose ) {
        if (getSensorType() == ST_BAYER) {
            printf("Demosaicing Bayer data: %s - %d usec\n", raw.bayersensor.method.c_str(), t2.etime(t1));
        } else if (getSensorType() == ST_FUJI_XTRANS) {
            printf("Demosaicing X-Trans data: %s - %d usec\n", raw.xtranssensor.method.c_str(), t2.etime(t1));
        }
    }
}


//void RawImageSource::retinexPrepareBuffers(ColorManagementParams cmp, RetinexParams retinexParams, multi_array2D<float, 3> &conversionBuffer, LUTu &lhist16RETI)
void RawImageSource::retinexPrepareBuffers(ColorManagementParams cmp, const RetinexParams &retinexParams, multi_array2D<float, 4> &conversionBuffer, LUTu &lhist16RETI)
{
    bool useHsl = (retinexParams.retinexcolorspace == "HSLLOG" || retinexParams.retinexcolorspace == "HSLLIN");
    conversionBuffer[0] (W - 2 * border, H - 2 * border);
    conversionBuffer[1] (W - 2 * border, H - 2 * border);
    conversionBuffer[2] (W - 2 * border, H - 2 * border);
    conversionBuffer[3] (W - 2 * border, H - 2 * border);

    LUTf *retinexgamtab = nullptr;//gamma before and after Retinex to restore tones
    LUTf lutTonereti;

    if(retinexParams.gammaretinex == "low") {
        retinexgamtab = &(Color::gammatab_115_2);
    } else if(retinexParams.gammaretinex == "mid") {
        retinexgamtab = &(Color::gammatab_13_2);
    } else if(retinexParams.gammaretinex == "hig") {
        retinexgamtab = &(Color::gammatab_145_3);
    } else if(retinexParams.gammaretinex == "fre") {
        GammaValues g_a;
        double pwr = 1.0 / retinexParams.gam;
        double gamm = retinexParams.gam;
        double ts = retinexParams.slope;
        double gamm2 = retinexParams.gam;

        if(gamm2 < 1.) {
            std::swap(pwr, gamm);
        }

        int mode = 0, imax = 0;
        Color::calcGamma(pwr, ts, mode, imax, g_a); // call to calcGamma with selected gamma and slope

        //    printf("g_a0=%f g_a1=%f g_a2=%f g_a3=%f g_a4=%f\n", g_a0,g_a1,g_a2,g_a3,g_a4);
        double start;
        double add;

        if(gamm2 < 1.) {
            start = g_a[2];
            add = g_a[4];
        } else {
            start = g_a[3];
            add = g_a[4];
        }

        double mul = 1. + g_a[4];

        lutTonereti(65536);

        for (int i = 0; i < 65536; i++) {
            double val = (i) / 65535.;
            double x;

            if(gamm2 < 1.) {
                x = Color::igammareti (val, gamm, start, ts, mul , add);
            } else {
                x = Color::gammareti (val, gamm, start, ts, mul , add);
            }

            lutTonereti[i] = CLIP(x * 65535.);// CLIP avoid in some case extra values
        }

        retinexgamtab = &lutTonereti;
    }

    /*
    //test with amsterdam.pef and other files
    float rr,gg,bb;
    rr=red[50][2300];
    gg=green[50][2300];
    bb=blue[50][2300];
    printf("rr=%f gg=%f bb=%f \n",rr,gg,bb);
    rr=red[1630][370];
    gg=green[1630][370];
    bb=blue[1630][370];
    printf("rr1=%f gg1=%f bb1=%f \n",rr,gg,bb);
    rr=red[380][1630];
    gg=green[380][1630];
    bb=blue[380][1630];
    printf("rr2=%f gg2=%f bb2=%f \n",rr,gg,bb);
    */
    /*
    if(retinexParams.highlig < 100 && retinexParams.retinexMethod == "highliplus") {//try to recover magenta...very difficult !
        float hig = ((float)retinexParams.highlig)/100.f;
        float higgb = ((float)retinexParams.grbl)/100.f;

    #ifdef _OPENMP
            #pragma omp parallel for
    #endif
            for (int i = border; i < H - border; i++ ) {
                for (int j = border; j < W - border; j++ ) {
                    float R_,G_,B_;
                    R_=red[i][j];
                    G_=green[i][j];
                    B_=blue[i][j];

                    //empirical method to find highlight magenta with no conversion RGB and no white balance
                    //red = master   Gr and Bl default higgb=0.5
         //           if(R_>65535.f*hig  && G_ > 65535.f*higgb && B_ > 65535.f*higgb) conversionBuffer[3][i - border][j - border] = R_;
          //          else conversionBuffer[3][i - border][j - border] = 0.f;
                }
            }
    }
    */
    if(retinexParams.gammaretinex != "none" && retinexParams.str != 0 && retinexgamtab) {//gamma

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = border; i < H - border; i++ ) {
            for (int j = border; j < W - border; j++ ) {
                float R_, G_, B_;
                R_ = red[i][j];
                G_ = green[i][j];
                B_ = blue[i][j];

                red[i][j] = (*retinexgamtab)[R_];
                green[i][j] = (*retinexgamtab)[G_];
                blue[i][j] = (*retinexgamtab)[B_];
            }
        }
    }

    if(useHsl) {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            // one LUT per thread
            LUTu lhist16RETIThr;

            if(lhist16RETI)
            {
                lhist16RETIThr(lhist16RETI.getSize());
                lhist16RETIThr.clear();
            }

#ifdef __SSE2__
            vfloat c32768 = F2V(32768.f);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = border; i < H - border; i++ )
            {
                int j = border;
#ifdef __SSE2__

                for (; j < W - border - 3; j += 4) {
                    vfloat H, S, L;
                    Color::rgb2hsl(LVFU(red[i][j]), LVFU(green[i][j]), LVFU(blue[i][j]), H, S, L);
                    STVFU(conversionBuffer[0][i - border][j - border], H);
                    STVFU(conversionBuffer[1][i - border][j - border], S);
                    L *= c32768;
                    STVFU(conversionBuffer[2][i - border][j - border], L);
                    STVFU(conversionBuffer[3][i - border][j - border], H);

                    if(lhist16RETI) {
                        for(int p = 0; p < 4; p++) {
                            int pos = ( conversionBuffer[2][i - border][j - border + p]);//histogram in curve HSL
                            lhist16RETIThr[pos]++;
                        }
                    }
                }

#endif

                for (; j < W - border; j++) {
                    float L;
                    //rgb=>lab
                    Color::rgb2hslfloat(red[i][j], green[i][j], blue[i][j], conversionBuffer[0][i - border][j - border], conversionBuffer[1][i - border][j - border], L);
                    L *= 32768.f;
                    conversionBuffer[2][i - border][j - border] = L;

                    if(lhist16RETI) {
                        int pos = L;
                        lhist16RETIThr[pos]++;
                    }
                }
            }

#ifdef _OPENMP
            #pragma omp critical
            {
                if(lhist16RETI)
                {
                    lhist16RETI += lhist16RETIThr; // Add per Thread LUT to global LUT
                }
            }
#endif

        }
    } else {
        TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (cmp.working);
        const float wp[3][3] = {
            {static_cast<float>(wprof[0][0]), static_cast<float>(wprof[0][1]), static_cast<float>(wprof[0][2])},
            {static_cast<float>(wprof[1][0]), static_cast<float>(wprof[1][1]), static_cast<float>(wprof[1][2])},
            {static_cast<float>(wprof[2][0]), static_cast<float>(wprof[2][1]), static_cast<float>(wprof[2][2])}
        };

        // Conversion rgb -> lab is hard to vectorize because it uses a lut (that's not the main problem)
        // and it uses a condition inside XYZ2Lab which is almost impossible to vectorize without making it slower...
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            // one LUT per thread
            LUTu lhist16RETIThr;

            if(lhist16RETI) {
                lhist16RETIThr(lhist16RETI.getSize());
                lhist16RETIThr.clear();
            }

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int i = border; i < H - border; i++ )
                for (int j = border; j < W - border; j++) {
                    float X, Y, Z, L, aa, bb;
                    //rgb=>lab
                    Color::rgbxyz(red[i][j], green[i][j], blue[i][j], X, Y, Z, wp);
                    //convert Lab
                    Color::XYZ2Lab(X, Y, Z, L, aa, bb);
                    conversionBuffer[0][i - border][j - border] = aa;
                    conversionBuffer[1][i - border][j - border] = bb;
                    conversionBuffer[2][i - border][j - border] = L;
                    conversionBuffer[3][i - border][j - border] = xatan2f(bb, aa);

//                   if(R_>40000.f  && G_ > 30000.f && B_ > 30000.f) conversionBuffer[3][i - border][j - border] = R_;
//                   else conversionBuffer[3][i - border][j - border] = 0.f;
                    if(lhist16RETI) {
                        int pos = L;
                        lhist16RETIThr[pos]++;//histogram in Curve Lab
                    }
                }

#ifdef _OPENMP
            #pragma omp critical
            {
                if(lhist16RETI) {
                    lhist16RETI += lhist16RETIThr; // Add per Thread LUT to global LUT
                }
            }
#endif

        }
    }



}

void RawImageSource::retinexPrepareCurves(const RetinexParams &retinexParams, LUTf &cdcurve, LUTf &mapcurve, RetinextransmissionCurve &retinextransmissionCurve, RetinexgaintransmissionCurve &retinexgaintransmissionCurve, bool &retinexcontlutili, bool &mapcontlutili, bool &useHsl, LUTu & lhist16RETI, LUTu & histLRETI)
{
    useHsl = (retinexParams.retinexcolorspace == "HSLLOG" || retinexParams.retinexcolorspace == "HSLLIN");

    if(useHsl) {
        CurveFactory::curveDehaContL (retinexcontlutili, retinexParams.cdHcurve, cdcurve, 1, lhist16RETI, histLRETI);
    } else {
        CurveFactory::curveDehaContL (retinexcontlutili, retinexParams.cdcurve, cdcurve, 1, lhist16RETI, histLRETI);
    }

    CurveFactory::mapcurve (mapcontlutili, retinexParams.mapcurve, mapcurve, 1, lhist16RETI, histLRETI);

    retinexParams.getCurves(retinextransmissionCurve, retinexgaintransmissionCurve);
}

void RawImageSource::retinex(ColorManagementParams cmp, const RetinexParams &deh, ToneCurveParams Tc, LUTf & cdcurve, LUTf & mapcurve, const RetinextransmissionCurve & dehatransmissionCurve, const RetinexgaintransmissionCurve & dehagaintransmissionCurve, multi_array2D<float, 4> &conversionBuffer, bool dehacontlutili, bool mapcontlutili, bool useHsl, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax, LUTu &histLRETI)
{
    MyTime t4, t5;
    t4.set();

    if (settings->verbose) {
        printf ("Applying Retinex\n");
    }

    LUTf lutToneireti;
    lutToneireti(65536);

    LUTf *retinexigamtab = nullptr;//gamma before and after Retinex to restore tones

    if(deh.gammaretinex == "low") {
        retinexigamtab = &(Color::igammatab_115_2);
    } else if(deh.gammaretinex == "mid") {
        retinexigamtab = &(Color::igammatab_13_2);
    } else if(deh.gammaretinex == "hig") {
        retinexigamtab = &(Color::igammatab_145_3);
    } else if(deh.gammaretinex == "fre") {
        GammaValues g_a;
        double pwr = 1.0 / deh.gam;
        double gamm = deh.gam;
        double gamm2 = gamm;
        double ts = deh.slope;
        int mode = 0, imax = 0;

        if(gamm2 < 1.) {
            std::swap(pwr, gamm);
        }

        Color::calcGamma(pwr, ts, mode, imax, g_a); // call to calcGamma with selected gamma and slope

        double mul = 1. + g_a[4];
        double add;
        double start;

        if(gamm2 < 1.) {
            start = g_a[3];
            add = g_a[3];
        } else {
            add = g_a[4];
            start = g_a[2];
        }

        //    printf("g_a0=%f g_a1=%f g_a2=%f g_a3=%f g_a4=%f\n", g_a0,g_a1,g_a2,g_a3,g_a4);
        for (int i = 0; i < 65536; i++) {
            double val = (i) / 65535.;
            double x;

            if(gamm2 < 1.) {
                x = Color::gammareti (val, gamm, start, ts, mul , add);
            } else {
                x = Color::igammareti (val, gamm, start, ts, mul , add);
            }

            lutToneireti[i] = CLIP(x * 65535.);
        }

        retinexigamtab = &lutToneireti;
    }

    // We need a buffer with original L data to allow correct blending
    // red, green and blue still have original size of raw, but we can't use the borders
    const int HNew = H - 2 * border;
    const int WNew = W - 2 * border;

    array2D<float> LBuffer (WNew, HNew);
    float **temp = conversionBuffer[2]; // one less dereference
    LUTf dLcurve;
    LUTu hist16RET;

    if(dehacontlutili && histLRETI) {
        hist16RET(32768);
        hist16RET.clear();
        histLRETI.clear();
        dLcurve(32768);
    }

    FlatCurve* chcurve = nullptr;//curve c=f(H)
    bool chutili = false;

    if (deh.enabled && deh.retinexMethod == "highli") {
        chcurve = new FlatCurve(deh.lhcurve);

        if (!chcurve || chcurve->isIdentity()) {
            if (chcurve) {
                delete chcurve;
                chcurve = nullptr;
            }
        } else {
            chutili = true;
        }
    }



#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        // one LUT per thread
        LUTu hist16RETThr;

        if(hist16RET) {
            hist16RETThr(hist16RET.getSize());
            hist16RETThr.clear();
        }

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H - 2 * border; i++ )
            if(dehacontlutili)
                for (int j = 0; j < W - 2 * border; j++) {
                    LBuffer[i][j] = cdcurve[2.f * temp[i][j]] / 2.f;

                    if(histLRETI) {
                        int pos = LBuffer[i][j];
                        hist16RETThr[pos]++; //histogram in Curve
                    }
                }
            else
                for (int j = 0; j < W - 2 * border; j++) {
                    LBuffer[i][j] = temp[i][j];
                }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            if(hist16RET) {
                hist16RET += hist16RETThr; // Add per Thread LUT to global LUT
            }
        }
    }

    if(hist16RET) {//update histogram
        // TODO : When rgbcurvesspeedup branch is merged into master, replace this by the following 1-liner
        // hist16RET.compressTo(histLRETI);
        // also remove declaration and init of dLcurve some lines above then and finally remove this comment :)
        for (int i = 0; i < 32768; i++) {
            float val = (double)i / 32767.0;
            dLcurve[i] = val;
        }

        for (int i = 0; i < 32768; i++) {
            float hval = dLcurve[i];
            int hi = (int)(255.0f * hval);
            histLRETI[hi] += hist16RET[i];
        }
    }

    MSR(LBuffer, conversionBuffer[2], conversionBuffer[3], mapcurve, mapcontlutili, WNew, HNew, deh, dehatransmissionCurve, dehagaintransmissionCurve, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);

    if(useHsl) {
        if(chutili) {
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int i = border; i < H - border; i++ ) {
                int j = border;

                for (; j < W - border; j++) {

                    float valp = (chcurve->getVal(conversionBuffer[3][i - border][j - border]) - 0.5f);
                    conversionBuffer[1][i - border][j - border] *= (1.f + 2.f * valp);

                }
            }
        }

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = border; i < H - border; i++ ) {
            int j = border;
#ifdef __SSE2__
            vfloat c32768 = F2V(32768.f);

            for (; j < W - border - 3; j += 4) {
                vfloat R, G, B;
                Color::hsl2rgb(LVFU(conversionBuffer[0][i - border][j - border]), LVFU(conversionBuffer[1][i - border][j - border]), LVFU(LBuffer[i - border][j - border]) / c32768, R, G, B);

                STVFU(red[i][j], R);
                STVFU(green[i][j], G);
                STVFU(blue[i][j], B);
            }

#endif

            for (; j < W - border; j++) {
                Color::hsl2rgbfloat(conversionBuffer[0][i - border][j - border], conversionBuffer[1][i - border][j - border], LBuffer[i - border][j - border] / 32768.f, red[i][j], green[i][j], blue[i][j]);
            }
        }

    } else {
        TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (cmp.working);

        double wip[3][3] = {
            {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
            {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
            {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
        };
        // gamut control only in Lab mode
        const bool highlight = Tc.hrenabled;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            // we need some line buffers to precalculate some expensive stuff using SSE
            float atan2Buffer[W] ALIGNED16;
            float sqrtBuffer[W] ALIGNED16;
            float sincosxBuffer[W] ALIGNED16;
            float sincosyBuffer[W] ALIGNED16;
            const vfloat c327d68v = F2V(327.68);
            const vfloat onev = F2V(1.f);
#endif // __SSE2__
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int i = border; i < H - border; i++ ) {
#ifdef __SSE2__
                // vectorized precalculation
                {
                    int j = border;

                    for (; j < W - border - 3; j += 4)
                    {
                        vfloat av = LVFU(conversionBuffer[0][i - border][j - border]);
                        vfloat bv = LVFU(conversionBuffer[1][i - border][j - border]);
                        vfloat chprovv = vsqrtf(SQRV(av) + SQRV(bv));
                        STVF(sqrtBuffer[j - border], chprovv / c327d68v);
                        vfloat HHv = xatan2f(bv, av);
                        STVF(atan2Buffer[j - border], HHv);
                        av /= chprovv;
                        bv /= chprovv;
                        vmask selMask = vmaskf_eq(chprovv, ZEROV);
                        STVF(sincosyBuffer[j - border], vself(selMask, onev, av));
                        STVF(sincosxBuffer[j - border], vselfnotzero(selMask, bv));
                    }

                    for (; j < W - border; j++)
                    {
                        float aa = conversionBuffer[0][i - border][j - border];
                        float bb = conversionBuffer[1][i - border][j - border];
                        float Chprov1 = sqrt(SQR(aa) + SQR(bb)) / 327.68f;
                        sqrtBuffer[j - border] = Chprov1;
                        float HH = xatan2f(bb, aa);
                        atan2Buffer[j - border] = HH;

                        if(Chprov1 == 0.0f) {
                            sincosyBuffer[j - border] = 1.f;
                            sincosxBuffer[j - border] = 0.0f;
                        } else {
                            sincosyBuffer[j - border] = aa / (Chprov1 * 327.68f);
                            sincosxBuffer[j - border] = bb / (Chprov1 * 327.68f);
                        }
                    }
                }
#endif // __SSE2__

                for (int j = border; j < W - border; j++) {
                    float Lprov1 = (LBuffer[i - border][j - border]) / 327.68f;
#ifdef __SSE2__
                    float Chprov1 = sqrtBuffer[j - border];
                    float  HH = atan2Buffer[j - border];
                    float2 sincosval;
                    sincosval.x = sincosxBuffer[j - border];
                    sincosval.y = sincosyBuffer[j - border];

#else
                    float aa = conversionBuffer[0][i - border][j - border];
                    float bb = conversionBuffer[1][i - border][j - border];
                    float Chprov1 = sqrt(SQR(aa) + SQR(bb)) / 327.68f;
                    float  HH = xatan2f(bb, aa);
                    float2 sincosval;// = xsincosf(HH);

                    if(Chprov1 == 0.0f) {
                        sincosval.y = 1.f;
                        sincosval.x = 0.0f;
                    } else {
                        sincosval.y = aa / (Chprov1 * 327.68f);
                        sincosval.x = bb / (Chprov1 * 327.68f);
                    }

#endif

                    if(chutili) {  // c=f(H)
                        float valp = float((chcurve->getVal(Color::huelab_to_huehsv2(HH)) - 0.5f));
                        Chprov1 *= (1.f + 2.f * valp);
                    }

                    float R, G, B;
#ifdef _DEBUG
                    bool neg = false;
                    bool more_rgb = false;
                    //gamut control : Lab values are in gamut
                    Color::gamutLchonly(HH, sincosval, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                    //gamut control : Lab values are in gamut
                    Color::gamutLchonly(HH, sincosval, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif



                    conversionBuffer[0][i - border][j - border] = 327.68f * Chprov1 * sincosval.y;
                    conversionBuffer[1][i - border][j - border] = 327.68f * Chprov1 * sincosval.x;
                    LBuffer[i - border][j - border] = Lprov1 * 327.68f;
                }
            }
        }
        //end gamut control
#ifdef __SSE2__
        vfloat wipv[3][3];

        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++) {
                wipv[i][j] = F2V(wiprof[i][j]);
            }

#endif // __SSE2__
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = border; i < H - border; i++ ) {
            int j = border;
#ifdef __SSE2__

            for (; j < W - border - 3; j += 4) {
                vfloat x_, y_, z_;
                vfloat R, G, B;
                Color::Lab2XYZ(LVFU(LBuffer[i - border][j - border]), LVFU(conversionBuffer[0][i - border][j - border]), LVFU(conversionBuffer[1][i - border][j - border]), x_, y_, z_) ;
                Color::xyz2rgb(x_, y_, z_, R, G, B, wipv);

                STVFU(red[i][j], R);
                STVFU(green[i][j], G);
                STVFU(blue[i][j], B);

            }

#endif

            for (; j < W - border; j++) {
                float x_, y_, z_;
                float R, G, B;
                Color::Lab2XYZ(LBuffer[i - border][j - border], conversionBuffer[0][i - border][j - border], conversionBuffer[1][i - border][j - border], x_, y_, z_) ;
                Color::xyz2rgb(x_, y_, z_, R, G, B, wip);
                red[i][j] = R;
                green[i][j] = G;
                blue[i][j] = B;
            }
        }
    }

    if (chcurve) {
        delete chcurve;
    }

    if(deh.gammaretinex != "none"  && deh.str != 0) { //inverse gamma
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = border; i < H - border; i++ ) {
            for (int j = border; j < W - border; j++ ) {
                float R_, G_, B_;
                R_ = red[i][j];
                G_ = green[i][j];
                B_ = blue[i][j];
                red[i][j] = (*retinexigamtab)[R_];
                green[i][j] = (*retinexigamtab)[G_];
                blue[i][j] = (*retinexigamtab)[B_];
            }
        }
    }

    rgbSourceModified = false; // tricky handling for Color propagation

    t5.set();

    if( settings->verbose ) {
        printf("Retinex=%d usec\n",  t5.etime(t4));
    }

}

void RawImageSource::flushRawData()
{
    if(cache) {
        delete [] cache;
        cache = nullptr;
    }

    if (rawData) {
        rawData(0, 0);
    }
}

void RawImageSource::flushRGB()
{
    if (green) {
        green(0, 0);
    }

    if (red) {
        red(0, 0);
    }

    if (blue) {
        blue(0, 0);
    }
}

void RawImageSource::HLRecovery_Global(ToneCurveParams hrp)
{
    if (hrp.hrenabled && hrp.method == "Color") {
        if(!rgbSourceModified) {
            if (settings->verbose) {
                printf ("Applying Highlight Recovery: Color propagation...\n");
            }

            HLRecovery_inpaint (red, green, blue);
            rgbSourceModified = true;
        }
    }

}


void RawImageSource::processFlatField(const RAWParams &raw, RawImage *riFlatFile, unsigned short black[4])
{
//    BENCHFUN
    float *cfablur = (float (*)) malloc (H * W * sizeof * cfablur);
    int BS = raw.ff_BlurRadius;
    BS += BS & 1;

    //function call to cfabloxblur
    if (raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::v_ff]) {
        cfaboxblur(riFlatFile, cfablur, 2 * BS, 0);
    } else if (raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::h_ff]) {
        cfaboxblur(riFlatFile, cfablur, 0, 2 * BS);
    } else if (raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::vh_ff]) {
        //slightly more complicated blur if trying to correct both vertical and horizontal anomalies
        cfaboxblur(riFlatFile, cfablur, BS, BS);    //first do area blur to correct vignette
    } else { //(raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::area_ff])
        cfaboxblur(riFlatFile, cfablur, BS, BS);
    }

    if(ri->getSensorType() == ST_BAYER) {
        float refcolor[2][2];

        //find centre average values by channel
        for (int m = 0; m < 2; m++)
            for (int n = 0; n < 2; n++) {
                int row = 2 * (H >> 2) + m;
                int col = 2 * (W >> 2) + n;
                int c  = FC(row, col);
                int c4 = ( c == 1 && !(row & 1) ) ? 3 : c;
                refcolor[m][n] = max(0.0f, cfablur[row * W + col] - black[c4]);
            }

        float limitFactor = 1.f;

        if(raw.ff_AutoClipControl) {
//            int clipControlGui = 0;

            for (int m = 0; m < 2; m++)
                for (int n = 0; n < 2; n++) {
                    float maxval = 0.f;
                    int c  = FC(m, n);
                    int c4 = ( c == 1 && !(m & 1) ) ? 3 : c;
#ifdef _OPENMP
                    #pragma omp parallel
#endif
                    {
                        float maxvalthr = 0.f;
#ifdef _OPENMP
                        #pragma omp for
#endif

                        for (int row = 0; row < H - m; row += 2) {
                            for (int col = 0; col < W - n; col += 2) {
                                float tempval = (rawData[row + m][col + n] - black[c4]) * ( refcolor[m][n] / max(1e-5f, cfablur[(row + m) * W + col + n] - black[c4]) );

                                if(tempval > maxvalthr) {
                                    maxvalthr = tempval;
                                }
                            }
                        }

#ifdef _OPENMP
                        #pragma omp critical
#endif
                        {

                            if(maxvalthr > maxval) {
                                maxval = maxvalthr;
                            }

                        }
                    }

                    // now we have the max value for the channel
                    // if it clips, calculate factor to avoid clipping
                    if(maxval + black[c4] >= ri->get_white(c4)) {
                        limitFactor = min(limitFactor, ri->get_white(c4) / (maxval + black[c4]));
                    }
                }

//            clipControlGui = (1.f - limitFactor) * 100.f;           // this value can be used to set the clip control slider in gui
        } else {
            limitFactor = max((float)(100 - raw.ff_clipControl) / 100.f, 0.01f);
        }

        for (int m = 0; m < 2; m++)
            for (int n = 0; n < 2; n++) {
                refcolor[m][n] *= limitFactor;
            }


        unsigned int c[2][2]  = {{FC(0, 0), FC(0, 1)}, {FC(1, 0), FC(1, 1)}};
        unsigned int c4[2][2];
        c4[0][0] = ( c[0][0] == 1) ? 3 : c[0][0];
        c4[0][1] = ( c[0][1] == 1) ? 3 : c[0][1];
        c4[1][0] = c[1][0];
        c4[1][1] = c[1][1];

#ifdef __SSE2__
        vfloat refcolorv[2] = {_mm_set_ps(refcolor[0][1], refcolor[0][0], refcolor[0][1], refcolor[0][0]),
                               _mm_set_ps(refcolor[1][1], refcolor[1][0], refcolor[1][1], refcolor[1][0])
                              };
        vfloat blackv[2] = {_mm_set_ps(black[c4[0][1]], black[c4[0][0]], black[c4[0][1]], black[c4[0][0]]),
                            _mm_set_ps(black[c4[1][1]], black[c4[1][0]], black[c4[1][1]], black[c4[1][0]])
                           };

        vfloat epsv = F2V(1e-5f);
#endif
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int row = 0; row < H; row ++) {
            int col = 0;
#ifdef __SSE2__
            vfloat rowBlackv = blackv[row & 1];
            vfloat rowRefcolorv = refcolorv[row & 1];

            for (; col < W - 3; col += 4) {
                vfloat vignettecorrv = rowRefcolorv / vmaxf(epsv, LVFU(cfablur[(row) * W + col]) - rowBlackv);
                vfloat valv = LVFU(rawData[row][col]);
                valv -= rowBlackv;
                STVFU(rawData[row][col], valv * vignettecorrv + rowBlackv);
            }

#endif

            for (; col < W; col ++) {
                float vignettecorr = refcolor[row & 1][col & 1] / max(1e-5f, cfablur[(row) * W + col] - black[c4[row & 1][col & 1]]);
                rawData[row][col] = (rawData[row][col] - black[c4[row & 1][col & 1]]) * vignettecorr + black[c4[row & 1][col & 1]];
            }
        }
    } else if(ri->getSensorType() == ST_FUJI_XTRANS) {
        float refcolor[3] = {0.f};
        int cCount[3] = {0};

        //find center ave values by channel
        for (int m = -3; m < 3; m++)
            for (int n = -3; n < 3; n++) {
                int row = 2 * (H >> 2) + m;
                int col = 2 * (W >> 2) + n;
                int c  = riFlatFile->XTRANSFC(row, col);
                refcolor[c] += max(0.0f, cfablur[row * W + col] - black[c]);
                cCount[c] ++;
            }

        for(int c = 0; c < 3; c++) {
            refcolor[c] = refcolor[c] / cCount[c];
        }

        float limitFactor = 1.f;

        if(raw.ff_AutoClipControl) {
            // determine maximum calculated value to avoid clipping
//            int clipControlGui = 0;
            float maxval = 0.f;
            // xtrans files have only one black level actually, so we can simplify the code a bit
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                float maxvalthr = 0.f;
#ifdef _OPENMP
                #pragma omp for schedule(dynamic,16) nowait
#endif

                for (int row = 0; row < H; row++) {
                    for (int col = 0; col < W; col++) {
                        float tempval = (rawData[row][col] - black[0]) * ( refcolor[ri->XTRANSFC(row, col)] / max(1e-5f, cfablur[(row) * W + col] - black[0]) );

                        if(tempval > maxvalthr) {
                            maxvalthr = tempval;
                        }
                    }
                }

#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    if(maxvalthr > maxval) {
                        maxval = maxvalthr;
                    }
                }
            }

            // there's only one white level for xtrans
            if(maxval + black[0] > ri->get_white(0)) {
                limitFactor = ri->get_white(0) / (maxval + black[0]);
//                clipControlGui = (1.f - limitFactor) * 100.f;           // this value can be used to set the clip control slider in gui
            }
        } else {
            limitFactor = max((float)(100 - raw.ff_clipControl) / 100.f, 0.01f);
        }


        for(int c = 0; c < 3; c++) {
            refcolor[c] *= limitFactor;
        }

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int row = 0; row < H; row++) {
            for (int col = 0; col < W; col++) {
                int c  = ri->XTRANSFC(row, col);
                float vignettecorr = ( refcolor[c] / max(1e-5f, cfablur[(row) * W + col] - black[c]) );
                rawData[row][col] = (rawData[row][col] - black[c]) * vignettecorr + black[c];
            }
        }
    }

    if (raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::vh_ff]) {
        float *cfablur1 = (float (*)) malloc (H * W * sizeof * cfablur1);
        float *cfablur2 = (float (*)) malloc (H * W * sizeof * cfablur2);
        //slightly more complicated blur if trying to correct both vertical and horizontal anomalies
        cfaboxblur(riFlatFile, cfablur1, 0, 2 * BS); //now do horizontal blur
        cfaboxblur(riFlatFile, cfablur2, 2 * BS, 0); //now do vertical blur

        if(ri->getSensorType() == ST_BAYER) {
            unsigned int c[2][2]  = {{FC(0, 0), FC(0, 1)}, {FC(1, 0), FC(1, 1)}};
            unsigned int c4[2][2];
            c4[0][0] = ( c[0][0] == 1) ? 3 : c[0][0];
            c4[0][1] = ( c[0][1] == 1) ? 3 : c[0][1];
            c4[1][0] = c[1][0];
            c4[1][1] = c[1][1];

#ifdef __SSE2__
            vfloat blackv[2] = {_mm_set_ps(black[c4[0][1]], black[c4[0][0]], black[c4[0][1]], black[c4[0][0]]),
                                _mm_set_ps(black[c4[1][1]], black[c4[1][0]], black[c4[1][1]], black[c4[1][0]])
                               };

            vfloat epsv = F2V(1e-5f);
#endif
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int row = 0; row < H; row ++) {
                int col = 0;
#ifdef __SSE2__
                vfloat rowBlackv = blackv[row & 1];

                for (; col < W - 3; col += 4) {
                    vfloat linecorrv = SQRV(vmaxf(epsv, LVFU(cfablur[row * W + col]) - rowBlackv)) /
                                       (vmaxf(epsv, LVFU(cfablur1[row * W + col]) - rowBlackv) * vmaxf(epsv, LVFU(cfablur2[row * W + col]) - rowBlackv));
                    vfloat valv = LVFU(rawData[row][col]);
                    valv -= rowBlackv;
                    STVFU(rawData[row][col], valv * linecorrv + rowBlackv);
                }

#endif

                for (; col < W; col ++) {
                    float linecorr = SQR(max(1e-5f, cfablur[row * W + col] - black[c4[row & 1][col & 1]])) /
                                     (max(1e-5f, cfablur1[row * W + col] - black[c4[row & 1][col & 1]]) * max(1e-5f, cfablur2[row * W + col] - black[c4[row & 1][col & 1]])) ;
                    rawData[row][col] = (rawData[row][col] - black[c4[row & 1][col & 1]]) * linecorr + black[c4[row & 1][col & 1]];
                }
            }
        } else if(ri->getSensorType() == ST_FUJI_XTRANS) {
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int row = 0; row < H; row++) {
                for (int col = 0; col < W; col++) {
                    int c  = ri->XTRANSFC(row, col);
                    float hlinecorr = (max(1e-5f, cfablur[(row) * W + col] - black[c]) / max(1e-5f, cfablur1[(row) * W + col] - black[c]) );
                    float vlinecorr = (max(1e-5f, cfablur[(row) * W + col] - black[c]) / max(1e-5f, cfablur2[(row) * W + col] - black[c]) );
                    rawData[row][col] = ((rawData[row][col] - black[c]) * hlinecorr * vlinecorr + black[c]);
                }
            }

        }

        free (cfablur1);
        free (cfablur2);
    }

    free (cfablur);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/* Copy original pixel data and
 * subtract dark frame (if present) from current image and apply flat field correction (if present)
 */
void RawImageSource::copyOriginalPixels(const RAWParams &raw, RawImage *src, RawImage *riDark, RawImage *riFlatFile, array2D<float> &rawData )
{
    // TODO: Change type of black[] to float to avoid conversions
    unsigned short black[4] = {
        (unsigned short)ri->get_cblack(0), (unsigned short)ri->get_cblack(1),
        (unsigned short)ri->get_cblack(2), (unsigned short)ri->get_cblack(3)
    };

    if (ri->getSensorType() == ST_BAYER || ri->getSensorType() == ST_FUJI_XTRANS) {
        if (!rawData) {
            rawData(W, H);
        }

        if (riDark && W == riDark->get_width() && H == riDark->get_height()) { // This works also for xtrans-sensors, because black[0] to black[4] are equal for these
            for (int row = 0; row < H; row++) {
                for (int col = 0; col < W; col++) {
                    int c  = FC(row, col);
                    int c4 = ( c == 1 && !(row & 1) ) ? 3 : c;
                    rawData[row][col] = max(src->data[row][col] + black[c4] - riDark->data[row][col], 0.0f);
                }
            }
        } else {
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int row = 0; row < H; row++) {
                for (int col = 0; col < W; col++) {
                    rawData[row][col] = src->data[row][col];
                }
            }
        }


        if (riFlatFile && W == riFlatFile->get_width() && H == riFlatFile->get_height()) {
            processFlatField(raw, riFlatFile, black);
        }  // flatfield
    } else if (ri->get_colors() == 1) {
        // Monochrome
        if (!rawData) {
            rawData(W, H);
        }

        if (riDark && W == riDark->get_width() && H == riDark->get_height()) {
            for (int row = 0; row < H; row++) {
                for (int col = 0; col < W; col++) {
                    rawData[row][col] = max(src->data[row][col] + black[0] - riDark->data[row][col], 0.0f);
                }
            }
        } else {
            for (int row = 0; row < H; row++) {
                for (int col = 0; col < W; col++) {
                    rawData[row][col] = src->data[row][col];
                }
            }
        }
    } else {
        // No bayer pattern
        // TODO: Is there a flat field correction possible?
        if (!rawData) {
            rawData(3 * W, H);
        }

        if (riDark && W == riDark->get_width() && H == riDark->get_height()) {
            for (int row = 0; row < H; row++) {
                for (int col = 0; col < W; col++) {
                    int c  = FC(row, col);
                    int c4 = ( c == 1 && !(row & 1) ) ? 3 : c;
                    rawData[row][3 * col + 0] = max(src->data[row][3 * col + 0] + black[c4] - riDark->data[row][3 * col + 0], 0.0f);
                    rawData[row][3 * col + 1] = max(src->data[row][3 * col + 1] + black[c4] - riDark->data[row][3 * col + 1], 0.0f);
                    rawData[row][3 * col + 2] = max(src->data[row][3 * col + 2] + black[c4] - riDark->data[row][3 * col + 2], 0.0f);
                }
            }
        } else {
            for (int row = 0; row < H; row++) {
                for (int col = 0; col < W; col++) {
                    rawData[row][3 * col + 0] = src->data[row][3 * col + 0];
                    rawData[row][3 * col + 1] = src->data[row][3 * col + 1];
                    rawData[row][3 * col + 2] = src->data[row][3 * col + 2];
                }
            }
        }
    }
}

SSEFUNCTION void RawImageSource::cfaboxblur(RawImage *riFlatFile, float* cfablur, const int boxH, const int boxW)
{

    if(boxW == 0 && boxH == 0) { // nothing to blur
        memcpy(cfablur, riFlatFile->data[0], W * H * sizeof(float));
        return;
    }

    float *tmpBuffer = nullptr;
    float *cfatmp = nullptr;
    float *srcVertical = nullptr;


    if(boxH > 0 && boxW > 0) {
        // we need a temporary buffer if we have to blur both directions
        tmpBuffer = (float (*)) calloc (H * W, sizeof * tmpBuffer);
    }

    if(boxH == 0) {
        // if boxH == 0 we can skip the vertical blur and process the horizontal blur from riFlatFile to cfablur without using a temporary buffer
        cfatmp = cfablur;
    } else {
        cfatmp = tmpBuffer;
    }

    if(boxW == 0) {
        // if boxW == 0 we can skip the horizontal blur and process the vertical blur from riFlatFile to cfablur without using a temporary buffer
        srcVertical = riFlatFile->data[0];
    } else {
        srcVertical = cfatmp;
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

        if(boxW > 0) {
            //box blur cfa image; box size = BS
            //horizontal blur
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 0; row < H; row++) {
                int len = boxW / 2 + 1;
                cfatmp[row * W + 0] = riFlatFile->data[row][0] / len;
                cfatmp[row * W + 1] = riFlatFile->data[row][1] / len;

                for (int j = 2; j <= boxW; j += 2) {
                    cfatmp[row * W + 0] += riFlatFile->data[row][j] / len;
                    cfatmp[row * W + 1] += riFlatFile->data[row][j + 1] / len;
                }

                for (int col = 2; col <= boxW; col += 2) {
                    cfatmp[row * W + col] = (cfatmp[row * W + col - 2] * len + riFlatFile->data[row][boxW + col]) / (len + 1);
                    cfatmp[row * W + col + 1] = (cfatmp[row * W + col - 1] * len + riFlatFile->data[row][boxW + col + 1]) / (len + 1);
                    len ++;
                }

                for (int col = boxW + 2; col < W - boxW; col++) {
                    cfatmp[row * W + col] = cfatmp[row * W + col - 2] + (riFlatFile->data[row][boxW + col] - cfatmp[row * W + col - boxW - 2]) / len;
                }

                for (int col = W - boxW; col < W; col += 2) {
                    cfatmp[row * W + col] = (cfatmp[row * W + col - 2] * len - cfatmp[row * W + col - boxW - 2]) / (len - 1);

                    if (col + 1 < W) {
                        cfatmp[row * W + col + 1] = (cfatmp[row * W + col - 1] * len - cfatmp[row * W + col - boxW - 1]) / (len - 1);
                    }

                    len --;
                }
            }
        }

        if(boxH > 0) {
            //vertical blur
#ifdef __SSE2__
            vfloat  leninitv = F2V(boxH / 2 + 1);
            vfloat  onev = F2V( 1.0f );
            vfloat  temp1v, temp2v, temp3v, temp4v, lenv, lenp1v, lenm1v;
            int row;
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int col = 0; col < W - 7; col += 8) {
                lenv = leninitv;
                temp1v = LVFU(srcVertical[0 * W + col]) / lenv;
                temp2v = LVFU(srcVertical[1 * W + col]) / lenv;
                temp3v = LVFU(srcVertical[0 * W + col + 4]) / lenv;
                temp4v = LVFU(srcVertical[1 * W + col + 4]) / lenv;

                for (int i = 2; i < boxH + 2; i += 2) {
                    temp1v += LVFU(srcVertical[i * W + col]) / lenv;
                    temp2v += LVFU(srcVertical[(i + 1) * W + col]) / lenv;
                    temp3v += LVFU(srcVertical[i * W + col + 4]) / lenv;
                    temp4v += LVFU(srcVertical[(i + 1) * W + col + 4]) / lenv;
                }

                STVFU(cfablur[0 * W + col], temp1v);
                STVFU(cfablur[1 * W + col], temp2v);
                STVFU(cfablur[0 * W + col + 4], temp3v);
                STVFU(cfablur[1 * W + col + 4], temp4v);

                for (row = 2; row < boxH + 2; row += 2) {
                    lenp1v = lenv + onev;
                    temp1v = (temp1v * lenv + LVFU(srcVertical[(row + boxH) * W + col])) / lenp1v;
                    temp2v = (temp2v * lenv + LVFU(srcVertical[(row + boxH + 1) * W + col])) / lenp1v;
                    temp3v = (temp3v * lenv + LVFU(srcVertical[(row + boxH) * W + col + 4])) / lenp1v;
                    temp4v = (temp4v * lenv + LVFU(srcVertical[(row + boxH + 1) * W + col + 4])) / lenp1v;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[(row + 1)*W + col], temp2v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                    STVFU(cfablur[(row + 1)*W + col + 4], temp4v);
                    lenv = lenp1v;
                }

                for (; row < H - boxH - 1; row += 2) {
                    temp1v = temp1v + (LVFU(srcVertical[(row + boxH) * W + col]) - LVFU(srcVertical[(row - boxH - 2) * W + col])) / lenv;
                    temp2v = temp2v + (LVFU(srcVertical[(row + 1 + boxH) * W + col]) - LVFU(srcVertical[(row + 1 - boxH - 2) * W + col])) / lenv;
                    temp3v = temp3v + (LVFU(srcVertical[(row + boxH) * W + col + 4]) - LVFU(srcVertical[(row - boxH - 2) * W + col + 4])) / lenv;
                    temp4v = temp4v + (LVFU(srcVertical[(row + 1 + boxH) * W + col + 4]) - LVFU(srcVertical[(row + 1 - boxH - 2) * W + col + 4])) / lenv;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[(row + 1)*W + col], temp2v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                    STVFU(cfablur[(row + 1)*W + col + 4], temp4v);
                }

                for(; row < H - boxH; row++) {
                    temp1v = temp1v + (LVFU(srcVertical[(row + boxH) * W + col]) - LVFU(srcVertical[(row - boxH - 2) * W + col])) / lenv;
                    temp3v = temp3v + (LVFU(srcVertical[(row + boxH) * W + col + 4]) - LVFU(srcVertical[(row - boxH - 2) * W + col + 4])) / lenv;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                    vfloat swapv = temp1v;
                    temp1v = temp2v;
                    temp2v = swapv;
                    swapv = temp3v;
                    temp3v = temp4v;
                    temp4v = swapv;
                }

                for (; row < H - 1; row += 2) {
                    lenm1v = lenv - onev;
                    temp1v = (temp1v * lenv - LVFU(srcVertical[(row - boxH - 2) * W + col])) / lenm1v;
                    temp2v = (temp2v * lenv - LVFU(srcVertical[(row - boxH - 1) * W + col])) / lenm1v;
                    temp3v = (temp3v * lenv - LVFU(srcVertical[(row - boxH - 2) * W + col + 4])) / lenm1v;
                    temp4v = (temp4v * lenv - LVFU(srcVertical[(row - boxH - 1) * W + col + 4])) / lenm1v;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[(row + 1)*W + col], temp2v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                    STVFU(cfablur[(row + 1)*W + col + 4], temp4v);
                    lenv = lenm1v;
                }

                for(; row < H; row++) {
                    lenm1v = lenv - onev;
                    temp1v = (temp1v * lenv - LVFU(srcVertical[(row - boxH - 2) * W + col])) / lenm1v;
                    temp3v = (temp3v * lenv - LVFU(srcVertical[(row - boxH - 2) * W + col + 4])) / lenm1v;
                    STVFU(cfablur[(row)*W + col], temp1v);
                    STVFU(cfablur[(row)*W + col + 4], temp3v);
                }

            }

            #pragma omp single

            for (int col = W - (W % 8); col < W; col++) {
                int len = boxH / 2 + 1;
                cfablur[0 * W + col] = srcVertical[0 * W + col] / len;
                cfablur[1 * W + col] = srcVertical[1 * W + col] / len;

                for (int i = 2; i < boxH + 2; i += 2) {
                    cfablur[0 * W + col] += srcVertical[i * W + col] / len;
                    cfablur[1 * W + col] += srcVertical[(i + 1) * W + col] / len;
                }

                for (int row = 2; row < boxH + 2; row += 2) {
                    cfablur[row * W + col] = (cfablur[(row - 2) * W + col] * len + srcVertical[(row + boxH) * W + col]) / (len + 1);
                    cfablur[(row + 1)*W + col] = (cfablur[(row - 1) * W + col] * len + srcVertical[(row + boxH + 1) * W + col]) / (len + 1);
                    len ++;
                }

                for (int row = boxH + 2; row < H - boxH; row++) {
                    cfablur[row * W + col] = cfablur[(row - 2) * W + col] + (srcVertical[(row + boxH) * W + col] - srcVertical[(row - boxH - 2) * W + col]) / len;
                }

                for (int row = H - boxH; row < H; row += 2) {
                    cfablur[row * W + col] = (cfablur[(row - 2) * W + col] * len - srcVertical[(row - boxH - 2) * W + col]) / (len - 1);

                    if (row + 1 < H) {
                        cfablur[(row + 1)*W + col] = (cfablur[(row - 1) * W + col] * len - srcVertical[(row - boxH - 1) * W + col]) / (len - 1);
                    }

                    len --;
                }
            }

#else
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int col = 0; col < W; col++) {
                int len = boxH / 2 + 1;
                cfablur[0 * W + col] = srcVertical[0 * W + col] / len;
                cfablur[1 * W + col] = srcVertical[1 * W + col] / len;

                for (int i = 2; i < boxH + 2; i += 2) {
                    cfablur[0 * W + col] += srcVertical[i * W + col] / len;
                    cfablur[1 * W + col] += srcVertical[(i + 1) * W + col] / len;
                }

                for (int row = 2; row < boxH + 2; row += 2) {
                    cfablur[row * W + col] = (cfablur[(row - 2) * W + col] * len + srcVertical[(row + boxH) * W + col]) / (len + 1);
                    cfablur[(row + 1)*W + col] = (cfablur[(row - 1) * W + col] * len + srcVertical[(row + boxH + 1) * W + col]) / (len + 1);
                    len ++;
                }

                for (int row = boxH + 2; row < H - boxH; row++) {
                    cfablur[row * W + col] = cfablur[(row - 2) * W + col] + (srcVertical[(row + boxH) * W + col] - srcVertical[(row - boxH - 2) * W + col]) / len;
                }

                for (int row = H - boxH; row < H; row += 2) {
                    cfablur[row * W + col] = (cfablur[(row - 2) * W + col] * len - srcVertical[(row - boxH - 2) * W + col]) / (len - 1);

                    if (row + 1 < H) {
                        cfablur[(row + 1)*W + col] = (cfablur[(row - 1) * W + col] * len - srcVertical[(row - boxH - 1) * W + col]) / (len - 1);
                    }

                    len --;
                }
            }

#endif
        }
    }

    if(tmpBuffer) {
        free (tmpBuffer);
    }
}


// Scale original pixels into the range 0 65535 using black offsets and multipliers
void RawImageSource::scaleColors(int winx, int winy, int winw, int winh, const RAWParams &raw, array2D<float> &rawData)
{
    chmax[0] = chmax[1] = chmax[2] = chmax[3] = 0; //channel maxima
    float black_lev[4] = {0.f};//black level

    //adjust black level  (eg Canon)
    bool isMono = false;

    if (getSensorType() == ST_BAYER || getSensorType() == ST_FOVEON ) {

        black_lev[0] = raw.bayersensor.black1; //R
        black_lev[1] = raw.bayersensor.black0; //G1
        black_lev[2] = raw.bayersensor.black2; //B
        black_lev[3] = raw.bayersensor.black3; //G2

        isMono = RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::mono] == raw.bayersensor.method;
    } else if (getSensorType() == ST_FUJI_XTRANS) {

        black_lev[0] = raw.xtranssensor.blackred; //R
        black_lev[1] = raw.xtranssensor.blackgreen; //G1
        black_lev[2] = raw.xtranssensor.blackblue; //B
        black_lev[3] = raw.xtranssensor.blackgreen; //G2  (set, only used with a Bayer filter)

        isMono = RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::mono] == raw.xtranssensor.method;
    }

    for(int i = 0; i < 4 ; i++) {
        cblacksom[i] = max( c_black[i] + black_lev[i], 0.0f );    // adjust black level
    }

    initialGain = calculate_scale_mul(scale_mul, ref_pre_mul, c_white, cblacksom, isMono, ri->get_colors()); // recalculate scale colors with adjusted levels

    //fprintf(stderr, "recalc: %f [%f %f %f %f]\n", initialGain, scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);
    for(int i = 0; i < 4 ; i++) {
        clmax[i] = (c_white[i] - cblacksom[i]) * scale_mul[i];    // raw clip level
    }

    // this seems strange, but it works

    // scale image colors

    if( ri->getSensorType() == ST_BAYER) {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            float tmpchmax[3];
            tmpchmax[0] = tmpchmax[1] = tmpchmax[2] = 0.0f;
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int row = winy; row < winy + winh; row ++)
            {
                for (int col = winx; col < winx + winw; col++) {
                    float val = rawData[row][col];
                    int c  = FC(row, col);                        // three colors,  0=R, 1=G,  2=B
                    int c4 = ( c == 1 && !(row & 1) ) ? 3 : c;    // four  colors,  0=R, 1=G1, 2=B, 3=G2
                    val -= cblacksom[c4];
                    val *= scale_mul[c4];
                    rawData[row][col] = (val);
                    tmpchmax[c] = max(tmpchmax[c], val);
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                chmax[0] = max(tmpchmax[0], chmax[0]);
                chmax[1] = max(tmpchmax[1], chmax[1]);
                chmax[2] = max(tmpchmax[2], chmax[2]);
            }
        }
    } else if ( ri->get_colors() == 1 ) {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            float tmpchmax = 0.0f;
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int row = winy; row < winy + winh; row ++)
            {
                for (int col = winx; col < winx + winw; col++) {
                    float val = rawData[row][col];
                    val -= cblacksom[0];
                    val *= scale_mul[0];
                    rawData[row][col] = (val);
                    tmpchmax = max(tmpchmax, val);
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                chmax[0] = chmax[1] = chmax[2] = chmax[3] = max(tmpchmax, chmax[0]);
            }
        }
    } else if(ri->getSensorType() == ST_FUJI_XTRANS) {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            float tmpchmax[3];
            tmpchmax[0] = tmpchmax[1] = tmpchmax[2] = 0.0f;
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int row = winy; row < winy + winh; row ++)
            {
                for (int col = winx; col < winx + winw; col++) {
                    float val = rawData[row][col];
                    int c = ri->XTRANSFC(row, col);
                    val -= cblacksom[c];
                    val *= scale_mul[c];

                    rawData[row][col] = (val);
                    tmpchmax[c] = max(tmpchmax[c], val);
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                chmax[0] = max(tmpchmax[0], chmax[0]);
                chmax[1] = max(tmpchmax[1], chmax[1]);
                chmax[2] = max(tmpchmax[2], chmax[2]);
            }
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            float tmpchmax[3];
            tmpchmax[0] = tmpchmax[1] = tmpchmax[2] = 0.0f;
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int row = winy; row < winy + winh; row ++)
            {
                for (int col = winx; col < winx + winw; col++) {
                    for (int c = 0; c < 3; c++) {                 // three colors,  0=R, 1=G,  2=B
                        float val = rawData[row][3 * col + c];
                        val -= cblacksom[c];
                        val *= scale_mul[c];
                        rawData[row][3 * col + c] = (val);
                        tmpchmax[c] = max(tmpchmax[c], val);
                    }
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                chmax[0] = max(tmpchmax[0], chmax[0]);
                chmax[1] = max(tmpchmax[1], chmax[1]);
                chmax[2] = max(tmpchmax[2], chmax[2]);
            }
        }
        chmax[3] = chmax[1];
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int RawImageSource::defTransform (int tran)
{

    int deg = ri->get_rotateDegree();

    if ((tran & TR_ROT) == TR_R180) {
        deg += 180;
    } else if ((tran & TR_ROT) == TR_R90) {
        deg += 90;
    } else if ((tran & TR_ROT) == TR_R270) {
        deg += 270;
    }

    deg %= 360;

    int ret = 0;

    if (deg == 90) {
        ret |= TR_R90;
    } else if (deg == 180) {
        ret |= TR_R180;
    } else if (deg == 270) {
        ret |= TR_R270;
    }

    if (tran & TR_HFLIP) {
        ret |= TR_HFLIP;
    }

    if (tran & TR_VFLIP) {
        ret |= TR_VFLIP;
    }

    return ret;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Thread called part
void RawImageSource::processFalseColorCorrectionThread  (Imagefloat* im, array2D<float> &rbconv_Y, array2D<float> &rbconv_I, array2D<float> &rbconv_Q, array2D<float> &rbout_I, array2D<float> &rbout_Q, const int row_from, const int row_to)
{

    const int W = im->getWidth();
    constexpr float onebynine = 1.f / 9.f;

#ifdef __SSE2__
    vfloat buffer[12];
    vfloat* pre1 = &buffer[0];
    vfloat* pre2 = &buffer[3];
    vfloat* post1 = &buffer[6];
    vfloat* post2 = &buffer[9];
#else
    float buffer[12];
    float* pre1 = &buffer[0];
    float* pre2 = &buffer[3];
    float* post1 = &buffer[6];
    float* post2 = &buffer[9];
#endif

    int px = (row_from - 1) % 3, cx = row_from % 3, nx = 0;

    convert_row_to_YIQ (im->r(row_from - 1), im->g(row_from - 1), im->b(row_from - 1), rbconv_Y[px], rbconv_I[px], rbconv_Q[px], W);
    convert_row_to_YIQ (im->r(row_from), im->g(row_from), im->b(row_from), rbconv_Y[cx], rbconv_I[cx], rbconv_Q[cx], W);

    for (int j = 0; j < W; j++) {
        rbout_I[px][j] = rbconv_I[px][j];
        rbout_Q[px][j] = rbconv_Q[px][j];
    }

    for (int i = row_from; i < row_to; i++) {

        px = (i - 1) % 3;
        cx = i % 3;
        nx = (i + 1) % 3;

        convert_row_to_YIQ (im->r(i + 1), im->g(i + 1), im->b(i + 1), rbconv_Y[nx], rbconv_I[nx], rbconv_Q[nx], W);

#ifdef __SSE2__
        pre1[0] = _mm_setr_ps(rbconv_I[px][0], rbconv_Q[px][0], 0, 0) , pre1[1] = _mm_setr_ps(rbconv_I[cx][0], rbconv_Q[cx][0], 0, 0), pre1[2] = _mm_setr_ps(rbconv_I[nx][0], rbconv_Q[nx][0], 0, 0);
        pre2[0] = _mm_setr_ps(rbconv_I[px][1], rbconv_Q[px][1], 0, 0) , pre2[1] = _mm_setr_ps(rbconv_I[cx][1], rbconv_Q[cx][1], 0, 0), pre2[2] = _mm_setr_ps(rbconv_I[nx][1], rbconv_Q[nx][1], 0, 0);

        // fill first element in rbout_I and rbout_Q
        rbout_I[cx][0] = rbconv_I[cx][0];
        rbout_Q[cx][0] = rbconv_Q[cx][0];

        // median I channel
        for (int j = 1; j < W - 2; j += 2) {
            post1[0] = _mm_setr_ps(rbconv_I[px][j + 1], rbconv_Q[px][j + 1], 0, 0), post1[1] = _mm_setr_ps(rbconv_I[cx][j + 1], rbconv_Q[cx][j + 1], 0, 0), post1[2] = _mm_setr_ps(rbconv_I[nx][j + 1], rbconv_Q[nx][j + 1], 0, 0);
            const auto middle = middle4of6(pre2[0], pre2[1], pre2[2], post1[0], post1[1], post1[2]);
            vfloat medianval = median(pre1[0], pre1[1], pre1[2], middle[0], middle[1], middle[2], middle[3]);
            rbout_I[cx][j] = medianval[0];
            rbout_Q[cx][j] = medianval[1];
            post2[0] = _mm_setr_ps(rbconv_I[px][j + 2], rbconv_Q[px][j + 2], 0, 0), post2[1] = _mm_setr_ps(rbconv_I[cx][j + 2], rbconv_Q[cx][j + 2], 0, 0), post2[2] = _mm_setr_ps(rbconv_I[nx][j + 2], rbconv_Q[nx][j + 2], 0, 0);
            medianval = median(post2[0], post2[1], post2[2], middle[0], middle[1], middle[2], middle[3]);
            rbout_I[cx][j + 1] = medianval[0];
            rbout_Q[cx][j + 1] = medianval[1];
            std::swap(pre1, post1);
            std::swap(pre2, post2);
        }

        // fill last elements in rbout_I and rbout_Q
        rbout_I[cx][W - 1] = rbconv_I[cx][W - 1];
        rbout_I[cx][W - 2] = rbconv_I[cx][W - 2];
        rbout_Q[cx][W - 1] = rbconv_Q[cx][W - 1];
        rbout_Q[cx][W - 2] = rbconv_Q[cx][W - 2];

#else
        pre1[0] = rbconv_I[px][0], pre1[1] = rbconv_I[cx][0], pre1[2] = rbconv_I[nx][0];
        pre2[0] = rbconv_I[px][1], pre2[1] = rbconv_I[cx][1], pre2[2] = rbconv_I[nx][1];

        // fill first element in rbout_I
        rbout_I[cx][0] = rbconv_I[cx][0];

        // median I channel
        for (int j = 1; j < W - 2; j += 2) {
            post1[0] = rbconv_I[px][j + 1], post1[1] = rbconv_I[cx][j + 1], post1[2] = rbconv_I[nx][j + 1];
            const auto middle = middle4of6(pre2[0], pre2[1], pre2[2], post1[0], post1[1], post1[2]);
            rbout_I[cx][j] = median(pre1[0], pre1[1], pre1[2], middle[0], middle[1], middle[2], middle[3]);
            post2[0] = rbconv_I[px][j + 2], post2[1] = rbconv_I[cx][j + 2], post2[2] = rbconv_I[nx][j + 2];
            rbout_I[cx][j + 1] = median(post2[0], post2[1], post2[2], middle[0], middle[1], middle[2], middle[3]);
            std::swap(pre1, post1);
            std::swap(pre2, post2);
        }

        // fill last elements in rbout_I
        rbout_I[cx][W - 1] = rbconv_I[cx][W - 1];
        rbout_I[cx][W - 2] = rbconv_I[cx][W - 2];

        pre1[0] = rbconv_Q[px][0], pre1[1] = rbconv_Q[cx][0], pre1[2] = rbconv_Q[nx][0];
        pre2[0] = rbconv_Q[px][1], pre2[1] = rbconv_Q[cx][1], pre2[2] = rbconv_Q[nx][1];

        // fill first element in rbout_Q
        rbout_Q[cx][0] = rbconv_Q[cx][0];

        // median Q channel
        for (int j = 1; j < W - 2; j += 2) {
            post1[0] = rbconv_Q[px][j + 1], post1[1] = rbconv_Q[cx][j + 1], post1[2] = rbconv_Q[nx][j + 1];
            const auto middle = middle4of6(pre2[0], pre2[1], pre2[2], post1[0], post1[1], post1[2]);
            rbout_Q[cx][j] = median(pre1[0], pre1[1], pre1[2], middle[0], middle[1], middle[2], middle[3]);
            post2[0] = rbconv_Q[px][j + 2], post2[1] = rbconv_Q[cx][j + 2], post2[2] = rbconv_Q[nx][j + 2];
            rbout_Q[cx][j + 1] = median(post2[0], post2[1], post2[2], middle[0], middle[1], middle[2], middle[3]);
            std::swap(pre1, post1);
            std::swap(pre2, post2);
        }

        // fill last elements in rbout_Q
        rbout_Q[cx][W - 1] = rbconv_Q[cx][W - 1];
        rbout_Q[cx][W - 2] = rbconv_Q[cx][W - 2];
#endif

        // blur i-1th row
        if (i > row_from) {
            convert_to_RGB (im->r(i - 1, 0), im->g(i - 1, 0), im->b(i - 1, 0), rbconv_Y[px][0], rbout_I[px][0], rbout_Q[px][0]);

#ifdef _OPENMP
            #pragma omp simd
#endif

            for (int j = 1; j < W - 1; j++) {
                float I = (rbout_I[px][j - 1] + rbout_I[px][j] + rbout_I[px][j + 1] + rbout_I[cx][j - 1] + rbout_I[cx][j] + rbout_I[cx][j + 1] + rbout_I[nx][j - 1] + rbout_I[nx][j] + rbout_I[nx][j + 1]) * onebynine;
                float Q = (rbout_Q[px][j - 1] + rbout_Q[px][j] + rbout_Q[px][j + 1] + rbout_Q[cx][j - 1] + rbout_Q[cx][j] + rbout_Q[cx][j + 1] + rbout_Q[nx][j - 1] + rbout_Q[nx][j] + rbout_Q[nx][j + 1]) * onebynine;
                convert_to_RGB (im->r(i - 1, j), im->g(i - 1, j), im->b(i - 1, j), rbconv_Y[px][j], I, Q);
            }

            convert_to_RGB (im->r(i - 1, W - 1), im->g(i - 1, W - 1), im->b(i - 1, W - 1), rbconv_Y[px][W - 1], rbout_I[px][W - 1], rbout_Q[px][W - 1]);
        }
    }

    // blur last 3 row and finalize H-1th row
    convert_to_RGB (im->r(row_to - 1, 0), im->g(row_to - 1, 0), im->b(row_to - 1, 0), rbconv_Y[cx][0], rbout_I[cx][0], rbout_Q[cx][0]);
#ifdef _OPENMP
    #pragma omp simd
#endif

    for (int j = 1; j < W - 1; j++) {
        float I = (rbout_I[px][j - 1] + rbout_I[px][j] + rbout_I[px][j + 1] + rbout_I[cx][j - 1] + rbout_I[cx][j] + rbout_I[cx][j + 1] + rbconv_I[nx][j - 1] + rbconv_I[nx][j] + rbconv_I[nx][j + 1]) * onebynine;
        float Q = (rbout_Q[px][j - 1] + rbout_Q[px][j] + rbout_Q[px][j + 1] + rbout_Q[cx][j - 1] + rbout_Q[cx][j] + rbout_Q[cx][j + 1] + rbconv_Q[nx][j - 1] + rbconv_Q[nx][j] + rbconv_Q[nx][j + 1]) * onebynine;
        convert_to_RGB (im->r(row_to - 1, j), im->g(row_to - 1, j), im->b(row_to - 1, j), rbconv_Y[cx][j], I, Q);
    }

    convert_to_RGB (im->r(row_to - 1, W - 1), im->g(row_to - 1, W - 1), im->b(row_to - 1, W - 1), rbconv_Y[cx][W - 1], rbout_I[cx][W - 1], rbout_Q[cx][W - 1]);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// correction_YIQ_LQ
void RawImageSource::processFalseColorCorrection  (Imagefloat* im, const int steps)
{

    if (im->getHeight() < 4 || steps < 1) {
        return;
    }

#ifdef _OPENMP
    #pragma omp parallel
    {
        multi_array2D<float, 5> buffer (W, 3);
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int blk = (im->getHeight() - 2) / nthreads;

        for (int t = 0; t < steps; t++) {

            if (tid < nthreads - 1) {
                processFalseColorCorrectionThread (im, buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], 1 + tid * blk, 1 + (tid + 1)*blk);
            } else {
                processFalseColorCorrectionThread (im, buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], 1 + tid * blk, im->getHeight() - 1);
            }

            #pragma omp barrier
        }
    }
#else
    multi_array2D<float, 5> buffer (W, 3);

    for (int t = 0; t < steps; t++) {
        processFalseColorCorrectionThread (im, buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], 1 , im->getHeight() - 1);
    }

#endif
}

// Some camera input profiles need gamma preprocessing
// gamma is applied before the CMS, correct line fac=lineFac*rawPixel+LineSum after the CMS
void RawImageSource::getProfilePreprocParams(cmsHPROFILE in, float& gammaFac, float& lineFac, float& lineSum)
{
    gammaFac = 0;
    lineFac = 1;
    lineSum = 0;

    char copyright[256];
    copyright[0] = 0;

    if (cmsGetProfileInfoASCII(in, cmsInfoCopyright, cmsNoLanguage, cmsNoCountry, copyright, 256) > 0) {
        if (strstr(copyright, "Phase One") != nullptr) {
            gammaFac = 0.55556;    // 1.8
        } else if (strstr(copyright, "Nikon Corporation") != nullptr) {
            gammaFac = 0.5;
            lineFac = -0.4;
            lineSum = 1.35; // determined in reverse by measuring NX an RT developed colorchecker PNGs
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

static void
lab2ProphotoRgbD50(float L, float A, float B, float& r, float& g, float& b)
{
    float X;
    float Y;
    float Z;
#define CLIP01(a) ((a)>0?((a)<1?(a):1):0)
    {
        // convert from Lab to XYZ
        float x, y, z, fx, fy, fz;

        fy = (L + 16.0f) / 116.0f;
        fx = A / 500.0f + fy;
        fz = fy - B / 200.0f;

        if (fy > 24.0f / 116.0f) {
            y = fy * fy * fy;
        } else {
            y = (fy - 16.0f / 116.0f) / 7.787036979f;
        }

        if (fx > 24.0f / 116.0f) {
            x = fx * fx * fx;
        } else {
            x = (fx - 16.0 / 116.0) / 7.787036979f;
        }

        if (fz > 24.0f / 116.0f) {
            z = fz * fz * fz;
        } else {
            z = (fz - 16.0f / 116.0f) / 7.787036979f;
        }

        //0.9642, 1.0000, 0.8249 D50
        X = x * 0.9642;
        Y = y;
        Z = z * 0.8249;
    }
    r = prophoto_xyz[0][0] * X + prophoto_xyz[0][1] * Y + prophoto_xyz[0][2] * Z;
    g = prophoto_xyz[1][0] * X + prophoto_xyz[1][1] * Y + prophoto_xyz[1][2] * Z;
    b = prophoto_xyz[2][0] * X + prophoto_xyz[2][1] * Y + prophoto_xyz[2][2] * Z;
    r = CLIP01(r);
    g = CLIP01(g);
    b = CLIP01(b);
}

// Converts raw image including ICC input profile to working space - floating point version
void RawImageSource::colorSpaceConversion_ (Imagefloat* im, ColorManagementParams &cmp, const ColorTemp &wb, double pre_mul[3], cmsHPROFILE embedded, cmsHPROFILE camprofile, double camMatrix[3][3], const std::string &camName)
{

//    MyTime t1, t2, t3;
//    t1.set ();
    cmsHPROFILE in;
    DCPProfile *dcpProf;

    if (!findInputProfile(cmp.input, embedded, camName, &dcpProf, in)) {
        return;
    }

    if (dcpProf != nullptr) {
        // DCP processing
        const DCPProfile::Triple pre_mul_row = {
            pre_mul[0],
            pre_mul[1],
            pre_mul[2]
        };
        const DCPProfile::Matrix cam_matrix = {{
                {camMatrix[0][0], camMatrix[0][1], camMatrix[0][2]},
                {camMatrix[1][0], camMatrix[1][1], camMatrix[1][2]},
                {camMatrix[2][0], camMatrix[2][1], camMatrix[2][2]}
            }
        };
        dcpProf->apply(im, cmp.dcpIlluminant, cmp.working, wb, pre_mul_row, cam_matrix, cmp.applyHueSatMap);
        return;
    }

    if (in == nullptr) {
        // use default camprofile, supplied by dcraw
        // in this case we avoid using the slllllooooooowwww lcms

        // Calculate matrix for direct conversion raw>working space
        TMatrix work = ICCStore::getInstance()->workingSpaceInverseMatrix (cmp.working);
        double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++) {
                    mat[i][j] += work[i][k] * camMatrix[k][j];    // rgb_xyz * imatrices.xyz_cam
                }

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < im->getHeight(); i++)
            for (int j = 0; j < im->getWidth(); j++) {

                float newr = mat[0][0] * im->r(i, j) + mat[0][1] * im->g(i, j) + mat[0][2] * im->b(i, j);
                float newg = mat[1][0] * im->r(i, j) + mat[1][1] * im->g(i, j) + mat[1][2] * im->b(i, j);
                float newb = mat[2][0] * im->r(i, j) + mat[2][1] * im->g(i, j) + mat[2][2] * im->b(i, j);

                im->r(i, j) = newr;
                im->g(i, j) = newg;
                im->b(i, j) = newb;

            }
    } else {
        const bool working_space_is_prophoto = (cmp.working == "ProPhoto");

        // use supplied input profile

        /*
          The goal here is to in addition to user-made custom ICC profiles also support profiles
          supplied with other popular raw converters. As curves affect color rendering and
          different raw converters deal with them differently (and few if any is as flexible
          as RawTherapee) we cannot really expect to get the *exact* same color rendering here.
          However we try hard to make the best out of it.

          Third-party input profiles that contain a LUT (usually A2B0 tag) often needs some preprocessing,
          as ICC LUTs are not really designed for dealing with linear camera data. Generally one
          must apply some sort of curve to get efficient use of the LUTs. Unfortunately how you
          should preprocess is not standardized so there are almost as many ways as there are
          software makers, and for each one we have to reverse engineer to find out how it has
          been done. (The ICC files made for RT has linear LUTs)

          ICC profiles which only contain the <r,g,b>XYZ tags (ie only a color matrix) should
          (hopefully) not require any pre-processing.

          Some LUT ICC profiles apply a contrast curve and desaturate highlights (to give a "film-like"
          behavior. These will generally work with RawTherapee, but will not produce good results when
          you enable highlight recovery/reconstruction, as that data is added linearly on top of the
          original range. RawTherapee works best with linear ICC profiles.
        */

        enum camera_icc_type {
            CAMERA_ICC_TYPE_GENERIC, // Generic, no special pre-processing required, RTs own is this way
            CAMERA_ICC_TYPE_PHASE_ONE, // Capture One profiles
            CAMERA_ICC_TYPE_LEAF, // Leaf profiles, former Leaf Capture now in Capture One, made for Leaf digital backs
            CAMERA_ICC_TYPE_NIKON // Nikon NX profiles
        } camera_icc_type = CAMERA_ICC_TYPE_GENERIC;

        float leaf_prophoto_mat[3][3];
        {
            // identify ICC type
            char copyright[256] = "";
            char description[256] = "";

            cmsGetProfileInfoASCII(in, cmsInfoCopyright, cmsNoLanguage, cmsNoCountry, copyright, 256);
            cmsGetProfileInfoASCII(in, cmsInfoDescription, cmsNoLanguage, cmsNoCountry, description, 256);
            camera_icc_type = CAMERA_ICC_TYPE_GENERIC;

            // Note: order the identification with the most detailed matching first since the more general ones may also match the more detailed
            if ((strstr(copyright, "Leaf") != nullptr ||
                    strstr(copyright, "Phase One A/S") != nullptr ||
                    strstr(copyright, "Kodak") != nullptr ||
                    strstr(copyright, "Creo") != nullptr) &&
                    (strstr(description, "LF2 ") == description ||
                     strstr(description, "LF3 ") == description ||
                     strstr(description, "LeafLF2") == description ||
                     strstr(description, "LeafLF3") == description ||
                     strstr(description, "LeafLF4") == description ||
                     strstr(description, "MamiyaLF2") == description ||
                     strstr(description, "MamiyaLF3") == description)) {
                camera_icc_type = CAMERA_ICC_TYPE_LEAF;
            } else if (strstr(copyright, "Phase One A/S") != nullptr) {
                camera_icc_type = CAMERA_ICC_TYPE_PHASE_ONE;
            } else if (strstr(copyright, "Nikon Corporation") != nullptr) {
                camera_icc_type = CAMERA_ICC_TYPE_NIKON;
            }
        }

        // Initialize transform
        cmsHTRANSFORM hTransform;
        cmsHPROFILE prophoto = ICCStore::getInstance()->workingSpace("ProPhoto"); // We always use Prophoto to apply the ICC profile to minimize problems with clipping in LUT conversion.
        bool transform_via_pcs_lab = false;
        bool separate_pcs_lab_highlights = false;
        lcmsMutex->lock ();

        switch (camera_icc_type) {
            case CAMERA_ICC_TYPE_PHASE_ONE:
            case CAMERA_ICC_TYPE_LEAF: {
                // These profiles have a RGB to Lab cLUT, gives gamma 1.8 output, and expects a "film-like" curve on input
                transform_via_pcs_lab = true;
                separate_pcs_lab_highlights = true;
                // We transform to Lab because we can and that we avoid getting an unnecessary unmatched gamma conversion which we would need to revert.
                hTransform = cmsCreateTransform (in, TYPE_RGB_FLT, nullptr, TYPE_Lab_FLT, INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );

                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        leaf_prophoto_mat[i][j] = 0;

                        for (int k = 0; k < 3; k++) {
                            leaf_prophoto_mat[i][j] += prophoto_xyz[i][k] * camMatrix[k][j];
                        }
                    }
                }

                break;
            }

            case CAMERA_ICC_TYPE_NIKON:
            case CAMERA_ICC_TYPE_GENERIC:
            default:
                hTransform = cmsCreateTransform (in, TYPE_RGB_FLT, prophoto, TYPE_RGB_FLT, INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety
                break;
        }

        lcmsMutex->unlock ();

        if (hTransform == nullptr) {
            // Fallback: create transform from camera profile. Should not happen normally.
            lcmsMutex->lock ();
            hTransform = cmsCreateTransform (camprofile, TYPE_RGB_FLT, prophoto, TYPE_RGB_FLT, INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );
            lcmsMutex->unlock ();
        }

        TMatrix toxyz = {}, torgb = {};

        if (!working_space_is_prophoto) {
            toxyz = ICCStore::getInstance()->workingSpaceMatrix ("ProPhoto");
            torgb = ICCStore::getInstance()->workingSpaceInverseMatrix (cmp.working); //sRGB .. Adobe...Wide...
        }

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            AlignedBuffer<float> buffer(im->getWidth() * 3);
            AlignedBuffer<float> hl_buffer(im->getWidth() * 3);
            AlignedBuffer<float> hl_scale(im->getWidth());
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif

            for ( int h = 0; h < im->getHeight(); ++h ) {
                float *p = buffer.data, *pR = im->r(h), *pG = im->g(h), *pB = im->b(h);

                // Apply pre-processing
                for ( int w = 0; w < im->getWidth(); ++w ) {
                    float r = *(pR++);
                    float g = *(pG++);
                    float b = *(pB++);

                    // convert to 0-1 range as LCMS expects that
                    r /= 65535.0f;
                    g /= 65535.0f;
                    b /= 65535.0f;

                    float maxc = max(r, g, b);

                    if (maxc <= 1.0) {
                        hl_scale.data[w] = 1.0;
                    } else {
                        // highlight recovery extend the range past the clip point, which means we can get values larger than 1.0 here.
                        // LUT ICC profiles only work in the 0-1 range so we scale down to fit and restore after conversion.
                        hl_scale.data[w] = 1.0 / maxc;
                        r *= hl_scale.data[w];
                        g *= hl_scale.data[w];
                        b *= hl_scale.data[w];
                    }

                    switch (camera_icc_type) {
                        case CAMERA_ICC_TYPE_PHASE_ONE:
                            // Here we apply a curve similar to Capture One's "Film Standard" + gamma, the reason is that the LUTs embedded in the
                            // ICCs are designed to work on such input, and if you provide it with a different curve you don't get as good result.
                            // We will revert this curve after we've made the color transform. However when we revert the curve, we'll notice that
                            // highlight rendering suffers due to that the LUT transform don't expand well, therefore we do a less compressed
                            // conversion too and mix them, this gives us the highest quality and most flexible result.
                            hl_buffer.data[3 * w + 0] = pow_F(r, 1.0 / 1.8);
                            hl_buffer.data[3 * w + 1] = pow_F(g, 1.0 / 1.8);
                            hl_buffer.data[3 * w + 2] = pow_F(b, 1.0 / 1.8);
                            r = phaseOneIccCurveInv->getVal(r);
                            g = phaseOneIccCurveInv->getVal(g);
                            b = phaseOneIccCurveInv->getVal(b);
                            break;

                        case CAMERA_ICC_TYPE_LEAF: {
                            // Leaf profiles expect that the camera native RGB has been converted to Prophoto RGB
                            float newr = leaf_prophoto_mat[0][0] * r + leaf_prophoto_mat[0][1] * g + leaf_prophoto_mat[0][2] * b;
                            float newg = leaf_prophoto_mat[1][0] * r + leaf_prophoto_mat[1][1] * g + leaf_prophoto_mat[1][2] * b;
                            float newb = leaf_prophoto_mat[2][0] * r + leaf_prophoto_mat[2][1] * g + leaf_prophoto_mat[2][2] * b;
                            hl_buffer.data[3 * w + 0] = pow_F(newr, 1.0 / 1.8);
                            hl_buffer.data[3 * w + 1] = pow_F(newg, 1.0 / 1.8);
                            hl_buffer.data[3 * w + 2] = pow_F(newb, 1.0 / 1.8);
                            r = phaseOneIccCurveInv->getVal(newr);
                            g = phaseOneIccCurveInv->getVal(newg);
                            b = phaseOneIccCurveInv->getVal(newb);
                            break;
                        }

                        case CAMERA_ICC_TYPE_NIKON:
                            // gamma 0.5
                            r = sqrtf(r);
                            g = sqrtf(g);
                            b = sqrtf(b);
                            break;

                        case CAMERA_ICC_TYPE_GENERIC:
                        default:
                            // do nothing
                            break;
                    }

                    *(p++) = r;
                    *(p++) = g;
                    *(p++) = b;
                }

                // Run icc transform
                cmsDoTransform (hTransform, buffer.data, buffer.data, im->getWidth());

                if (separate_pcs_lab_highlights) {
                    cmsDoTransform (hTransform, hl_buffer.data, hl_buffer.data, im->getWidth());
                }

                // Apply post-processing
                p = buffer.data;
                pR = im->r(h);
                pG = im->g(h);
                pB = im->b(h);

                for ( int w = 0; w < im->getWidth(); ++w ) {

                    float r, g, b, hr, hg, hb;

                    if (transform_via_pcs_lab) {
                        float L = *(p++);
                        float A = *(p++);
                        float B = *(p++);
                        // profile connection space CIELAB should have D50 illuminant
                        lab2ProphotoRgbD50(L, A, B, r, g, b);

                        if (separate_pcs_lab_highlights) {
                            lab2ProphotoRgbD50(hl_buffer.data[3 * w + 0], hl_buffer.data[3 * w + 1], hl_buffer.data[3 * w + 2], hr, hg, hb);
                        }
                    } else {
                        r = *(p++);
                        g = *(p++);
                        b = *(p++);
                    }

                    // restore pre-processing and/or add post-processing for the various ICC types
                    switch (camera_icc_type) {
                        default:
                            break;

                        case CAMERA_ICC_TYPE_PHASE_ONE:
                        case CAMERA_ICC_TYPE_LEAF: {
                            // note the 1/1.8 gamma, it's the gamma that the profile has applied, which we must revert before we can revert the curve
                            r = phaseOneIccCurve->getVal(pow_F(r, 1.0 / 1.8));
                            g = phaseOneIccCurve->getVal(pow_F(g, 1.0 / 1.8));
                            b = phaseOneIccCurve->getVal(pow_F(b, 1.0 / 1.8));
                            const float mix = 0.25; // may seem a low number, but remember this is linear space, mixing starts 2 stops from clipping
                            const float maxc = max(r, g, b);

                            if (maxc > mix) {
                                float fac = (maxc - mix) / (1.0 - mix);
                                fac = sqrtf(sqrtf(fac)); // gamma 0.25 to mix in highlight render relatively quick
                                r = (1.0 - fac) * r + fac * hr;
                                g = (1.0 - fac) * g + fac * hg;
                                b = (1.0 - fac) * b + fac * hb;
                            }

                            break;
                        }

                        case CAMERA_ICC_TYPE_NIKON: {
                            const float lineFac = -0.4;
                            const float lineSum = 1.35;
                            r *= r * lineFac + lineSum;
                            g *= g * lineFac + lineSum;
                            b *= b * lineFac + lineSum;
                            break;
                        }
                    }

                    // restore highlight scaling if any
                    if (hl_scale.data[w] != 1.0) {
                        float fac = 1.0 / hl_scale.data[w];
                        r *= fac;
                        g *= fac;
                        b *= fac;
                    }

                    // If we don't have ProPhoto as chosen working profile, convert. This conversion is clipless, ie if we convert
                    // to a small space such as sRGB we may end up with negative values and values larger than max.
                    if (!working_space_is_prophoto) {
                        //convert from Prophoto to XYZ
                        float x = (toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b ) ;
                        float y = (toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b ) ;
                        float z = (toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b ) ;
                        //convert from XYZ to cmp.working  (sRGB...Adobe...Wide..)
                        r = ((torgb[0][0] * x + torgb[0][1] * y + torgb[0][2] * z)) ;
                        g = ((torgb[1][0] * x + torgb[1][1] * y + torgb[1][2] * z)) ;
                        b = ((torgb[2][0] * x + torgb[2][1] * y + torgb[2][2] * z)) ;
                    }

                    // return to the 0.0 - 65535.0 range (with possible negative and > max values present)
                    r *= 65535.0;
                    g *= 65535.0;
                    b *= 65535.0;

                    *(pR++) = r;
                    *(pG++) = g;
                    *(pB++) = b;
                }
            }
        } // End of parallelization
        cmsDeleteTransform(hTransform);
    }

//t3.set ();
//        printf ("ICM TIME: %d usec\n", t3.etime(t1));
}


// Determine RAW input and output profiles. Returns TRUE on success
bool RawImageSource::findInputProfile(Glib::ustring inProfile, cmsHPROFILE embedded, std::string camName, DCPProfile **dcpProf, cmsHPROFILE& in)
{
    in = nullptr; // cam will be taken on NULL
    *dcpProf = nullptr;

    if (inProfile == "(none)") {
        return false;
    }

    if (inProfile == "(embedded)" && embedded) {
        in = embedded;
    } else if (inProfile == "(cameraICC)") {
        // DCPs have higher quality, so use them first
        *dcpProf = DCPStore::getInstance()->getStdProfile(camName);

        if (*dcpProf == nullptr) {
            in = ICCStore::getInstance()->getStdProfile(camName);
        }
    } else if (inProfile != "(camera)" && inProfile != "") {
        Glib::ustring normalName = inProfile;

        if (!inProfile.compare (0, 5, "file:")) {
            normalName = inProfile.substr(5);
        }

        if (DCPStore::getInstance()->isValidDCPFileName(normalName)) {
            *dcpProf = DCPStore::getInstance()->getProfile(normalName);
        }

        if (*dcpProf == nullptr) {
            in = ICCStore::getInstance()->getProfile (inProfile);
        }
    }

    // "in" might be NULL because of "not found". That's ok, we take the cam profile then

    return true;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// derived from Dcraw "blend_highlights()"
//  very effective to reduce (or remove) the magenta, but with levels of grey !
void RawImageSource::HLRecovery_blend(float* rin, float* gin, float* bin, int width, float maxval, float* hlmax)
{
    const int ColorCount = 3;

    // Transform matrixes rgb>lab and back
    static const float trans[2][ColorCount][ColorCount] = {
        { { 1, 1, 1 }, { 1.7320508, -1.7320508, 0 }, { -1, -1, 2 } },
        { { 1, 1, 1 }, { 1, -1, 1 }, { 1, 1, -1 } }
    };
    static const float itrans[2][ColorCount][ColorCount] = {
        { { 1, 0.8660254, -0.5 }, { 1, -0.8660254, -0.5 }, { 1, 0, 1 } },
        { { 1, 1, 1 }, { 1, -1, 1 }, { 1, 1, -1 } }
    };

#define FOREACHCOLOR for (int c=0; c < ColorCount; c++)

    float minpt = min(hlmax[0], hlmax[1], hlmax[2]); //min of the raw clip points
    //float maxpt=max(hlmax[0],hlmax[1],hlmax[2]);//max of the raw clip points
    //float medpt=hlmax[0]+hlmax[1]+hlmax[2]-minpt-maxpt;//median of the raw clip points
    float maxave = (hlmax[0] + hlmax[1] + hlmax[2]) / 3; //ave of the raw clip points
    //some thresholds:
    const float clipthresh = 0.95;
    const float fixthresh = 0.5;
    const float satthresh = 0.5;

    float clip[3];
    FOREACHCOLOR clip[c] = min(maxave, hlmax[c]);

    // Determine the maximum level (clip) of all channels
    const float clippt = clipthresh * maxval;
    const float fixpt = fixthresh * minpt;
    const float desatpt = satthresh * maxave + (1 - satthresh) * maxval;

    for (int col = 0; col < width; col++) {
        float rgb[ColorCount], cam[2][ColorCount], lab[2][ColorCount], sum[2], chratio, lratio = 0;
        float L, C, H;

        // Copy input pixel to rgb so it's easier to access in loops
        rgb[0] = rin[col];
        rgb[1] = gin[col];
        rgb[2] = bin[col];

        // If no channel is clipped, do nothing on pixel
        int c;

        for (c = 0; c < ColorCount; c++) {
            if (rgb[c] > clippt) {
                break;
            }
        }

        if (c == ColorCount) {
            continue;
        }

        // Initialize cam with raw input [0] and potentially clipped input [1]
        FOREACHCOLOR {
            lratio += min(rgb[c], clip[c]);
            cam[0][c] = rgb[c];
            cam[1][c] = min(cam[0][c], maxval);
        }

        // Calculate the lightness correction ratio (chratio)
        for (int i = 0; i < 2; i++) {
            FOREACHCOLOR {
                lab[i][c] = 0;

                for (int j = 0; j < ColorCount; j++)
                {
                    lab[i][c] += trans[ColorCount - 3][c][j] * cam[i][j];
                }
            }

            sum[i] = 0;

            for (int c = 1; c < ColorCount; c++) {
                sum[i] += SQR(lab[i][c]);
            }
        }

        chratio = (sqrt(sum[1] / sum[0]));

        // Apply ratio to lightness in LCH space
        for (int c = 1; c < ColorCount; c++) {
            lab[0][c] *= chratio;
        }

        // Transform back from LCH to RGB
        FOREACHCOLOR {
            cam[0][c] = 0;

            for (int j = 0; j < ColorCount; j++)
            {
                cam[0][c] += itrans[ColorCount - 3][c][j] * lab[0][j];
            }
        }
        FOREACHCOLOR rgb[c] = cam[0][c] / ColorCount;

        // Copy converted pixel back
        if (rin[col] > fixpt) {
            float rfrac = SQR((min(clip[0], rin[col]) - fixpt) / (clip[0] - fixpt));
            rin[col] = min(maxave, rfrac * rgb[0] + (1 - rfrac) * rin[col]);
        }

        if (gin[col] > fixpt) {
            float gfrac = SQR((min(clip[1], gin[col]) - fixpt) / (clip[1] - fixpt));
            gin[col] = min(maxave, gfrac * rgb[1] + (1 - gfrac) * gin[col]);
        }

        if (bin[col] > fixpt) {
            float bfrac = SQR((min(clip[2], bin[col]) - fixpt) / (clip[2] - fixpt));
            bin[col] = min(maxave, bfrac * rgb[2] + (1 - bfrac) * bin[col]);
        }

        lratio /= (rin[col] + gin[col] + bin[col]);
        L = (rin[col] + gin[col] + bin[col]) / 3;
        C = lratio * 1.732050808 * (rin[col] - gin[col]);
        H = lratio * (2 * bin[col] - rin[col] - gin[col]);
        rin[col] = L - H / 6.0 + C / 3.464101615;
        gin[col] = L - H / 6.0 - C / 3.464101615;
        bin[col] = L + H / 3.0;

        if ((L = (rin[col] + gin[col] + bin[col]) / 3) > desatpt) {
            float Lfrac = max(0.0f, (maxave - L) / (maxave - desatpt));
            C = Lfrac * 1.732050808 * (rin[col] - gin[col]);
            H = Lfrac * (2 * bin[col] - rin[col] - gin[col]);
            rin[col] = L - H / 6.0 + C / 3.464101615;
            gin[col] = L - H / 6.0 - C / 3.464101615;
            bin[col] = L + H / 3.0;
        }
    }
}

void RawImageSource::HLRecovery_Luminance (float* rin, float* gin, float* bin, float* rout, float* gout, float* bout, int width, float maxval)
{

    for (int i = 0; i < width; i++) {
        float r = rin[i], g = gin[i], b = bin[i];

        if (r > maxval || g > maxval || b > maxval) {
            float ro = min(r, maxval);
            float go = min(g, maxval);
            float bo = min(b, maxval);
            double L = r + g + b;
            double C = 1.732050808 * (r - g);
            double H = 2 * b - r - g;
            double Co = 1.732050808 * (ro - go);
            double Ho = 2 * bo - ro - go;

            if (r != g && g != b) {
                double ratio = sqrt ((Co * Co + Ho * Ho) / (C * C + H * H));
                C *= ratio;
                H *= ratio;
            }

            float rr = L / 3.0 - H / 6.0 + C / 3.464101615;
            float gr = L / 3.0 - H / 6.0 - C / 3.464101615;
            float br = L / 3.0 + H / 3.0;
            rout[i] = rr;
            gout[i] = gr;
            bout[i] = br;
        } else {
            rout[i] = rin[i];
            gout[i] = gin[i];
            bout[i] = bin[i];
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::HLRecovery_CIELab (float* rin, float* gin, float* bin, float* rout, float* gout, float* bout,
                                        int width, float maxval, double xyz_cam[3][3], double cam_xyz[3][3])
{

    //static bool crTableReady = false;

    // lookup table for Lab conversion
    // perhaps should be centralized, universally defined so we don't keep remaking it???
    /*for (int ix=0; ix < 0x10000; ix++) {
            float rx = ix / 65535.0;
            fv[ix] = rx > 0.008856 ? exp(1.0/3 * log(rx)) : 7.787*rx + 16/116.0;
        }*/
    //crTableReady = true;


    for (int i = 0; i < width; i++) {
        float r = rin[i], g = gin[i], b = bin[i];

        if (r > maxval || g > maxval || b > maxval) {
            float ro = min(r, maxval);
            float go = min(g, maxval);
            float bo = min(b, maxval);
            float yy = xyz_cam[1][0] * r + xyz_cam[1][1] * g + xyz_cam[1][2] * b;
            float fy = (yy < 65535.0 ? Color::cachef[yy] / 327.68 : std::cbrt(yy / MAXVALD));
            // compute LCH decompostion of the clipped pixel (only color information, thus C and H will be used)
            float x = xyz_cam[0][0] * ro + xyz_cam[0][1] * go + xyz_cam[0][2] * bo;
            float y = xyz_cam[1][0] * ro + xyz_cam[1][1] * go + xyz_cam[1][2] * bo;
            float z = xyz_cam[2][0] * ro + xyz_cam[2][1] * go + xyz_cam[2][2] * bo;
            x = (x < 65535.0 ? Color::cachef[x] / 327.68 : std::cbrt(x / MAXVALD));
            y = (y < 65535.0 ? Color::cachef[y] / 327.68 : std::cbrt(y / MAXVALD));
            z = (z < 65535.0 ? Color::cachef[z] / 327.68 : std::cbrt(z / MAXVALD));
            // convert back to rgb
            double fz = fy - y + z;
            double fx = fy + x - y;

            double zr = Color::f2xyz(fz);
            double xr = Color::f2xyz(fx);

            x = xr * 65535.0 ;
            y = yy;
            z = zr * 65535.0 ;
            float rr = cam_xyz[0][0] * x + cam_xyz[0][1] * y + cam_xyz[0][2] * z;
            float gr = cam_xyz[1][0] * x + cam_xyz[1][1] * y + cam_xyz[1][2] * z;
            float br = cam_xyz[2][0] * x + cam_xyz[2][1] * y + cam_xyz[2][2] * z;
            rout[i] = (rr);
            gout[i] = (gr);
            bout[i] = (br);
        } else {
            rout[i] = (rin[i]);
            gout[i] = (gin[i]);
            bout[i] = (bin[i]);
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hlRecovery (const std::string &method, float* red, float* green, float* blue, int width, float* hlmax )
{

    if (method == "Luminance") {
        HLRecovery_Luminance (red, green, blue, red, green, blue, width, 65535.0);
    } else if (method == "CIELab blending") {
        HLRecovery_CIELab (red, green, blue, red, green, blue, width, 65535.0, imatrices.xyz_cam, imatrices.cam_xyz);
    }
    else if (method == "Blend") { // derived from Dcraw
        HLRecovery_blend(red, green, blue, width, 65535.0, hlmax);
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::getAutoExpHistogram (LUTu & histogram, int& histcompr)
{
//    BENCHFUN
    histcompr = 3;

    histogram(65536 >> histcompr);
    histogram.clear();
    const float refwb[3] = {static_cast<float>(refwb_red  / (1 << histcompr)), static_cast<float>(refwb_green / (1 << histcompr)), static_cast<float>(refwb_blue / (1 << histcompr))};

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        LUTu tmphistogram(histogram.getSize());
        tmphistogram.clear();
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16) nowait
#endif

        for (int i = border; i < H - border; i++) {
            int start, end;
            getRowStartEnd (i, start, end);

            if (ri->getSensorType() == ST_BAYER) {
                // precalculate factors to avoid expensive per pixel calculations
                float refwb0 =  refwb[ri->FC(i, start)];
                float refwb1 =  refwb[ri->FC(i, start + 1)];
                int j;

                for (j = start; j < end - 1; j += 2) {
                    tmphistogram[(int)(refwb0 * rawData[i][j])] += 4;
                    tmphistogram[(int)(refwb1 * rawData[i][j + 1])] += 4;
                }

                if(j < end) {
                    tmphistogram[(int)(refwb0 * rawData[i][j])] += 4;
                }
            } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
                // precalculate factors to avoid expensive per pixel calculations
                float refwb0 =  refwb[ri->XTRANSFC(i, start)];
                float refwb1 =  refwb[ri->XTRANSFC(i, start + 1)];
                float refwb2 =  refwb[ri->XTRANSFC(i, start + 2)];
                float refwb3 =  refwb[ri->XTRANSFC(i, start + 3)];
                float refwb4 =  refwb[ri->XTRANSFC(i, start + 4)];
                float refwb5 =  refwb[ri->XTRANSFC(i, start + 5)];
                int j;

                for (j = start; j < end - 5; j += 6) {
                    tmphistogram[(int)(refwb0 * rawData[i][j])] += 4;
                    tmphistogram[(int)(refwb1 * rawData[i][j + 1])] += 4;
                    tmphistogram[(int)(refwb2 * rawData[i][j + 2])] += 4;
                    tmphistogram[(int)(refwb3 * rawData[i][j + 3])] += 4;
                    tmphistogram[(int)(refwb4 * rawData[i][j + 4])] += 4;
                    tmphistogram[(int)(refwb5 * rawData[i][j + 5])] += 4;
                }

                for (; j < end; j++) {
                    tmphistogram[(int)(refwb[ri->XTRANSFC(i, j)] * rawData[i][j])] += 4;
                }
            } else if (ri->get_colors() == 1) {
                for (int j = start; j < end; j++) {
                    tmphistogram[(int)(refwb[0] * rawData[i][j])]++;
                }
            } else {
                for (int j = start; j < end; j++) {
                    tmphistogram[(int)(refwb[0] * rawData[i][3 * j + 0])]++;
                    tmphistogram[(int)(refwb[1] * rawData[i][3 * j + 1])]++;
                    tmphistogram[(int)(refwb[2] * rawData[i][3 * j + 2])]++;
                }
            }
        }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            histogram += tmphistogram;
        }
    }
}

// Histogram MUST be 256 in size; gamma is applied, blackpoint and gain also
void RawImageSource::getRAWHistogram (LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw)
{
//    BENCHFUN
    histRedRaw.clear();
    histGreenRaw.clear();
    histBlueRaw.clear();
    const float mult[4] = { 65535.0f / ri->get_white(0),
                            65535.0f / ri->get_white(1),
                            65535.0f / ri->get_white(2),
                            65535.0f / ri->get_white(3)
                          };

    const bool fourColours = ri->getSensorType() == ST_BAYER && ((mult[1] != mult[3] || cblacksom[1] != cblacksom[3]) || FC(0, 0) == 3 || FC(0, 1) == 3 || FC(1, 0) == 3 || FC(1, 1) == 3);

    constexpr int histoSize = 65536;
    LUTu hist[4];
    hist[0](histoSize);
    hist[0].clear();

    if (ri->get_colors() > 1) {
        hist[1](histoSize);
        hist[1].clear();
        hist[2](histoSize);
        hist[2].clear();
    }

    if (fourColours) {
        hist[3](histoSize);
        hist[3].clear();
    }

#ifdef _OPENMP
    int numThreads;
    // reduce the number of threads under certain conditions to avoid overhead of too many critical regions
    numThreads = sqrt((((H - 2 * border) * (W - 2 * border)) / 262144.f));
    numThreads = std::min(std::max(numThreads, 1), omp_get_max_threads());

    #pragma omp parallel num_threads(numThreads)
#endif
    {
        // we need one LUT per color and thread, which corresponds to 1 MB per thread
        LUTu tmphist[4];
        tmphist[0](histoSize);
        tmphist[0].clear();

        if (ri->get_colors() > 1) {
            tmphist[1](histoSize);
            tmphist[1].clear();
            tmphist[2](histoSize);
            tmphist[2].clear();

            if (fourColours) {
                tmphist[3](histoSize);
                tmphist[3].clear();
            }
        }

#ifdef _OPENMP
        #pragma omp for nowait
#endif

        for (int i = border; i < H - border; i++) {
            int start, end;
            getRowStartEnd (i, start, end);

            if (ri->getSensorType() == ST_BAYER) {
                int j;
                int c1 = FC(i, start);
                c1 = ( fourColours && c1 == 1 && !(i & 1) ) ? 3 : c1;
                int c2 = FC(i, start + 1);
                c2 = ( fourColours && c2 == 1 && !(i & 1) ) ? 3 : c2;

                for (j = start; j < end - 1; j += 2) {
                    tmphist[c1][(int)ri->data[i][j]]++;
                    tmphist[c2][(int)ri->data[i][j + 1]]++;
                }

                if(j < end) { // last pixel of row if width is odd
                    tmphist[c1][(int)ri->data[i][j]]++;
                }
            } else if (ri->get_colors() == 1) {
                for (int j = start; j < end; j++) {
                    tmphist[0][(int)ri->data[i][j]]++;
                }
            } else if(ri->getSensorType() == ST_FUJI_XTRANS) {
                for (int j = start; j < end - 1; j += 2) {
                    int c = ri->XTRANSFC(i, j);
                    tmphist[c][(int)ri->data[i][j]]++;
                }
            } else {
                for (int j = start; j < end; j++) {
                    for (int c = 0; c < 3; c++) {
                        tmphist[c][(int)ri->data[i][3 * j + c]]++;
                    }
                }
            }
        }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            hist[0] += tmphist[0];

            if (ri->get_colors() > 1) {
                hist[1] += tmphist[1];
                hist[2] += tmphist[2];

                if (fourColours) {
                    hist[3] += tmphist[3];
                }
            }
        } // end of critical region
    } // end of parallel region

    for(int i = 0; i < 65536; i++) {
        int idx;
        idx = CLIP((int)Color::gamma(mult[0] * (i - (cblacksom[0]/*+black_lev[0]*/))));
        histRedRaw[idx >> 8] += hist[0][i];

        if (ri->get_colors() > 1) {
            idx = CLIP((int)Color::gamma(mult[1] * (i - (cblacksom[1]/*+black_lev[1]*/))));
            histGreenRaw[idx >> 8] += hist[1][i];

            if (fourColours) {
                idx = CLIP((int)Color::gamma(mult[3] * (i - (cblacksom[3]/*+black_lev[3]*/))));
                histGreenRaw[idx >> 8] += hist[3][i];
            }

            idx = CLIP((int)Color::gamma(mult[2] * (i - (cblacksom[2]/*+black_lev[2]*/))));
            histBlueRaw[idx >> 8] += hist[2][i];
        }
    }

    if (ri->getSensorType() == ST_BAYER)    // since there are twice as many greens, correct for it
        for (int i = 0; i < 256; i++) {
            histGreenRaw[i] >>= 1;
        }
    else if(ri->getSensorType() == ST_FUJI_XTRANS)  // since Xtrans has 2.5 as many greens, correct for it
        for (int i = 0; i < 256; i++) {
            histGreenRaw[i] = (histGreenRaw[i] * 2) / 5;
        }
    else if(ri->get_colors() == 1) { // monochrome sensor => set all histograms equal
        histGreenRaw += histRedRaw;
        histBlueRaw += histRedRaw;
    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::getRowStartEnd (int x, int &start, int &end)
{
    if (fuji) {
        int fw = ri->get_FujiWidth();
        start = ABS(fw - x) + border;
        end = min(H + W - fw - x, fw + x) - border;
    } else {
        start = border;
        end = W - border;
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void RawImageSource::getAutoWBMultipliers (double &rm, double &gm, double &bm)
{
//    BENCHFUN
    constexpr double clipHigh = 64000.0;

    if (ri->get_colors() == 1) {
        rm = gm = bm = 1;
        return;
    }

    if (redAWBMul != -1.) {
        rm = redAWBMul;
        gm = greenAWBMul;
        bm = blueAWBMul;
        return;
    }

    if (!isWBProviderReady()) {
        rm = -1.0;
        gm = -1.0;
        bm = -1.0;
        return;
    }

    double avg_r = 0;
    double avg_g = 0;
    double avg_b = 0;
    int rn = 0, gn = 0, bn = 0;

    if (fuji) {
        for (int i = 32; i < H - 32; i++) {
            int fw = ri->get_FujiWidth();
            int start = ABS(fw - i) + 32;
            int end = min(H + W - fw - i, fw + i) - 32;

            for (int j = start; j < end; j++) {
                if (ri->getSensorType() != ST_BAYER) {
                    double dr = CLIP(initialGain * (rawData[i][3 * j]  ));
                    double dg = CLIP(initialGain * (rawData[i][3 * j + 1]));
                    double db = CLIP(initialGain * (rawData[i][3 * j + 2]));

                    if (dr > clipHigh || dg > clipHigh || db > clipHigh) {
                        continue;
                    }

                    avg_r += dr;
                    avg_g += dg;
                    avg_b += db;
                    rn = gn = ++bn;
                } else {
                    int c = FC( i, j);
                    double d = CLIP(initialGain * (rawData[i][j]));

                    if (d > clipHigh) {
                        continue;
                    }

                    // Let's test green first, because they are more numerous
                    if (c == 1) {
                        avg_g += d;
                        gn++;
                    } else if (c == 0) {
                        avg_r += d;
                        rn++;
                    } else { /*if (c==2)*/
                        avg_b += d;
                        bn++;
                    }
                }
            }
        }
    } else {
        if (ri->getSensorType() != ST_BAYER) {
            if(ri->getSensorType() == ST_FUJI_XTRANS) {
                const double compval = clipHigh / initialGain;
#ifdef _OPENMP
                #pragma omp parallel
#endif
                {
                    double avg_c[3] = {0.0};
                    int cn[3] = {0};
#ifdef _OPENMP
                    #pragma omp for schedule(dynamic,16) nowait
#endif

                    for (int i = 32; i < H - 32; i++) {
                        for (int j = 32; j < W - 32; j++) {
                            // each loop read 1 rgb triplet value
                            double d = rawData[i][j];

                            if (d > compval) {
                                continue;
                            }

                            int c = ri->XTRANSFC(i, j);
                            avg_c[c] += d;
                            cn[c]++;
                        }
                    }

#ifdef _OPENMP
                    #pragma omp critical
#endif
                    {
                        avg_r += avg_c[0];
                        avg_g += avg_c[1];
                        avg_b += avg_c[2];
                        rn += cn[0];
                        gn += cn[1];
                        bn += cn[2];
                    }
                }
                avg_r *= initialGain;
                avg_g *= initialGain;
                avg_b *= initialGain;
            } else {
                for (int i = 32; i < H - 32; i++)
                    for (int j = 32; j < W - 32; j++) {
                        // each loop read 1 rgb triplet value

                        double dr = CLIP(initialGain * (rawData[i][3 * j]  ));
                        double dg = CLIP(initialGain * (rawData[i][3 * j + 1]));
                        double db = CLIP(initialGain * (rawData[i][3 * j + 2]));

                        if (dr > clipHigh || dg > clipHigh || db > clipHigh) {
                            continue;
                        }

                        avg_r += dr;
                        rn++;
                        avg_g += dg;
                        avg_b += db;
                    }

                gn = rn;
                bn = rn;
            }
        } else {
            //determine GRBG coset; (ey,ex) is the offset of the R subarray
            int ey, ex;

            if (ri->ISGREEN(0, 0)) { //first pixel is G
                if (ri->ISRED(0, 1)) {
                    ey = 0;
                    ex = 1;
                } else {
                    ey = 1;
                    ex = 0;
                }
            } else {//first pixel is R or B
                if (ri->ISRED(0, 0)) {
                    ey = 0;
                    ex = 0;
                } else {
                    ey = 1;
                    ex = 1;
                }
            }

            const double compval = clipHigh / initialGain;
#ifdef _OPENMP
            #pragma omp parallel for reduction(+:avg_r,avg_g,avg_b,rn,gn,bn) schedule(dynamic,8)
#endif

            for (int i = 32; i < H - 32; i += 2)
                for (int j = 32; j < W - 32; j += 2) {
                    //average each Bayer quartet component individually if non-clipped
                    double d[2][2];
                    d[0][0] = rawData[i][j];
                    d[0][1] = rawData[i][j + 1];
                    d[1][0] = rawData[i + 1][j];
                    d[1][1] = rawData[i + 1][j + 1];

                    if (d[ey][ex] <= compval) {
                        avg_r += d[ey][ex];
                        rn++;
                    }

                    if (d[1 - ey][ex] <= compval) {
                        avg_g += d[1 - ey][ex];
                        gn++;
                    }

                    if (d[ey][1 - ex] <= compval) {
                        avg_g += d[ey][1 - ex];
                        gn++;
                    }

                    if (d[1 - ey][1 - ex] <= compval) {
                        avg_b += d[1 - ey][1 - ex];
                        bn++;
                    }
                }

            avg_r *= initialGain;
            avg_g *= initialGain;
            avg_b *= initialGain;

        }
    }

    if( settings->verbose ) {
        printf ("AVG: %g %g %g\n", avg_r / std::max(1, rn), avg_g / std::max(1, gn), avg_b / std::max(1, bn));
    }

    //    return ColorTemp (pow(avg_r/rn, 1.0/6.0)*img_r, pow(avg_g/gn, 1.0/6.0)*img_g, pow(avg_b/bn, 1.0/6.0)*img_b);

    double reds   = avg_r / std::max(1, rn) * refwb_red;
    double greens = avg_g / std::max(1, gn) * refwb_green;
    double blues  = avg_b / std::max(1, bn) * refwb_blue;

    redAWBMul   = rm = imatrices.rgb_cam[0][0] * reds + imatrices.rgb_cam[0][1] * greens + imatrices.rgb_cam[0][2] * blues;
    greenAWBMul = gm = imatrices.rgb_cam[1][0] * reds + imatrices.rgb_cam[1][1] * greens + imatrices.rgb_cam[1][2] * blues;
    blueAWBMul  = bm = imatrices.rgb_cam[2][0] * reds + imatrices.rgb_cam[2][1] * greens + imatrices.rgb_cam[2][2] * blues;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ColorTemp RawImageSource::getSpotWB (std::vector<Coord2D> &red, std::vector<Coord2D> &green, std::vector<Coord2D> &blue, int tran, double equal)
{

    int x;
    int y;
    double reds = 0, greens = 0, blues = 0;
    unsigned int rn = 0;

    if (ri->getSensorType() != ST_BAYER) {
        if(ri->getSensorType() == ST_FUJI_XTRANS) {
            int d[9][2] = {{0, 0}, { -1, -1}, { -1, 0}, { -1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

            for (size_t i = 0; i < red.size(); i++) {
                transformPosition (red[i].x, red[i].y, tran, x, y);
                double rloc, gloc, bloc;
                int rnbrs, gnbrs, bnbrs;
                rloc = gloc = bloc = rnbrs = gnbrs = bnbrs = 0;

                for (int k = 0; k < 9; k++) {
                    int xv = x + d[k][0];
                    int yv = y + d[k][1];

                    if(xv >= 0 && yv >= 0 && xv < W && yv < H) {
                        if (ri->ISXTRANSRED(yv, xv)) { //RED
                            rloc += (rawData[yv][xv]);
                            rnbrs++;
                            continue;
                        } else if (ri->ISXTRANSBLUE(yv, xv)) { //BLUE
                            bloc += (rawData[yv][xv]);
                            bnbrs++;
                            continue;
                        } else { // GREEN
                            gloc += (rawData[yv][xv]);
                            gnbrs++;
                            continue;
                        }
                    }
                }

                rloc /= rnbrs;
                gloc /= gnbrs;
                bloc /= bnbrs;

                if (rloc * initialGain < 64000. && gloc * initialGain < 64000. && bloc * initialGain < 64000.) {
                    reds += rloc;
                    greens += gloc;
                    blues += bloc;
                    rn++;
                }
            }

        } else {
            int xmin, xmax, ymin, ymax;
            int xr, xg, xb, yr, yg, yb;

            for (size_t i = 0; i < red.size(); i++) {
                transformPosition (red[i].x, red[i].y, tran, xr, yr);
                transformPosition (green[i].x, green[i].y, tran, xg, yg);
                transformPosition (blue[i].x, blue[i].y, tran, xb, yb);

                if (initialGain * (rawData[yr][3 * xr]  ) > 52500 ||
                        initialGain * (rawData[yg][3 * xg + 1]) > 52500 ||
                        initialGain * (rawData[yb][3 * xb + 2]) > 52500) {
                    continue;
                }

                xmin = min(xr, xg, xb);
                xmax = max(xr, xg, xb);
                ymin = min(yr, yg, yb);
                ymax = max(yr, yg, yb);

                if (xmin >= 0 && ymin >= 0 && xmax < W && ymax < H) {
                    reds    += (rawData[yr][3 * xr]  );
                    greens  += (rawData[yg][3 * xg + 1]);
                    blues   += (rawData[yb][3 * xb + 2]);
                    rn++;
                }
            }
        }

    } else {

        int d[9][2] = {{0, 0}, { -1, -1}, { -1, 0}, { -1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

        for (size_t i = 0; i < red.size(); i++) {
            transformPosition (red[i].x, red[i].y, tran, x, y);
            double rloc, gloc, bloc;
            int rnbrs, gnbrs, bnbrs;
            rloc = gloc = bloc = rnbrs = gnbrs = bnbrs = 0;

            for (int k = 0; k < 9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                int c = FC(yv, xv);

                if(xv >= 0 && yv >= 0 && xv < W && yv < H) {
                    if (c == 0) { //RED
                        rloc += (rawData[yv][xv]);
                        rnbrs++;
                        continue;
                    } else if (c == 2) { //BLUE
                        bloc += (rawData[yv][xv]);
                        bnbrs++;
                        continue;
                    } else { // GREEN
                        gloc += (rawData[yv][xv]);
                        gnbrs++;
                        continue;
                    }
                }
            }

            rloc /= std::max(1, rnbrs);
            gloc /= std::max(1, gnbrs);
            bloc /= std::max(1, bnbrs);

            if (rloc * initialGain < 64000. && gloc * initialGain < 64000. && bloc * initialGain < 64000.) {
                reds += rloc;
                greens += gloc;
                blues += bloc;
                rn++;
            }

            transformPosition (green[i].x, green[i].y, tran, x, y);//these are redundant now ??? if not, repeat for these blocks same as for red[]
            rloc = gloc = bloc = rnbrs = gnbrs = bnbrs = 0;

            for (int k = 0; k < 9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                int c = FC(yv, xv);

                if(xv >= 0 && yv >= 0 && xv < W && yv < H) {
                    if (c == 0) { //RED
                        rloc += (rawData[yv][xv]);
                        rnbrs++;
                        continue;
                    } else if (c == 2) { //BLUE
                        bloc += (rawData[yv][xv]);
                        bnbrs++;
                        continue;
                    } else { // GREEN
                        gloc += (rawData[yv][xv]);
                        gnbrs++;
                        continue;
                    }
                }
            }

            rloc /= std::max(rnbrs, 1);
            gloc /= std::max(gnbrs, 1);
            bloc /= std::max(bnbrs, 1);

            if (rloc * initialGain < 64000. && gloc * initialGain < 64000. && bloc * initialGain < 64000.) {
                reds += rloc;
                greens += gloc;
                blues += bloc;
                rn++;
            }

            transformPosition (blue[i].x, blue[i].y, tran, x, y);
            rloc = gloc = bloc = rnbrs = gnbrs = bnbrs = 0;

            for (int k = 0; k < 9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                int c = FC(yv, xv);

                if(xv >= 0 && yv >= 0 && xv < W && yv < H) {
                    if (c == 0) { //RED
                        rloc += (rawData[yv][xv]);
                        rnbrs++;
                        continue;
                    } else if (c == 2) { //BLUE
                        bloc += (rawData[yv][xv]);
                        bnbrs++;
                        continue;
                    } else { // GREEN
                        gloc += (rawData[yv][xv]);
                        gnbrs++;
                        continue;
                    }
                }
            }

            rloc /= std::max(rnbrs, 1);
            gloc /= std::max(gnbrs, 1);
            bloc /= std::max(bnbrs, 1);

            if (rloc * initialGain < 64000. && gloc * initialGain < 64000. && bloc * initialGain < 64000.) {
                reds += rloc;
                greens += gloc;
                blues += bloc;
                rn++;
            }
        }
    }

    if (2u * rn < red.size()) {
        return ColorTemp (equal);
    } else {
        reds = reds / std::max(1u, rn) * refwb_red;
        greens = greens / std::max(1u, rn) * refwb_green;
        blues = blues / std::max(1u, rn) * refwb_blue;

        double rm = imatrices.rgb_cam[0][0] * reds + imatrices.rgb_cam[0][1] * greens + imatrices.rgb_cam[0][2] * blues;
        double gm = imatrices.rgb_cam[1][0] * reds + imatrices.rgb_cam[1][1] * greens + imatrices.rgb_cam[1][2] * blues;
        double bm = imatrices.rgb_cam[2][0] * reds + imatrices.rgb_cam[2][1] * greens + imatrices.rgb_cam[2][2] * blues;

        return ColorTemp (rm, gm, bm, equal);
    }
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::transformPosition (int x, int y, int tran, int& ttx, int& tty)
{

    tran = defTransform (tran);

    x += border;
    y += border;

    if (d1x) {
        if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
            x /= 2;
        } else {
            y /= 2;
        }
    }

    int w = W, h = H;

    if (fuji) {
        w = ri->get_FujiWidth() * 2 + 1;
        h = (H - ri->get_FujiWidth()) * 2 + 1;
    }

    int sw = w, sh = h;

    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = h;
        sh = w;
    }

    int ppx = x, ppy = y;

    if (tran & TR_HFLIP) {
        ppx = sw - 1 - x ;
    }

    if (tran & TR_VFLIP) {
        ppy = sh - 1 - y;
    }

    int tx = ppx;
    int ty = ppy;

    if ((tran & TR_ROT) == TR_R180) {
        tx = w - 1 - ppx;
        ty = h - 1 - ppy;
    } else if ((tran & TR_ROT) == TR_R90) {
        tx = ppy;
        ty = h - 1 - ppx;
    } else if ((tran & TR_ROT) == TR_R270) {
        tx = w - 1 - ppy;
        ty = ppx;
    }

    if (fuji) {
        ttx = (tx + ty) / 2;
        tty = (ty - tx) / 2 + ri->get_FujiWidth();
    } else {
        ttx = tx;
        tty = ty;
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::inverse33 (const double (*rgb_cam)[3], double (*cam_rgb)[3])
{
    double nom = (rgb_cam[0][2] * rgb_cam[1][1] * rgb_cam[2][0] - rgb_cam[0][1] * rgb_cam[1][2] * rgb_cam[2][0] -
                  rgb_cam[0][2] * rgb_cam[1][0] * rgb_cam[2][1] + rgb_cam[0][0] * rgb_cam[1][2] * rgb_cam[2][1] +
                  rgb_cam[0][1] * rgb_cam[1][0] * rgb_cam[2][2] - rgb_cam[0][0] * rgb_cam[1][1] * rgb_cam[2][2] );
    cam_rgb[0][0] = (rgb_cam[1][2] * rgb_cam[2][1] - rgb_cam[1][1] * rgb_cam[2][2]) / nom;
    cam_rgb[0][1] = -(rgb_cam[0][2] * rgb_cam[2][1] - rgb_cam[0][1] * rgb_cam[2][2]) / nom;
    cam_rgb[0][2] = (rgb_cam[0][2] * rgb_cam[1][1] - rgb_cam[0][1] * rgb_cam[1][2]) / nom;
    cam_rgb[1][0] = -(rgb_cam[1][2] * rgb_cam[2][0] - rgb_cam[1][0] * rgb_cam[2][2]) / nom;
    cam_rgb[1][1] = (rgb_cam[0][2] * rgb_cam[2][0] - rgb_cam[0][0] * rgb_cam[2][2]) / nom;
    cam_rgb[1][2] = -(rgb_cam[0][2] * rgb_cam[1][0] - rgb_cam[0][0] * rgb_cam[1][2]) / nom;
    cam_rgb[2][0] = (rgb_cam[1][1] * rgb_cam[2][0] - rgb_cam[1][0] * rgb_cam[2][1]) / nom;
    cam_rgb[2][1] = -(rgb_cam[0][1] * rgb_cam[2][0] - rgb_cam[0][0] * rgb_cam[2][1]) / nom;
    cam_rgb[2][2] = (rgb_cam[0][1] * rgb_cam[1][0] - rgb_cam[0][0] * rgb_cam[1][1]) / nom;
}

DiagonalCurve* RawImageSource::phaseOneIccCurve;
DiagonalCurve* RawImageSource::phaseOneIccCurveInv;

void RawImageSource::init ()
{

    {
        // Initialize Phase One ICC curves

        /* This curve is derived from TIFFTAG_TRANSFERFUNCTION of a Capture One P25+ image with applied film curve,
           exported to TIFF with embedded camera ICC. It's assumed to be similar to most standard curves in
           Capture One. It's not necessary to be exactly the same, it's just to be close to a typical curve to
           give the Phase One ICC files a good working space. */
        const double phase_one_forward[] = {
            0.0000000000, 0.0000000000, 0.0152590219, 0.0029602502, 0.0305180438, 0.0058899825, 0.0457770657, 0.0087739376, 0.0610360876, 0.0115968566,
            0.0762951095, 0.0143587396, 0.0915541314, 0.0171969177, 0.1068131533, 0.0201876860, 0.1220721752, 0.0232852674, 0.1373311971, 0.0264744030,
            0.1525902190, 0.0297245747, 0.1678492409, 0.0330205234, 0.1831082628, 0.0363775082, 0.1983672847, 0.0397802701, 0.2136263066, 0.0432593271,
            0.2288853285, 0.0467841611, 0.2441443503, 0.0503700313, 0.2594033722, 0.0540474556, 0.2746623941, 0.0577859159, 0.2899214160, 0.0616159304,
            0.3051804379, 0.0655222400, 0.3204394598, 0.0695353628, 0.3356984817, 0.0736552987, 0.3509575036, 0.0778973068, 0.3662165255, 0.0822461280,
            0.3814755474, 0.0867170214, 0.3967345693, 0.0913252461, 0.4119935912, 0.0960860609, 0.4272526131, 0.1009994659, 0.4425116350, 0.1060654612,
            0.4577706569, 0.1113298238, 0.4730296788, 0.1167925536, 0.4882887007, 0.1224841688, 0.5035477226, 0.1284046693, 0.5188067445, 0.1345540551,
            0.5340657664, 0.1409781033, 0.5493247883, 0.1476615549, 0.5645838102, 0.1546501869, 0.5798428321, 0.1619287404, 0.5951018540, 0.1695277333,
            0.6103608759, 0.1774776837, 0.6256198978, 0.1858091096, 0.6408789197, 0.1945525292, 0.6561379416, 0.2037384604, 0.6713969635, 0.2134279393,
            0.6866559854, 0.2236667430, 0.7019150072, 0.2345159075, 0.7171740291, 0.2460517281, 0.7324330510, 0.2583047227, 0.7476920729, 0.2714122225,
            0.7629510948, 0.2854352636, 0.7782101167, 0.3004959182, 0.7934691386, 0.3167620356, 0.8087281605, 0.3343862058, 0.8239871824, 0.3535820554,
            0.8392462043, 0.3745937285, 0.8545052262, 0.3977111467, 0.8697642481, 0.4232547494, 0.8850232700, 0.4515754940, 0.9002822919, 0.4830701152,
            0.9155413138, 0.5190966659, 0.9308003357, 0.5615320058, 0.9460593576, 0.6136263066, 0.9613183795, 0.6807965209, 0.9765774014, 0.7717402914,
            0.9918364233, 0.9052109560, 1.0000000000, 1.0000000000
        };
        std::vector<double> cForwardPoints;
        cForwardPoints.push_back(double(DCT_Spline));  // The first value is the curve type
        std::vector<double> cInversePoints;
        cInversePoints.push_back(double(DCT_Spline));  // The first value is the curve type

        for (unsigned int i = 0; i < sizeof(phase_one_forward) / sizeof(phase_one_forward[0]); i += 2) {
            cForwardPoints.push_back(phase_one_forward[i + 0]);
            cForwardPoints.push_back(phase_one_forward[i + 1]);
            cInversePoints.push_back(phase_one_forward[i + 1]);
            cInversePoints.push_back(phase_one_forward[i + 0]);
        }

        phaseOneIccCurve = new DiagonalCurve(cForwardPoints, CURVES_MIN_POLY_POINTS);
        phaseOneIccCurveInv = new DiagonalCurve(cInversePoints, CURVES_MIN_POLY_POINTS);
    }
}

void RawImageSource::getRawValues(int x, int y, int rotate, int &R, int &G, int &B)
{
    int xnew = x + border;
    int ynew = y + border;
    rotate += ri->get_rotateDegree();
    rotate %= 360;
    if (rotate == 90) {
        std::swap(xnew,ynew);
        ynew = H - 1 - ynew;
    } else if (rotate == 180) {
        xnew = W - 1 - xnew;
        ynew = H - 1 - ynew;
    } else if (rotate == 270) {
        std::swap(xnew,ynew);
        ynew = H - 1 - ynew;
        xnew = W - 1 - xnew;
        ynew = H - 1 - ynew;
    }

    int c = ri->getSensorType() == ST_FUJI_XTRANS ? ri->XTRANSFC(ynew,xnew) : ri->FC(ynew,xnew);
    int val = round(rawData[ynew][xnew] / scale_mul[c]);
    if(c == 0) {
        R = val; G = 0; B = 0;
    } else if(c == 2) {
        R = 0; G = 0; B = val;
    } else {
        R = 0; G = val; B = 0;
    }
}

void RawImageSource::cleanup ()
{
    delete phaseOneIccCurve;
    delete phaseOneIccCurveInv;
}

} /* namespace */
