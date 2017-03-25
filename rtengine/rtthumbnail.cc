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
#include "rtthumbnail.h"
#include "../rtgui/options.h"
#include "image8.h"
#include <lcms2.h>
#include "curves.h"
#include <glibmm.h>
#include "improcfun.h"
#include "colortemp.h"
#include "mytime.h"
#include "utils.h"
#include "iccstore.h"
#include "iccmatrices.h"
#include "rawimagesource.h"
#include "stdimagesource.h"
#include <glib/gstdio.h>
#include "rawimage.h"
#include "jpeg.h"
#include "../rtgui/ppversion.h"
#include "improccoordinator.h"
#include <locale.h>
//#define BENCHMARK
#include "StopWatch.h"


namespace
{

    bool checkRawImageThumb(const rtengine::RawImage& raw_image)
    {
        if (!raw_image.is_supportedThumb()) {
            return false;
        }

        const std::size_t length =
            fdata(raw_image.get_thumbOffset(), raw_image.get_file())[1] != 0xD8 && raw_image.is_ppmThumb()
                ? raw_image.get_thumbWidth() * raw_image.get_thumbHeight() * (raw_image.get_thumbBPS() / 8) * 3
                : raw_image.get_thumbLength();

        return raw_image.get_thumbOffset() + length <= raw_image.get_file()->size;
    }

}

extern Options options;

namespace rtengine
{
using namespace procparams;

Thumbnail* Thumbnail::loadFromImage (const Glib::ustring& fname, int &w, int &h, int fixwh, double wbEq, bool inspectorMode)
{

    StdImageSource imgSrc;

    if (imgSrc.load(fname)) {
        return nullptr;
    }

    ImageIO* img = imgSrc.getImageIO();

    Thumbnail* tpp = new Thumbnail ();

    unsigned char* data;
    img->getEmbeddedProfileData (tpp->embProfileLength, data);

    if (data && tpp->embProfileLength) {
        tpp->embProfileData = new unsigned char [tpp->embProfileLength];
        memcpy (tpp->embProfileData, data, tpp->embProfileLength);
    }

    tpp->scaleForSave = 8192;
    tpp->defGain = 1.0;
    tpp->gammaCorrected = false;
    tpp->isRaw = 0;
    memset (tpp->colorMatrix, 0, sizeof(tpp->colorMatrix));
    tpp->colorMatrix[0][0] = 1.0;
    tpp->colorMatrix[1][1] = 1.0;
    tpp->colorMatrix[2][2] = 1.0;

    if (inspectorMode) {
        // Special case, meaning that we want a full sized thumbnail image (e.g. for the Inspector feature)
        w = img->getWidth();
        h = img->getHeight();
        tpp->scale = 1.;
    } else {
        if (fixwh == 1) {
            w = h * img->getWidth() / img->getHeight();
            tpp->scale = (double)img->getHeight() / h;
        } else {
            h = w * img->getHeight() / img->getWidth();
            tpp->scale = (double)img->getWidth() / w;
        }
    }

    // bilinear interpolation
    if (tpp->thumbImg) {
        delete tpp->thumbImg;
        tpp->thumbImg = nullptr;
    }

    if (inspectorMode) {
        // we want an Image8
        if (img->getType() == rtengine::sImage8) {
            // copy the image
            Image8 *srcImg = static_cast<Image8*>(img);
            Image8 *thImg = new Image8 (w, h);
            srcImg->copyData(thImg);
            tpp->thumbImg = thImg;
        } else {
            // copy the image with a conversion
            tpp->thumbImg = resizeTo<Image8>(w, h, TI_Bilinear, img);
        }
    } else {
        // we want the same image type than the source file
        tpp->thumbImg = resizeToSameType(w, h, TI_Bilinear, img);

        // histogram computation
        tpp->aeHistCompression = 3;
        tpp->aeHistogram(65536 >> tpp->aeHistCompression);

        double avg_r = 0;
        double avg_g = 0;
        double avg_b = 0;
        int n = 0;

        if (img->getType() == rtengine::sImage8) {
            Image8 *image = static_cast<Image8*>(img);
            image->computeHistogramAutoWB(avg_r, avg_g, avg_b, n, tpp->aeHistogram, tpp->aeHistCompression);
        } else if (img->getType() == sImage16) {
            Image16 *image = static_cast<Image16*>(img);
            image->computeHistogramAutoWB(avg_r, avg_g, avg_b, n, tpp->aeHistogram, tpp->aeHistCompression);
        } else if (img->getType() == sImagefloat) {
            Imagefloat *image = static_cast<Imagefloat*>(img);
            image->computeHistogramAutoWB(avg_r, avg_g, avg_b, n, tpp->aeHistogram, tpp->aeHistCompression);
        } else {
            printf("loadFromImage: Unsupported image type \"%s\"!\n", img->getType());
        }

        if (n > 0) {
            ColorTemp cTemp;

            tpp->redAWBMul   = avg_r / double(n);
            tpp->greenAWBMul = avg_g / double(n);
            tpp->blueAWBMul  = avg_b / double(n);
            tpp->wbEqual = wbEq;
            tpp->wbTempBias = 0.0;

            cTemp.mul2temp (tpp->redAWBMul, tpp->greenAWBMul, tpp->blueAWBMul, tpp->wbEqual, tpp->autoWBTemp, tpp->autoWBGreen);
        }

        tpp->init ();
    }

    return tpp;
}

Thumbnail* Thumbnail::loadQuickFromRaw (const Glib::ustring& fname, RawMetaDataLocation& rml, int &w, int &h, int fixwh, bool rotate, bool inspectorMode)
{
    RawImage *ri = new RawImage(fname);
    unsigned int imageNum = 0;
    int r = ri->loadRaw(false, imageNum, false);

    if( r ) {
        delete ri;
        return nullptr;
    }

    rml.exifBase = ri->get_exifBase();
    rml.ciffBase = ri->get_ciffBase();
    rml.ciffLength = ri->get_ciffLen();

    Image8* img = new Image8 ();
    // No sample format detection occurred earlier, so we set them here,
    // as they are mandatory for the setScanline method
    img->setSampleFormat(IIOSF_UNSIGNED_CHAR);
    img->setSampleArrangement(IIOSA_CHUNKY);

    int err = 1;

    // See if it is something we support
    if (checkRawImageThumb(*ri)) {
        const char* data((const char*)fdata(ri->get_thumbOffset(), ri->get_file()));

        if ( (unsigned char)data[1] == 0xd8 ) {
            err = img->loadJPEGFromMemory(data, ri->get_thumbLength());
        } else if (ri->is_ppmThumb()) {
            err = img->loadPPMFromMemory(data, ri->get_thumbWidth(), ri->get_thumbHeight(), ri->get_thumbSwap(), ri->get_thumbBPS());
        }
    }

    // did we succeed?
    if ( err ) {
        printf("Could not extract thumb from %s\n", fname.data());
        delete img;
        delete ri;
        return nullptr;
    }

    Thumbnail* tpp = new Thumbnail ();

    tpp->isRaw = 1;
    memset (tpp->colorMatrix, 0, sizeof(tpp->colorMatrix));
    tpp->colorMatrix[0][0] = 1.0;
    tpp->colorMatrix[1][1] = 1.0;
    tpp->colorMatrix[2][2] = 1.0;

    if (inspectorMode) {
        // Special case, meaning that we want a full sized thumbnail image (e.g. for the Inspector feature)
        w = img->getWidth();
        h = img->getHeight();
        tpp->scale = 1.;
    } else {
        if (fixwh == 1) {
            w = h * img->getWidth() / img->getHeight();
            tpp->scale = (double)img->getHeight() / h;
        } else {
            h = w * img->getHeight() / img->getWidth();
            tpp->scale = (double)img->getWidth() / w;
        }
    }

    if (tpp->thumbImg) {
        delete tpp->thumbImg;
        tpp->thumbImg = nullptr;
    }

    if (inspectorMode) {
        tpp->thumbImg = img;
    } else {
        tpp->thumbImg = resizeTo<Image8>(w, h, TI_Nearest, img);
        delete img;
    }

    if (rotate && ri->get_rotateDegree() > 0) {
        std::string fname = ri->get_filename();
        std::string suffix = fname.length() > 4 ? fname.substr(fname.length() - 3) : "";

        for (int i = 0; i < suffix.length(); i++) {
            suffix[i] = std::tolower(suffix[i]);
        }

        // Leaf .mos, Mamiya .mef and Phase One .iiq files have thumbnails already rotated.
        if (suffix != "mos" && suffix != "mef" && suffix != "iiq")  {
            tpp->thumbImg->rotate(ri->get_rotateDegree());
            // width/height may have changed after rotating
            w = tpp->thumbImg->getWidth();
            h = tpp->thumbImg->getHeight();
        }
    }

    if (!inspectorMode) {
        tpp->init ();
    }

    delete ri;

    return tpp;
}

#define FISRED(filter,row,col) \
    ((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==0 || !filter)
#define FISGREEN(filter,row,col) \
    ((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==1 || !filter)
#define FISBLUE(filter,row,col) \
    ((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==2 || !filter)

RawMetaDataLocation Thumbnail::loadMetaDataFromRaw (const Glib::ustring& fname)
{
    RawMetaDataLocation rml;
    rml.exifBase = -1;
    rml.ciffBase = -1;
    rml.ciffLength = -1;

    RawImage ri(fname);
    unsigned int imageNum = 0;

    int r = ri.loadRaw(false, imageNum);

    if( !r ) {
        rml.exifBase = ri.get_exifBase();
        rml.ciffBase = ri.get_ciffBase();
        rml.ciffLength = ri.get_ciffLen();
    }

    return rml;
}

Thumbnail* Thumbnail::loadFromRaw (const Glib::ustring& fname, RawMetaDataLocation& rml, int &w, int &h, int fixwh, double wbEq, bool rotate, int imageNum)
{
    RawImage *ri = new RawImage (fname);
    unsigned int tempImageNum = 0;

    int r = ri->loadRaw(1, tempImageNum, 0);

    if( r ) {
        delete ri;
        return nullptr;
    }

    int width = ri->get_width();
    int height = ri->get_height();
    rtengine::Thumbnail* tpp = new rtengine::Thumbnail;

    tpp->isRaw = true;
    tpp->embProfile = nullptr;
    tpp->embProfileData = nullptr;
    tpp->embProfileLength = ri->get_profileLen();

    if (ri->get_profileLen())
        tpp->embProfile = cmsOpenProfileFromMem(ri->get_profile(),
                                                ri->get_profileLen()); //\ TODO check if mutex is needed

    tpp->redMultiplier = ri->get_pre_mul(0);
    tpp->greenMultiplier = ri->get_pre_mul(1);
    tpp->blueMultiplier = ri->get_pre_mul(2);

    ri->scale_colors();
    ri->pre_interpolate();

    rml.exifBase = ri->get_exifBase();
    rml.ciffBase = ri->get_ciffBase();
    rml.ciffLength = ri->get_ciffLen();

    tpp->camwbRed = tpp->redMultiplier / ri->get_pre_mul(0);
    tpp->camwbGreen = tpp->greenMultiplier / ri->get_pre_mul(1);
    tpp->camwbBlue = tpp->blueMultiplier / ri->get_pre_mul(2);
    tpp->defGain = 1.0 / min(ri->get_pre_mul(0), ri->get_pre_mul(1), ri->get_pre_mul(2));
    tpp->gammaCorrected = true;

    unsigned filter = ri->get_filters();
    int firstgreen = 1;

    // locate first green location in the first row
    if(ri->getSensorType() == ST_BAYER)
        while (!FISGREEN(filter, 1, firstgreen) && firstgreen < 3) {
            firstgreen++;
        }

    int skip = 1;

    if (ri->get_FujiWidth() != 0) {
        if (fixwh == 1) { // fix height, scale width
            skip = ((ri->get_height() - ri->get_FujiWidth()) / sqrt(0.5) - firstgreen - 1) / h;
        } else {
            skip = (ri->get_FujiWidth() / sqrt(0.5) - firstgreen - 1) / w;
        }
    } else {
        if (fixwh == 1) { // fix height, scale width
            skip = (ri->get_height() - firstgreen - 1) / h;
        } else {
            skip = (ri->get_width() - firstgreen - 1) / w;
        }
    }

    if (skip % 2) {
        skip--;
    }

    if (skip < 1) {
        skip = 1;
    }

    int hskip = skip, vskip = skip;

    if (!ri->get_model().compare("D1X")) {
        hskip *= 2;
    }

    int rofs = 0;
    int tmpw = (width - 2) / hskip;
    int tmph = (height - 2) / vskip;

    DCraw::dcrawImage_t image = ri->get_image();

    Imagefloat* tmpImg = new Imagefloat(tmpw, tmph);

    if (ri->getSensorType() == ST_BAYER) {
        // demosaicing! (sort of)
        for (int row = 1, y = 0; row < height - 1 && y < tmph; row += vskip, y++) {
            rofs = row * width;

            for (int col = firstgreen, x = 0; col < width - 1 && x < tmpw; col += hskip, x++) {
                int ofs = rofs + col;
                int g = image[ofs][1];
                int r, b;

                if (FISRED(filter, row, col + 1)) {
                    r = (image[ofs + 1    ][0] + image[ofs - 1    ][0]) >> 1;
                    b = (image[ofs + width][2] + image[ofs - width][2]) >> 1;
                } else {
                    b = (image[ofs + 1    ][2] + image[ofs - 1    ][2]) >> 1;
                    r = (image[ofs + width][0] + image[ofs - width][0]) >> 1;
                }

                tmpImg->r(y, x) = r;
                tmpImg->g(y, x) = g;
                tmpImg->b(y, x) = b;
            }
        }
    } else if (ri->get_colors() == 1) {
        for (int row = 1, y = 0; row < height - 1 && y < tmph; row += vskip, y++) {
            rofs = row * width;

            for (int col = firstgreen, x = 0; col < width - 1 && x < tmpw; col
                    += hskip, x++) {
                int ofs = rofs + col;
                tmpImg->r(y, x) = tmpImg->g(y, x) = tmpImg->b(y, x) = image[ofs][0];
            }
        }
    } else {
        if(ri->getSensorType() == ST_FUJI_XTRANS) {
            for( int row = 1, y = 0; row < height - 1 && y < tmph; row += vskip, y++) {
                rofs = row * width;

                for( int col = 1, x = 0; col < width - 1 && x < tmpw; col += hskip, x++ ) {
                    int ofs = rofs + col;
                    float sum[3] = {};
                    int c;

                    for(int v = -1; v <= 1; v++) {
                        for(int h = -1; h <= 1; h++) {
                            c = ri->XTRANSFC(row + v, col + h);
                            sum[c] += image[ofs + v * width + h][c];
                        }
                    }

                    c = ri->XTRANSFC(row, col);

                    switch (c) {
                        case 0:
                            tmpImg->r(y, x) = image[ofs][0];
                            tmpImg->g(y, x) = sum[1] / 5.f;
                            tmpImg->b(y, x) = sum[2] / 3.f;
                            break;

                        case 1:
                            tmpImg->r(y, x) = sum[0] / 2.f;
                            tmpImg->g(y, x) = image[ofs][1];
                            tmpImg->b(y, x) = sum[2] / 2.f;
                            break;

                        case 2:
                            tmpImg->r(y, x) = sum[0] / 3.f;
                            tmpImg->g(y, x) = sum[1] / 5.f;
                            tmpImg->b(y, x) = image[ofs][2];
                            break;
                    }
                }
            }
        } else {
            int iwidth = ri->get_iwidth();
            int iheight = ri->get_iheight();
            int left_margin = ri->get_leftmargin();
            firstgreen += left_margin;
            int top_margin = ri->get_topmargin();
            if(ri->get_maker() == "Sigma" && ri->DNGVERSION()) { // Hack to prevent sigma dng files from crashing
                tmpw = (width - 2 - left_margin) / hskip;
                tmph = (height - 2 - top_margin) / vskip;
            }

            for (int row = 1 + top_margin, y = 0; row < iheight + top_margin  - 1 && y < tmph; row += vskip, y++) {
                rofs = row * iwidth;

                for (int col = firstgreen, x = 0; col < iwidth + left_margin - 1 && x < tmpw; col += hskip, x++) {
                    int ofs = rofs + col;
                    tmpImg->r(y, x) = image[ofs][0];
                    tmpImg->g(y, x) = image[ofs][1];
                    tmpImg->b(y, x) = image[ofs][2];
                }
            }
        }
    }

    if (ri->get_FujiWidth() != 0) {
        int fw = ri->get_FujiWidth() / hskip;
        double step = sqrt(0.5);
        int wide = fw / step;
        int high = (tmph - fw) / step;
        Imagefloat* fImg = new Imagefloat(wide, high);
        float r, c;

        for (int row = 0; row < high; row++)
            for (int col = 0; col < wide; col++) {
                unsigned ur = r = fw + (row - col) * step;
                unsigned uc = c = (row + col) * step;

                if (ur > tmph - 2 || uc > tmpw - 2) {
                    continue;
                }

                double fr = r - ur;
                double fc = c - uc;
                fImg->r(row, col) = (tmpImg->r(ur, uc) * (1 - fc) + tmpImg->r(ur, uc + 1) * fc) * (1 - fr)   +   (tmpImg->r(ur + 1, uc) * (1 - fc) + tmpImg->r(ur + 1, uc + 1) * fc) * fr;
                fImg->g(row, col) = (tmpImg->g(ur, uc) * (1 - fc) + tmpImg->g(ur, uc + 1) * fc) * (1 - fr)   +   (tmpImg->g(ur + 1, uc) * (1 - fc) + tmpImg->g(ur + 1, uc + 1) * fc) * fr;
                fImg->b(row, col) = (tmpImg->b(ur, uc) * (1 - fc) + tmpImg->b(ur, uc + 1) * fc) * (1 - fr)   +   (tmpImg->b(ur + 1, uc) * (1 - fc) + tmpImg->b(ur + 1, uc + 1) * fc) * fr;
            }

        delete tmpImg;
        tmpImg = fImg;
        tmpw = wide;
        tmph = high;
    }

    const bool rotate_90 =
        rotate
        && (
            ri->get_rotateDegree() == 90
            || ri->get_rotateDegree() == 270
        );

    if (rotate_90) {
        std::swap(tmpw, tmph);
    }

    if (fixwh == 1) { // fix height, scale width
        w = tmpw * h / tmph;
    } else {
        h = tmph * w / tmpw;
    }

    if (tpp->thumbImg) {
        delete tpp->thumbImg;
    }

    if (rotate_90) {
        tpp->thumbImg = resizeTo<Image16>(h, w, TI_Bilinear, tmpImg);
    } else {
        tpp->thumbImg = resizeTo<Image16>(w, h, TI_Bilinear, tmpImg);
    }
    delete tmpImg;


    if (ri->get_FujiWidth() != 0) {
        tpp->scale = (double) (height - ri->get_FujiWidth()) * 2.0 / (rotate_90 ? w : h);
    } else {
        tpp->scale = (double) height / (rotate_90 ? w : h);
    }

    // generate histogram for auto exposure
    tpp->aeHistCompression = 3;
    tpp->aeHistogram(65536 >> tpp->aeHistCompression);
    tpp->aeHistogram.clear();
    int radd = 4;
    int gadd = 4;
    int badd = 4;

    if (!filter) {
        radd = gadd = badd = 1;
    }

    for (int i = 8; i < height - 8; i++) {
        int start, end;

        if (ri->get_FujiWidth() != 0) {
            int fw = ri->get_FujiWidth();
            start = ABS(fw - i) + 8;
            end = min(height + width - fw - i, fw + i) - 8;
        } else {
            start = 8;
            end = width - 8;
        }

        if (ri->get_colors() == 1) {
            for (int j = start; j < end; j++) {
                tpp->aeHistogram[((int)(image[i * width + j][0])) >> tpp->aeHistCompression] += radd;
                tpp->aeHistogram[((int)(image[i * width + j][0])) >> tpp->aeHistCompression] += gadd;
                tpp->aeHistogram[((int)(image[i * width + j][0])) >> tpp->aeHistCompression] += badd;
            }
        } else if(ri->getSensorType() == ST_BAYER) {
            for (int j = start; j < end; j++)
                if (FISGREEN(filter, i, j)) {
                    tpp->aeHistogram[((int)(tpp->camwbGreen * image[i * width + j][1])) >> tpp->aeHistCompression] += gadd;
                } else if (FISRED(filter, i, j)) {
                    tpp->aeHistogram[((int)(tpp->camwbRed * image[i * width + j][0])) >> tpp->aeHistCompression] += radd;
                } else if (FISBLUE(filter, i, j)) {
                    tpp->aeHistogram[((int)(tpp->camwbBlue * image[i * width + j][2])) >> tpp->aeHistCompression] += badd;
                }
        } else if(ri->getSensorType() == ST_FUJI_XTRANS) {
            for (int j = start; j < end; j++)
                if (ri->ISXTRANSGREEN(i, j)) {
                    tpp->aeHistogram[((int)(tpp->camwbGreen * image[i * width + j][1])) >> tpp->aeHistCompression] += gadd;
                } else if (ri->ISXTRANSRED(i, j)) {
                    tpp->aeHistogram[((int)(tpp->camwbRed * image[i * width + j][0])) >> tpp->aeHistCompression] += radd;
                } else if (ri->ISXTRANSBLUE(i, j)) {
                    tpp->aeHistogram[((int)(tpp->camwbBlue * image[i * width + j][2])) >> tpp->aeHistCompression] += badd;
                }
        } else { /* if(ri->getSensorType()==ST_FOVEON) */
            for (int j = start; j < end; j++) {
                tpp->aeHistogram[((int)(image[i * width + j][0] * 2.f)) >> tpp->aeHistCompression] += radd;
                tpp->aeHistogram[((int)(image[i * width + j][1])) >> tpp->aeHistCompression] += gadd;
                tpp->aeHistogram[((int)(image[i * width + j][2] * 0.5f)) >> tpp->aeHistCompression] += badd;
            }
        }
    }

    // generate autoWB
    double avg_r = 0;
    double avg_g = 0;
    double avg_b = 0;
    const float eps = 1e-5; //tolerance to avoid dividing by zero

    float rn = eps, gn = eps, bn = eps;

    for (int i = 32; i < height - 32; i++) {
        int start, end;

        if (ri->get_FujiWidth() != 0) {
            int fw = ri->get_FujiWidth();
            start = ABS(fw - i) + 32;
            end = min(height + width - fw - i, fw + i) - 32;
        } else {
            start = 32;
            end = width - 32;
        }

        if(ri->getSensorType() == ST_BAYER) {
            for (int j = start; j < end; j++) {
                if (!filter) {
                    double d = tpp->defGain * image[i * width + j][0];

                    if (d > 64000.) {
                        continue;
                    }

                    avg_g += d;
                    avg_r += d;
                    avg_b += d;
                    rn++;
                    gn++;
                    bn++;
                } else if (FISGREEN(filter, i, j)) {
                    double d = tpp->defGain * image[i * width + j][1];

                    if (d > 64000.) {
                        continue;
                    }

                    avg_g += d;
                    gn++;
                } else if (FISRED(filter, i, j)) {
                    double d = tpp->defGain * image[i * width + j][0];

                    if (d > 64000.) {
                        continue;
                    }

                    avg_r += d;
                    rn++;
                } else if (FISBLUE(filter, i, j)) {
                    double d = tpp->defGain * image[i * width + j][2];

                    if (d > 64000.) {
                        continue;
                    }

                    avg_b += d;
                    bn++;
                }
            }
        } else if(ri->getSensorType() == ST_FUJI_XTRANS) {
            for (int j = start; j < end; j++) {
                if (ri->ISXTRANSGREEN(i, j)) {
                    double d = tpp->defGain * image[i * width + j][1];

                    if (d > 64000.) {
                        continue;
                    }

                    avg_g += d;
                    gn++;
                } else if (ri->ISXTRANSRED(i, j)) {
                    double d = tpp->defGain * image[i * width + j][0];

                    if (d > 64000.) {
                        continue;
                    }

                    avg_r += d;
                    rn++;
                } else if (ri->ISXTRANSBLUE(i, j)) {
                    double d = tpp->defGain * image[i * width + j][2];

                    if (d > 64000.) {
                        continue;
                    }

                    avg_b += d;
                    bn++;
                }
            }
        } else { /* if(ri->getSensorType()==ST_FOVEON) */
            for (int j = start; j < end; j++) {
                double d = tpp->defGain * image[i * width + j][0];

                if (d <= 64000.) {
                    avg_r += d;
                    rn++;
                }

                d = tpp->defGain * image[i * width + j][1];

                if (d <= 64000.) {
                    avg_g += d;
                    gn++;
                }

                d = tpp->defGain * image[i * width + j][2];

                if (d <= 64000.) {
                    avg_b += d;
                    bn++;
                }
            }
        }
    }

    double reds = avg_r / rn * tpp->camwbRed;
    double greens = avg_g / gn * tpp->camwbGreen;
    double blues = avg_b / bn * tpp->camwbBlue;

    tpp->redAWBMul   = ri->get_rgb_cam(0, 0) * reds + ri->get_rgb_cam(0, 1) * greens + ri->get_rgb_cam(0, 2) * blues;
    tpp->greenAWBMul = ri->get_rgb_cam(1, 0) * reds + ri->get_rgb_cam(1, 1) * greens + ri->get_rgb_cam(1, 2) * blues;
    tpp->blueAWBMul  = ri->get_rgb_cam(2, 0) * reds + ri->get_rgb_cam(2, 1) * greens + ri->get_rgb_cam(2, 2) * blues;
    tpp->wbEqual = wbEq;
    tpp->wbTempBias = 0.0;

    ColorTemp cTemp;
    cTemp.mul2temp(tpp->redAWBMul, tpp->greenAWBMul, tpp->blueAWBMul, tpp->wbEqual, tpp->autoWBTemp, tpp->autoWBGreen);

    if (rotate && ri->get_rotateDegree() > 0) {
        tpp->thumbImg->rotate(ri->get_rotateDegree());
    }

    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            tpp->colorMatrix[a][b] = ri->get_rgb_cam(a, b);
        }

    tpp->init();
    delete ri;
    return tpp;
}
#undef FISRED
#undef FISGREEN
#undef FISBLUE

void Thumbnail::init ()
{
    RawImageSource::inverse33 (colorMatrix, iColorMatrix);
    //colorMatrix is rgb_cam
    memset (cam2xyz, 0, sizeof(cam2xyz));

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) {
                cam2xyz[i][j] += xyz_sRGB[i][k] * colorMatrix[k][j];
            }

    camProfile = ICCStore::getInstance()->createFromMatrix (cam2xyz, false, "Camera");
}

Thumbnail::Thumbnail () :
    iColorMatrix{}, cam2xyz{}, scale(1.0), colorMatrix{}, isRaw(true),
    camProfile(nullptr), thumbImg(nullptr),
    camwbRed(1.0), camwbGreen(1.0), camwbBlue(1.0),
    redAWBMul(-1.0), greenAWBMul(-1.0), blueAWBMul(-1.0),
    autoWBTemp(2700), autoWBGreen(1.0), wbEqual(-1.0), wbTempBias(0.0),
    embProfileLength(0), embProfileData(nullptr), embProfile(nullptr),
    redMultiplier(1.0), greenMultiplier(1.0), blueMultiplier(1.0),
    defGain(1.0),
    scaleForSave(8192),
    gammaCorrected(false),
    aeHistCompression(3)
{
}

Thumbnail::~Thumbnail ()
{

    delete thumbImg;
    //delete [] aeHistogram;
    delete [] embProfileData;

    if (embProfile) {
        cmsCloseProfile(embProfile);
    }

    if (camProfile) {
        cmsCloseProfile(camProfile);
    }
}

// Simple processing of RAW internal JPGs
IImage8* Thumbnail::quickProcessImage (const procparams::ProcParams& params, int rheight, rtengine::TypeInterpolation interp, double& myscale)
{

    int rwidth;

    if (params.coarse.rotate == 90 || params.coarse.rotate == 270) {
        rwidth = rheight;
        rheight = thumbImg->getHeight() * rwidth / thumbImg->getWidth();
    } else {
        rwidth = thumbImg->getWidth() * rheight / thumbImg->getHeight();
    }

    Image8* baseImg = resizeTo<Image8>(rwidth, rheight, interp, thumbImg);

    if (params.coarse.rotate) {
        baseImg->rotate (params.coarse.rotate);
    }

    if (params.coarse.hflip) {
        baseImg->hflip ();
    }

    if (params.coarse.vflip) {
        baseImg->vflip ();
    }

    return baseImg;
}

// Full thumbnail processing, second stage if complete profile exists
IImage8* Thumbnail::processImage (const procparams::ProcParams& params, int rheight, TypeInterpolation interp, std::string camName,
                                  double focalLen, double focalLen35mm, float focusDist, float shutter, float fnumber, float iso, std::string expcomp_, double& myscale)
{
    BENCHFUN

    // check if the WB's equalizer value has changed
    if (wbEqual < (params.wb.equal - 5e-4) || wbEqual > (params.wb.equal + 5e-4) || wbTempBias < (params.wb.tempBias - 5e-4) || wbTempBias > (params.wb.tempBias + 5e-4)) {
        wbEqual = params.wb.equal;
        wbTempBias = params.wb.tempBias;
        // recompute the autoWB
        ColorTemp cTemp;
        cTemp.mul2temp (redAWBMul, greenAWBMul, blueAWBMul, wbEqual, autoWBTemp, autoWBGreen);
        autoWBTemp += autoWBTemp * wbTempBias;
    }

    // compute WB multipliers
    ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

    if (params.wb.method == "Camera") {
        //recall colorMatrix is rgb_cam
        double cam_r = colorMatrix[0][0] * camwbRed + colorMatrix[0][1] * camwbGreen + colorMatrix[0][2] * camwbBlue;
        double cam_g = colorMatrix[1][0] * camwbRed + colorMatrix[1][1] * camwbGreen + colorMatrix[1][2] * camwbBlue;
        double cam_b = colorMatrix[2][0] * camwbRed + colorMatrix[2][1] * camwbGreen + colorMatrix[2][2] * camwbBlue;
        currWB = ColorTemp (cam_r, cam_g, cam_b, params.wb.equal);
    } else if (params.wb.method == "Auto") {
        currWB = ColorTemp (autoWBTemp, autoWBGreen, wbEqual, "Custom");
    }

    double r, g, b;
    currWB.getMultipliers (r, g, b);
    //iColorMatrix is cam_rgb
    double rm = iColorMatrix[0][0] * r + iColorMatrix[0][1] * g + iColorMatrix[0][2] * b;
    double gm = iColorMatrix[1][0] * r + iColorMatrix[1][1] * g + iColorMatrix[1][2] * b;
    double bm = iColorMatrix[2][0] * r + iColorMatrix[2][1] * g + iColorMatrix[2][2] * b;
    rm = camwbRed / rm;
    gm = camwbGreen / gm;
    bm = camwbBlue / bm;
    double mul_lum = 0.299 * rm + 0.587 * gm + 0.114 * bm;
    double logDefGain = 0.0;
    float rmi, gmi, bmi;

    rmi = rm * defGain / mul_lum;
    gmi = gm * defGain / mul_lum;
    bmi = bm * defGain / mul_lum;

    // The RAW exposure is not reflected since it's done in preprocessing. If we only have e.g. the chached thumb,
    // that is already preprocessed. So we simulate the effect here roughly my modifying the exposure accordingly
    if (isRaw) {
        rmi *= params.raw.expos;
        gmi *= params.raw.expos;
        bmi *= params.raw.expos;
    }

    // resize to requested width and perform coarse transformation
    int rwidth;

    if (params.coarse.rotate == 90 || params.coarse.rotate == 270) {
        rwidth = rheight;
        rheight = int(size_t(thumbImg->getHeight()) * size_t(rwidth) / size_t(thumbImg->getWidth()));
    } else {
        rwidth = int(size_t(thumbImg->getWidth()) * size_t(rheight) / size_t(thumbImg->getHeight()));
    }


    Imagefloat* baseImg = resizeTo<Imagefloat>(rwidth, rheight, interp, thumbImg);

    if (params.coarse.rotate) {
        baseImg->rotate (params.coarse.rotate);
        rwidth = baseImg->getWidth();
        rheight = baseImg->getHeight();
    }

    if (params.coarse.hflip) {
        baseImg->hflip ();
    }

    if (params.coarse.vflip) {
        baseImg->vflip ();
    }


    // apply white balance and raw white point (simulated)
    for (int i = 0; i < rheight; i++) {
#ifdef _OPENMP
        #pragma omp simd
#endif

        for (int j = 0; j < rwidth; j++) {
            float red = baseImg->r(i, j) * rmi;
            baseImg->r(i, j) = CLIP(red);
            float green = baseImg->g(i, j) * gmi;
            baseImg->g(i, j) = CLIP(green);
            float blue = baseImg->b(i, j) * bmi;
            baseImg->b(i, j) = CLIP(blue);

        }
    }

    // if luma denoise has to be done for thumbnails, it should be right here

    // perform color space transformation

    if (isRaw) {
        double pre_mul[3] = { redMultiplier, greenMultiplier, blueMultiplier };
        RawImageSource::colorSpaceConversion (baseImg, params.icm, currWB, pre_mul, embProfile, camProfile, cam2xyz, camName );
    } else {
        StdImageSource::colorSpaceConversion (baseImg, params.icm, embProfile, thumbImg->getSampleFormat());
    }

    int fw = baseImg->getWidth();
    int fh = baseImg->getHeight();
    //ColorTemp::CAT02 (baseImg, &params)   ;//perhaps not good!

    ImProcFunctions ipf (&params, false);
    ipf.setScale (sqrt(double(fw * fw + fh * fh)) / sqrt(double(thumbImg->getWidth() * thumbImg->getWidth() + thumbImg->getHeight() * thumbImg->getHeight()))*scale);
    ipf.updateColorProfiles (options.rtSettings.monitorProfile, options.rtSettings.monitorIntent, false, false);

    LUTu hist16 (65536);

    double gamma = isRaw ? Color::sRGBGamma : 0;  // usually in ImageSource, but we don't have that here
    ipf.firstAnalysis (baseImg, params, hist16);

    // perform transform
    if (ipf.needsTransform()) {
        Imagefloat* trImg = new Imagefloat (fw, fh);
        int origFW;
        int origFH;
        double tscale;
        getDimensions(origFW, origFH, tscale);
        ipf.transform (baseImg, trImg, 0, 0, 0, 0, fw, fh, origFW * tscale + 0.5, origFH * tscale + 0.5, focalLen, focalLen35mm, focusDist, 0, true); // Raw rotate degree not detectable here
        delete baseImg;
        baseImg = trImg;
    }

    // update blurmap
    SHMap* shmap = nullptr;

    if (params.sh.enabled) {
        shmap = new SHMap (fw, fh, false);
        double radius = sqrt (double(fw * fw + fh * fh)) / 2.0;
        double shradius = params.sh.radius;

        if (!params.sh.hq) {
            shradius *= radius / 1800.0;
        }

        shmap->update (baseImg, shradius, ipf.lumimul, params.sh.hq, 16);
    }

    // RGB processing
    double  expcomp = params.toneCurve.expcomp;
    int     bright = params.toneCurve.brightness;
    int     contr = params.toneCurve.contrast;
    int     black = params.toneCurve.black;
    int     hlcompr = params.toneCurve.hlcompr;
    int     hlcomprthresh = params.toneCurve.hlcomprthresh;

    if (params.toneCurve.autoexp && aeHistogram) {
        ipf.getAutoExp (aeHistogram, aeHistCompression, logDefGain, params.toneCurve.clip, expcomp, bright, contr, black, hlcompr, hlcomprthresh);
    }

    LUTf curve1 (65536);
    LUTf curve2 (65536);
    LUTf curve (65536);

    LUTf satcurve (65536);
    LUTf lhskcurve (65536);
    LUTf lumacurve(32770, 0); // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
    LUTf clcurve (65536);
    LUTf clToningcurve;
    LUTf cl2Toningcurve;
    LUTu dummy;

    ToneCurve customToneCurve1, customToneCurve2;
    ColorGradientCurve ctColorCurve;
    OpacityCurve ctOpacityCurve;

    ColorAppearance customColCurve1;
    ColorAppearance customColCurve2;
    ColorAppearance customColCurve3;
    ToneCurve customToneCurvebw1;
    ToneCurve customToneCurvebw2;
    CurveFactory::complexCurve (expcomp, black / 65535.0, hlcompr, hlcomprthresh,
                                params.toneCurve.shcompr, bright, contr,
                                params.toneCurve.curveMode, params.toneCurve.curve,
                                params.toneCurve.curveMode2, params.toneCurve.curve2,
                                hist16, curve1, curve2, curve, dummy, customToneCurve1, customToneCurve2, 16);

    LUTf rCurve;
    LUTf gCurve;
    LUTf bCurve;
    CurveFactory::RGBCurve (params.rgbCurves.rcurve, rCurve, 16);
    CurveFactory::RGBCurve (params.rgbCurves.gcurve, gCurve, 16);
    CurveFactory::RGBCurve (params.rgbCurves.bcurve, bCurve, 16);

    bool opautili = false;

    if(params.colorToning.enabled) {
        TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params.icm.working);
        double wp[3][3] = {
            {wprof[0][0], wprof[0][1], wprof[0][2]},
            {wprof[1][0], wprof[1][1], wprof[1][2]},
            {wprof[2][0], wprof[2][1], wprof[2][2]}
        };
        TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params.icm.working);
        double wip[3][3] = {
            {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
            {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
            {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
        };
        params.colorToning.getCurves(ctColorCurve, ctOpacityCurve, wp, wip, opautili);

        clToningcurve (65536);
        CurveFactory::curveToning(params.colorToning.clcurve, clToningcurve, scale == 1 ? 1 : 16);

        cl2Toningcurve (65536);
        CurveFactory::curveToning(params.colorToning.cl2curve, cl2Toningcurve, scale == 1 ? 1 : 16);
    }

    if(params.blackwhite.enabled) {
        CurveFactory::curveBW (params.blackwhite.beforeCurve, params.blackwhite.afterCurve, hist16, dummy, customToneCurvebw1, customToneCurvebw2, 16);
    }

    double rrm, ggm, bbm;
    float autor, autog, autob;
    float satLimit = float(params.colorToning.satProtectionThreshold) / 100.f * 0.7f + 0.3f;
    float satLimitOpacity = 1.f - (float(params.colorToning.saturatedOpacity) / 100.f);

    if(params.colorToning.enabled  && params.colorToning.autosat) { //for colortoning evaluation of saturation settings
        float moyS = 0.f;
        float eqty = 0.f;
        ipf.moyeqt (baseImg, moyS, eqty);//return image : mean saturation and standard dev of saturation
        //printf("moy=%f ET=%f\n", moyS,eqty);
        float satp = ((moyS + 1.5f * eqty) - 0.3f) / 0.7f; //1.5 sigma ==> 93% pixels with high saturation -0.3 / 0.7 convert to Hombre scale

        if(satp >= 0.92f) {
            satp = 0.92f;    //avoid values too high (out of gamut)
        }

        if(satp <= 0.15f) {
            satp = 0.15f;    //avoid too low values
        }

        satLimit = 100.f * satp;

        satLimitOpacity = 100.f * (moyS - 0.85f * eqty); //-0.85 sigma==>20% pixels with low saturation
    }

    autor = autog = autob = -9000.f; // This will ask to compute the "auto" values for the B&W tool

    LabImage* labView = new LabImage (fw, fh);
    DCPProfile *dcpProf = nullptr;
    DCPProfile::ApplyState as;

    if (isRaw) {
        cmsHPROFILE dummy;
        RawImageSource::findInputProfile(params.icm.input, nullptr, camName, &dcpProf, dummy);

        if (dcpProf) {
            dcpProf->setStep2ApplyState(params.icm.working, params.icm.toneCurve, params.icm.applyLookTable, params.icm.applyBaselineExposureOffset, as);
        }
    }

    LUTu histToneCurve;
    ipf.rgbProc (baseImg, labView, nullptr, curve1, curve2, curve, shmap, params.toneCurve.saturation, rCurve, gCurve, bCurve, satLimit , satLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2, customToneCurvebw1, customToneCurvebw2, rrm, ggm, bbm, autor, autog, autob, expcomp, hlcompr, hlcomprthresh, dcpProf, as, histToneCurve);

    // freeing up some memory
    customToneCurve1.Reset();
    customToneCurve2.Reset();
    ctColorCurve.Reset();
    ctOpacityCurve.Reset();
    customToneCurvebw1.Reset();
    customToneCurvebw2.Reset();

    if (shmap) {
        delete shmap;
    }

    // luminance histogram update
    if(params.labCurve.contrast != 0) {
        hist16.clear();

        for (int i = 0; i < fh; i++)
            for (int j = 0; j < fw; j++) {
                hist16[(int)((labView->L[i][j]))]++;
            }
    }

    // luminance processing
//  ipf.EPDToneMap(labView,0,6);

    bool utili;
    CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve,
                                 hist16, lumacurve, dummy, 16, utili);

    bool clcutili;
    CurveFactory::curveCL(clcutili, params.labCurve.clcurve, clcurve, 16);

    bool autili, butili, ccutili, cclutili;
    CurveFactory::complexsgnCurve (autili, butili, ccutili, cclutili, params.labCurve.acurve, params.labCurve.bcurve, params.labCurve.cccurve,
                                   params.labCurve.lccurve, curve1, curve2, satcurve, lhskcurve, 16);


	ipf.chromiLuminanceCurve (nullptr, 1, labView, labView, curve1, curve2, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, dummy, dummy);

    ipf.vibrance(labView);

    if((params.colorappearance.enabled && !params.colorappearance.tonecie) || !params.colorappearance.enabled) {
        ipf.EPDToneMap(labView, 5, 6);
    }

    if(params.colorappearance.enabled) {
        CurveFactory::curveLightBrightColor (
            params.colorappearance.curve,
            params.colorappearance.curve2,
            params.colorappearance.curve3,
            hist16, dummy,
            dummy, dummy,
            customColCurve1,
            customColCurve2,
            customColCurve3,
            16);
        int begh = 0, endh = labView->H;
        bool execsharp = false;
        float d;
        float fnum = fnumber;// F number
        float fiso = iso;// ISO
        float fspeed = shutter;//speed
        char * writ = new char[expcomp_.size() + 1];//convert expcomp_ to char
        std::copy(expcomp_.begin(), expcomp_.end(), writ);
        writ[expcomp_.size()] = '\0';
        float fcomp = atof(writ); //compensation + -
        delete[] writ;
        float adap;

        if(fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f)
            //if no exif data or wrong
        {
            adap = 2000.f;
        } else {
            float E_V = fcomp + log2 ((fnum * fnum) / fspeed / (fiso / 100.f));
            float expo2 = params.toneCurve.expcomp; // exposure compensation in tonecurve ==> direct EV
            E_V += expo2;
            float expo1;//exposure raw white point
            expo1 = log2(params.raw.expos); //log2 ==>linear to EV
            E_V += expo1;
            adap = powf(2.f, E_V - 3.f); //cd / m2
            //end calculation adaptation scene luminosity
        }

        LUTf CAMBrightCurveJ;
        LUTf CAMBrightCurveQ;
        float CAMMean;
        int sk;
        int scale;
        sk = 16;
        int rtt = 0;
        CieImage* cieView = new CieImage (fw, fh);
        CAMMean = NAN;
        CAMBrightCurveJ.dirty = true;
        CAMBrightCurveQ.dirty = true;
        ipf.ciecam_02float (cieView, adap, begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 6, execsharp, d, sk, rtt);
        delete cieView;
    }

    // color processing
    //ipf.colorCurve (labView, labView);

    // obtain final image
    Image8* readyImg = new Image8 (fw, fh);
    ipf.lab2monitorRgb (labView, readyImg);
    delete labView;
    delete baseImg;

    // calculate scale
    if (params.coarse.rotate == 90 || params.coarse.rotate == 270) {
        myscale = scale * thumbImg->getWidth() / fh;
    } else {
        myscale = scale * thumbImg->getHeight() / fh;
    }

    myscale = 1.0 / myscale;
    /*    // apply crop
        if (params.crop.enabled) {
            int ix = 0;
            for (int i=0; i<fh; i++)
                for (int j=0; j<fw; j++)
                    if (i<params.crop.y/myscale || i>(params.crop.y+params.crop.h)/myscale || j<params.crop.x/myscale || j>(params.crop.x+params.crop.w)/myscale) {
                        readyImg->data[ix++] /= 3;
                        readyImg->data[ix++] /= 3;
                        readyImg->data[ix++] /= 3;
                    }
                    else
                        ix += 3;
        }*/

    return readyImg;
}

int Thumbnail::getImageWidth (const procparams::ProcParams& params, int rheight, float &ratio)
{
    if (!thumbImg) {
        return 0;    // Can happen if thumb is just building and GUI comes in with resize wishes
    }

    int rwidth;

    if (params.coarse.rotate == 90 || params.coarse.rotate == 270) {
        ratio = (float)(thumbImg->getHeight()) / (float)(thumbImg->getWidth());
    } else {
        ratio = (float)(thumbImg->getWidth()) / (float)(thumbImg->getHeight());
    }

    rwidth = (int)(ratio * (float)rheight);

    return rwidth;
}

void Thumbnail::getDimensions (int& w, int& h, double& scaleFac)
{
    if (thumbImg) {
        w = thumbImg->getWidth();
        h = thumbImg->getHeight();
        scaleFac = scale;
    } else {
        w = 0;
        h = 0;
        scale = 1;
    }
}

void Thumbnail::getCamWB (double& temp, double& green)
{

    double cam_r = colorMatrix[0][0] * camwbRed + colorMatrix[0][1] * camwbGreen + colorMatrix[0][2] * camwbBlue;
    double cam_g = colorMatrix[1][0] * camwbRed + colorMatrix[1][1] * camwbGreen + colorMatrix[1][2] * camwbBlue;
    double cam_b = colorMatrix[2][0] * camwbRed + colorMatrix[2][1] * camwbGreen + colorMatrix[2][2] * camwbBlue;
    ColorTemp currWB = ColorTemp (cam_r, cam_g, cam_b, 1.0);  // we do not take the equalizer into account here, because we want camera's WB
    temp = currWB.getTemp ();
    green = currWB.getGreen ();
}

void Thumbnail::getAutoWB (double& temp, double& green, double equal, double tempBias)
{

    if (equal != wbEqual || tempBias != wbTempBias) {
        // compute the values depending on equal
        ColorTemp cTemp;
        wbEqual = equal;
        wbTempBias = tempBias;
        // compute autoWBTemp and autoWBGreen
        cTemp.mul2temp(redAWBMul, greenAWBMul, blueAWBMul, wbEqual, autoWBTemp, autoWBGreen);
        autoWBTemp += autoWBTemp * tempBias;
    }

    temp = autoWBTemp;
    green = autoWBGreen;
}

void Thumbnail::getAutoWBMultipliers (double& rm, double& gm, double& bm)
{
    rm = redAWBMul;
    gm = greenAWBMul;
    bm = blueAWBMul;
}

void Thumbnail::applyAutoExp (procparams::ProcParams& params)
{

    if (params.toneCurve.autoexp && aeHistogram) {
        ImProcFunctions ipf (&params, false);
        ipf.getAutoExp (aeHistogram, aeHistCompression, log(defGain) / log(2.0), params.toneCurve.clip, params.toneCurve.expcomp,
                        params.toneCurve.brightness, params.toneCurve.contrast, params.toneCurve.black, params.toneCurve.hlcompr, params.toneCurve.hlcomprthresh);
    }
}

void Thumbnail::getSpotWB (const procparams::ProcParams& params, int xp, int yp, int rect, double& rtemp, double& rgreen)
{

    std::vector<Coord2D> points, red, green, blue;

    for (int i = yp - rect; i <= yp + rect; i++)
        for (int j = xp - rect; j <= xp + rect; j++) {
            points.push_back (Coord2D (j, i));
        }

    int fw = thumbImg->getWidth(), fh = thumbImg->getHeight();

    if (params.coarse.rotate == 90 || params.coarse.rotate == 270) {
        fw = thumbImg->getHeight();
        fh = thumbImg->getWidth();
    }

    ImProcFunctions ipf (&params, false);
    ipf.transCoord (fw, fh, points, red, green, blue);
    int tr = getCoarseBitMask(params.coarse);
    // calculate spot wb (copy & pasted from stdimagesource)
    double reds = 0, greens = 0, blues = 0;
    int rn = 0, gn = 0, bn = 0;
    thumbImg->getSpotWBData(reds, greens, blues, rn, gn, bn, red, green, blue, tr);
    reds = reds / rn * camwbRed;
    greens = greens / gn * camwbGreen;
    blues = blues / bn * camwbBlue;

    double rm = colorMatrix[0][0] * reds + colorMatrix[0][1] * greens + colorMatrix[0][2] * blues;
    double gm = colorMatrix[1][0] * reds + colorMatrix[1][1] * greens + colorMatrix[1][2] * blues;
    double bm = colorMatrix[2][0] * reds + colorMatrix[2][1] * greens + colorMatrix[2][2] * blues;

    ColorTemp ct (rm, gm, bm, params.wb.equal);
    rtemp = ct.getTemp ();
    rgreen = ct.getGreen ();
}
void Thumbnail::transformPixel (int x, int y, int tran, int& tx, int& ty)
{

    int W = thumbImg->getWidth();
    int H = thumbImg->getHeight();
    int sw = W, sh = H;

    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = H;
        sh = W;
    }

    int ppx = x, ppy = y;

    if (tran & TR_HFLIP) {
        ppx = sw - 1 - x ;
    }

    if (tran & TR_VFLIP) {
        ppy = sh - 1 - y;
    }

    tx = ppx;
    ty = ppy;

    if ((tran & TR_ROT) == TR_R180) {
        tx = W - 1 - ppx;
        ty = H - 1 - ppy;
    } else if ((tran & TR_ROT) == TR_R90) {
        tx = ppy;
        ty = H - 1 - ppx;
    } else if ((tran & TR_ROT) == TR_R270) {
        tx = W - 1 - ppy;
        ty = ppx;
    }

    tx /= scale;
    ty /= scale;
}

unsigned char* Thumbnail::getGrayscaleHistEQ (int trim_width)
{
    if (!thumbImg) {
        return nullptr;
    }

    if (thumbImg->getWidth() < trim_width) {
        return nullptr;
    }

    // to utilize the 8 bit color range of the thumbnail we brighten it and apply gamma correction
    unsigned char* tmpdata = new unsigned char[thumbImg->getHeight() * trim_width];
    int ix = 0, max;

    if (gammaCorrected) {
        // if it's gamma correct (usually a RAW), we have the problem that there is a lot noise etc. that makes the maximum way too high.
        // Strategy is limit a certain percent of pixels so the overall picture quality when scaling to 8 bit is way better
        const double BurnOffPct = 0.03; // *100 = percent pixels that may be clipped

        // Calc the histogram
        unsigned int* hist16 = new unsigned int [65536];
        memset(hist16, 0, sizeof(int) * 65536);

        if (thumbImg->getType() == sImage8) {
            Image8 *image = static_cast<Image8*>(thumbImg);
            image->calcGrayscaleHist(hist16);
        } else if (thumbImg->getType() == sImage16) {
            Image16 *image = static_cast<Image16*>(thumbImg);
            image->calcGrayscaleHist(hist16);
        } else if (thumbImg->getType() == sImagefloat) {
            Imagefloat *image = static_cast<Imagefloat*>(thumbImg);
            image->calcGrayscaleHist(hist16);
        } else {
            printf("getGrayscaleHistEQ #1: Unsupported image type \"%s\"!\n", thumbImg->getType());
        }

        // Go down till we cut off that many pixels
        unsigned long cutoff = thumbImg->getHeight() * thumbImg->getHeight() * 4 * BurnOffPct;

        int max_;
        unsigned long sum = 0;

        for (max_ = 65535; max_ > 16384 && sum < cutoff; max_--) {
            sum += hist16[max_];
        }

        delete[] hist16;

        scaleForSave = 65535 * 8192 / max_;

        // Correction and gamma to 8 Bit
        if (thumbImg->getType() == sImage8) {
            Image8 *image = static_cast<Image8*>(thumbImg);

            for (int i = 0; i < thumbImg->getHeight(); i++)
                for (int j = (thumbImg->getWidth() - trim_width) / 2; j < trim_width + (thumbImg->getWidth() - trim_width) / 2; j++) {
                    unsigned short r_, g_, b_;
                    image->convertTo(image->r(i, j), r_);
                    image->convertTo(image->g(i, j), g_);
                    image->convertTo(image->b(i, j), b_);
                    int r = Color::gammatabThumb[min(r_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    int g = Color::gammatabThumb[min(g_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    int b = Color::gammatabThumb[min(b_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    tmpdata[ix++] = (r * 19595 + g * 38469 + b * 7472) >> 16;
                }
        } else if (thumbImg->getType() == sImage16) {
            Image16 *image = static_cast<Image16*>(thumbImg);

            for (int i = 0; i < thumbImg->getHeight(); i++)
                for (int j = (thumbImg->getWidth() - trim_width) / 2; j < trim_width + (thumbImg->getWidth() - trim_width) / 2; j++) {
                    unsigned short r_, g_, b_;
                    image->convertTo(image->r(i, j), r_);
                    image->convertTo(image->g(i, j), g_);
                    image->convertTo(image->b(i, j), b_);
                    int r = Color::gammatabThumb[min(r_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    int g = Color::gammatabThumb[min(g_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    int b = Color::gammatabThumb[min(b_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    tmpdata[ix++] = (r * 19595 + g * 38469 + b * 7472) >> 16;
                }
        } else if (thumbImg->getType() == sImagefloat) {
            Imagefloat *image = static_cast<Imagefloat*>(thumbImg);

            for (int i = 0; i < thumbImg->getHeight(); i++)
                for (int j = (thumbImg->getWidth() - trim_width) / 2; j < trim_width + (thumbImg->getWidth() - trim_width) / 2; j++) {
                    unsigned short r_, g_, b_;
                    image->convertTo(image->r(i, j), r_);
                    image->convertTo(image->g(i, j), g_);
                    image->convertTo(image->b(i, j), b_);
                    int r = Color::gammatabThumb[min(r_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    int g = Color::gammatabThumb[min(g_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    int b = Color::gammatabThumb[min(b_, static_cast<unsigned short>(max_)) * scaleForSave >> 13];
                    tmpdata[ix++] = (r * 19595 + g * 38469 + b * 7472) >> 16;
                }
        }
    } else {
        // If it's not gamma corrected (usually a JPG) we take the normal maximum
        max = 0;

        if (thumbImg->getType() == sImage8) {
            Image8 *image = static_cast<Image8*>(thumbImg);
            unsigned char max_ = 0;

            for (int row = 0; row < image->getHeight(); row++)
                for (int col = 0; col < image->getWidth(); col++) {
                    if (image->r(row, col) > max_) {
                        max_ = image->r(row, col);
                    }

                    if (image->g(row, col) > max_) {
                        max_ = image->g(row, col);
                    }

                    if (image->b(row, col) > max_) {
                        max_ = image->b(row, col);
                    }
                }

            image->convertTo(max_, max);

            if (max < 16384) {
                max = 16384;
            }

            scaleForSave = 65535 * 8192 / max;

            // Correction and gamma to 8 Bit
            for (int i = 0; i < image->getHeight(); i++)
                for (int j = (image->getWidth() - trim_width) / 2; j < trim_width + (image->getWidth() - trim_width) / 2; j++) {
                    unsigned short rtmp, gtmp, btmp;
                    image->convertTo(image->r(i, j), rtmp);
                    image->convertTo(image->g(i, j), gtmp);
                    image->convertTo(image->b(i, j), btmp);
                    int r = rtmp * scaleForSave >> 21;
                    int g = gtmp * scaleForSave >> 21;
                    int b = btmp * scaleForSave >> 21;
                    tmpdata[ix++] = (r * 19595 + g * 38469 + b * 7472) >> 16;
                }
        } else if (thumbImg->getType() == sImage16) {
            Image16 *image = static_cast<Image16*>(thumbImg);
            unsigned short max_ = 0;

            for (int row = 0; row < image->getHeight(); row++)
                for (int col = 0; col < image->getWidth(); col++) {
                    if (image->r(row, col) > max_) {
                        max_ = image->r(row, col);
                    }

                    if (image->g(row, col) > max_) {
                        max_ = image->g(row, col);
                    }

                    if (image->b(row, col) > max_) {
                        max_ = image->b(row, col);
                    }
                }

            image->convertTo(max_, max);

            if (max < 16384) {
                max = 16384;
            }

            scaleForSave = 65535 * 8192 / max;

            // Correction and gamma to 8 Bit
            for (int i = 0; i < image->getHeight(); i++)
                for (int j = (image->getWidth() - trim_width) / 2; j < trim_width + (image->getWidth() - trim_width) / 2; j++) {
                    unsigned short rtmp, gtmp, btmp;
                    image->convertTo(image->r(i, j), rtmp);
                    image->convertTo(image->g(i, j), gtmp);
                    image->convertTo(image->b(i, j), btmp);
                    int r = rtmp * scaleForSave >> 21;
                    int g = gtmp * scaleForSave >> 21;
                    int b = btmp * scaleForSave >> 21;
                    tmpdata[ix++] = (r * 19595 + g * 38469 + b * 7472) >> 16;
                }
        } else if (thumbImg->getType() == sImagefloat) {
            Imagefloat *image = static_cast<Imagefloat*>(thumbImg);
            float max_ = 0.f;

            for (int row = 0; row < image->getHeight(); row++)
                for (int col = 0; col < image->getWidth(); col++) {
                    if (image->r(row, col) > max_) {
                        max_ = image->r(row, col);
                    }

                    if (image->g(row, col) > max_) {
                        max_ = image->g(row, col);
                    }

                    if (image->b(row, col) > max_) {
                        max_ = image->b(row, col);
                    }
                }

            image->convertTo(max_, max);

            if (max < 16384) {
                max = 16384;
            }

            scaleForSave = 65535 * 8192 / max;

            // Correction and gamma to 8 Bit
            for (int i = 0; i < image->getHeight(); i++)
                for (int j = (image->getWidth() - trim_width) / 2; j < trim_width + (image->getWidth() - trim_width) / 2; j++) {
                    unsigned short rtmp, gtmp, btmp;
                    image->convertTo(image->r(i, j), rtmp);
                    image->convertTo(image->g(i, j), gtmp);
                    image->convertTo(image->b(i, j), btmp);
                    int r = rtmp * scaleForSave >> 21;
                    int g = gtmp * scaleForSave >> 21;
                    int b = btmp * scaleForSave >> 21;
                    tmpdata[ix++] = (r * 19595 + g * 38469 + b * 7472) >> 16;
                }
        } else {
            printf("getGrayscaleHistEQ #2: Unsupported image type \"%s\"!\n", thumbImg->getType());
        }
    }

    // histogram equalization
    unsigned int hist[256] = {0};

    for (int i = 0; i < ix; i++) {
        hist[tmpdata[i]]++;
    }

    int cdf = 0, cdf_min = -1;

    for (int i = 0; i < 256; i++) {
        cdf += hist[i];

        if (cdf > 0 && cdf_min == -1) {
            cdf_min = cdf;
        }

        if (cdf_min != -1) {
            hist[i] = (cdf - cdf_min) * 255 / ((thumbImg->getHeight() * trim_width) - cdf_min);
        }
    }

    for (int i = 0; i < ix; i++) {
        tmpdata[i] = hist[tmpdata[i]];
    }

    return tmpdata;
}

bool Thumbnail::writeImage (const Glib::ustring& fname, int format)
{

    if (!thumbImg) {
        return false;
    }

    Glib::ustring fullFName = fname + ".rtti";

    FILE* f = g_fopen (fullFName.c_str (), "wb");

    if (!f) {
        return false;
    }

    fwrite (thumbImg->getType(), sizeof (char), strlen(thumbImg->getType()), f);
    fputc ('\n', f);
    guint32 w = guint32(thumbImg->getWidth());
    guint32 h = guint32(thumbImg->getHeight());
    fwrite (&w, sizeof (guint32), 1, f);
    fwrite (&h, sizeof (guint32), 1, f);

    if (thumbImg->getType() == sImage8) {
        Image8 *image = static_cast<Image8*>(thumbImg);
        image->writeData(f);
    } else if (thumbImg->getType() == sImage16) {
        Image16 *image = static_cast<Image16*>(thumbImg);
        image->writeData(f);
    } else if (thumbImg->getType() == sImagefloat) {
        Imagefloat *image = static_cast<Imagefloat*>(thumbImg);
        image->writeData(f);
    }

    //thumbImg->writeData(f);
    fclose (f);
    return true;
}

bool Thumbnail::readImage (const Glib::ustring& fname)
{

    if (thumbImg) {
        delete thumbImg;
        thumbImg = nullptr;
    }

    Glib::ustring fullFName = fname + ".rtti";

    if (!Glib::file_test (fullFName, Glib::FILE_TEST_EXISTS)) {
        return false;
    }

    FILE* f = g_fopen (fullFName.c_str (), "rb");

    if (!f) {
        return false;
    }

    char imgType[31];  // 30 -> arbitrary size, but should be enough for all image type's name
    fgets(imgType, 30, f);
    imgType[strlen(imgType) - 1] = '\0'; // imgType has a \n trailing character, so we overwrite it by the \0 char

    guint32 width, height;
    fread (&width, 1, sizeof (guint32), f);
    fread (&height, 1, sizeof (guint32), f);

    bool success = false;

    if (!strcmp(imgType, sImage8)) {
        Image8 *image = new Image8(width, height);
        image->readData(f);
        thumbImg = image;
        success = true;
    } else if (!strcmp(imgType, sImage16)) {
        Image16 *image = new Image16(width, height);
        image->readData(f);
        thumbImg = image;
        success = true;
    } else if (!strcmp(imgType, sImagefloat)) {
        Imagefloat *image = new Imagefloat(width, height);
        image->readData(f);
        thumbImg = image;
        success = true;
    } else {
        printf("readImage: Unsupported image type \"%s\"!\n", imgType);
    }

    fclose(f);
    return success;
}

bool Thumbnail::readData  (const Glib::ustring& fname)
{
    setlocale(LC_NUMERIC, "C"); // to set decimal point to "."
    Glib::KeyFile keyFile;

    try {
        MyMutex::MyLock thmbLock(thumbMutex);

        try {
            keyFile.load_from_file (fname);
        } catch (Glib::Error&) {
            return false;
        }

        if (keyFile.has_group ("LiveThumbData")) {
            if (keyFile.has_key ("LiveThumbData", "CamWBRed")) {
                camwbRed            = keyFile.get_double ("LiveThumbData", "CamWBRed");
            }

            if (keyFile.has_key ("LiveThumbData", "CamWBGreen")) {
                camwbGreen          = keyFile.get_double ("LiveThumbData", "CamWBGreen");
            }

            if (keyFile.has_key ("LiveThumbData", "CamWBBlue")) {
                camwbBlue           = keyFile.get_double ("LiveThumbData", "CamWBBlue");
            }

            if (keyFile.has_key ("LiveThumbData", "RedAWBMul")) {
                redAWBMul           = keyFile.get_double ("LiveThumbData", "RedAWBMul");
            }

            if (keyFile.has_key ("LiveThumbData", "GreenAWBMul")) {
                greenAWBMul         = keyFile.get_double ("LiveThumbData", "GreenAWBMul");
            }

            if (keyFile.has_key ("LiveThumbData", "BlueAWBMul")) {
                blueAWBMul          = keyFile.get_double ("LiveThumbData", "BlueAWBMul");
            }

            if (keyFile.has_key ("LiveThumbData", "AEHistCompression")) {
                aeHistCompression   = keyFile.get_integer ("LiveThumbData", "AEHistCompression");
            }

            if (keyFile.has_key ("LiveThumbData", "RedMultiplier")) {
                redMultiplier       = keyFile.get_double ("LiveThumbData", "RedMultiplier");
            }

            if (keyFile.has_key ("LiveThumbData", "GreenMultiplier")) {
                greenMultiplier     = keyFile.get_double ("LiveThumbData", "GreenMultiplier");
            }

            if (keyFile.has_key ("LiveThumbData", "BlueMultiplier")) {
                blueMultiplier      = keyFile.get_double ("LiveThumbData", "BlueMultiplier");
            }

            if (keyFile.has_key ("LiveThumbData", "Scale")) {
                scale               = keyFile.get_double ("LiveThumbData", "Scale");
            }

            if (keyFile.has_key ("LiveThumbData", "DefaultGain")) {
                defGain             = keyFile.get_double ("LiveThumbData", "DefaultGain");
            }

            if (keyFile.has_key ("LiveThumbData", "ScaleForSave")) {
                scaleForSave        = keyFile.get_integer ("LiveThumbData", "ScaleForSave");
            }

            if (keyFile.has_key ("LiveThumbData", "GammaCorrected")) {
                gammaCorrected      = keyFile.get_boolean ("LiveThumbData", "GammaCorrected");
            }

            if (keyFile.has_key ("LiveThumbData", "ColorMatrix")) {
                std::vector<double> cm = keyFile.get_double_list ("LiveThumbData", "ColorMatrix");
                int ix = 0;

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
                        colorMatrix[i][j] = cm[ix++];
                    }
            }
        }

        return true;
    } catch (Glib::Error &err) {
        if (options.rtSettings.verbose) {
            printf("Thumbnail::readData / Error code %d while reading values from \"%s\":\n%s\n", err.code(), fname.c_str(), err.what().c_str());
        }
    } catch (...) {
        if (options.rtSettings.verbose) {
            printf("Thumbnail::readData / Unknown exception while trying to load \"%s\"!\n", fname.c_str());
        }
    }

    return false;
}

bool Thumbnail::writeData  (const Glib::ustring& fname)
{
    MyMutex::MyLock thmbLock(thumbMutex);

    Glib::ustring keyData;

    try {

        Glib::KeyFile keyFile;

        try {
            keyFile.load_from_file (fname);
        } catch (Glib::Error&) {}

        keyFile.set_double  ("LiveThumbData", "CamWBRed", camwbRed);
        keyFile.set_double  ("LiveThumbData", "CamWBGreen", camwbGreen);
        keyFile.set_double  ("LiveThumbData", "CamWBBlue", camwbBlue);
        keyFile.set_double  ("LiveThumbData", "RedAWBMul", redAWBMul);
        keyFile.set_double  ("LiveThumbData", "GreenAWBMul", greenAWBMul);
        keyFile.set_double  ("LiveThumbData", "BlueAWBMul", blueAWBMul);
        keyFile.set_integer ("LiveThumbData", "AEHistCompression", aeHistCompression);
        keyFile.set_double  ("LiveThumbData", "RedMultiplier", redMultiplier);
        keyFile.set_double  ("LiveThumbData", "GreenMultiplier", greenMultiplier);
        keyFile.set_double  ("LiveThumbData", "BlueMultiplier", blueMultiplier);
        keyFile.set_double  ("LiveThumbData", "Scale", scale);
        keyFile.set_double  ("LiveThumbData", "DefaultGain", defGain);
        keyFile.set_integer ("LiveThumbData", "ScaleForSave", scaleForSave);
        keyFile.set_boolean ("LiveThumbData", "GammaCorrected", gammaCorrected);
        Glib::ArrayHandle<double> cm ((double*)colorMatrix, 9, Glib::OWNERSHIP_NONE);
        keyFile.set_double_list ("LiveThumbData", "ColorMatrix", cm);

        keyData = keyFile.to_data ();

    } catch (Glib::Error& err) {
        if (options.rtSettings.verbose) {
            printf("Thumbnail::writeData / Error code %d while reading values from \"%s\":\n%s\n", err.code(), fname.c_str(), err.what().c_str());
        }
    } catch (...) {
        if (options.rtSettings.verbose) {
            printf("Thumbnail::writeData / Unknown exception while trying to save \"%s\"!\n", fname.c_str());
        }
    }

    if (keyData.empty ()) {
        return false;
    }

    FILE *f = g_fopen (fname.c_str (), "wt");

    if (!f) {
        if (options.rtSettings.verbose) {
            printf("Thumbnail::writeData / Error: unable to open file \"%s\" with write access!\n", fname.c_str());
        }

        return false;
    } else {
        fprintf (f, "%s", keyData.c_str ());
        fclose (f);
    }

    return true;
}

bool Thumbnail::readEmbProfile  (const Glib::ustring& fname)
{

    FILE* f = g_fopen (fname.c_str (), "rb");

    if (!f) {
        embProfileData = nullptr;
        embProfile = nullptr;
        embProfileLength = 0;
    } else {
        fseek (f, 0, SEEK_END);
        embProfileLength = ftell (f);
        fseek (f, 0, SEEK_SET);
        embProfileData = new unsigned char[embProfileLength];
        fread (embProfileData, 1, embProfileLength, f);
        fclose (f);
        embProfile = cmsOpenProfileFromMem (embProfileData, embProfileLength);
        return true;
    }

    return false;
}

bool Thumbnail::writeEmbProfile (const Glib::ustring& fname)
{

    if (embProfileData) {
        FILE* f = g_fopen(fname.c_str (), "wb");

        if (f) {
            fwrite (embProfileData, 1, embProfileLength, f);
            fclose (f);
            return true;
        }
    }

    return false;
}

bool Thumbnail::readAEHistogram  (const Glib::ustring& fname)
{

    FILE* f = g_fopen (fname.c_str (), "rb");

    if (!f) {
        aeHistogram(0);
    } else {
        aeHistogram(65536 >> aeHistCompression);
        fread (&aeHistogram[0], 1, (65536 >> aeHistCompression)*sizeof(aeHistogram[0]), f);
        fclose (f);
        return true;
    }

    return false;
}

bool Thumbnail::writeAEHistogram (const Glib::ustring& fname)
{

    if (aeHistogram) {
        FILE* f = g_fopen (fname.c_str (), "wb");

        if (f) {
            fwrite (&aeHistogram[0], 1, (65536 >> aeHistCompression)*sizeof(aeHistogram[0]), f);
            fclose (f);
            return true;
        }
    }

    return false;
}

unsigned char* Thumbnail::getImage8Data()
{
    if (thumbImg && thumbImg->getType() == rtengine::sImage8) {
        Image8* img8 = static_cast<Image8*>(thumbImg);
        return img8->data;
    }

    return nullptr;
}



}
