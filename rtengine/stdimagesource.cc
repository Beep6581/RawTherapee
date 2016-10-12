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
#include "stdimagesource.h"
#include "mytime.h"
#include "iccstore.h"
#include "imageio.h"
#include "curves.h"
#include "color.h"

#undef THREAD_PRIORITY_NORMAL

namespace rtengine
{

extern const Settings* settings;

template<class T> void freeArray (T** a, int H)
{
    for (int i = 0; i < H; i++) {
        delete [] a[i];
    }

    delete [] a;
}
template<class T> T** allocArray (int W, int H)
{

    T** t = new T*[H];

    for (int i = 0; i < H; i++) {
        t[i] = new T[W];
    }

    return t;
}

#define HR_SCALE 2
StdImageSource::StdImageSource () : ImageSource(), img(NULL), plistener(NULL), full(false), max{}, rgbSourceModified(false)
{

    hrmap[0] = NULL;
    hrmap[1] = NULL;
    hrmap[2] = NULL;
    needhr = NULL;
    embProfile = NULL;
    idata = NULL;
}

StdImageSource::~StdImageSource ()
{

    delete idata;

    if (hrmap[0] != NULL) {
        int dh = img->getH() / HR_SCALE;
        freeArray<float>(hrmap[0], dh);
        freeArray<float>(hrmap[1], dh);
        freeArray<float>(hrmap[2], dh);
    }

    if (needhr) {
        freeArray<char>(needhr, img->getH());
    }

    if (img) {
        delete img;
    }
}

void StdImageSource::getSampleFormat (const Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement)
{

    sFormat = IIOSF_UNKNOWN;
    sArrangement = IIOSA_UNKNOWN;

    size_t lastdot = fname.find_last_of ('.');

    if( Glib::ustring::npos == lastdot ) {
        return;
    }

    if (!fname.casefold().compare (lastdot, 4, ".jpg") ||
            !fname.casefold().compare (lastdot, 5, ".jpeg")) {
        // For now, png and jpeg files are converted to unsigned short by the loader itself,
        // but there should be functions that read the sample format first, like the TIFF case below
        sFormat = IIOSF_UNSIGNED_CHAR;
        sArrangement = IIOSA_CHUNKY;
        return;
    } else if (!fname.casefold().compare (lastdot, 4, ".png")) {
        int result = ImageIO::getPNGSampleFormat (fname, sFormat, sArrangement);

        if (result == IMIO_SUCCESS) {
            return;
        }
    } else if (!fname.casefold().compare (lastdot, 4, ".tif") ||
               !fname.casefold().compare (lastdot, 5, ".tiff")) {
        int result = ImageIO::getTIFFSampleFormat (fname, sFormat, sArrangement);

        if (result == IMIO_SUCCESS) {
            return;
        }
    }

    return;
}

/*
 * This method make define the correspondence between the input image type
 * and RT's image data type (Image8, Image16 and Imagefloat), then it will
 * load the image into it
 */
int StdImageSource::load (const Glib::ustring &fname, bool batch)
{

    fileName = fname;

    // First let's find out the input image's type

    IIOSampleFormat sFormat;
    IIOSampleArrangement sArrangement;
    getSampleFormat(fname, sFormat, sArrangement);

    // Then create the appropriate object

    switch (sFormat) {
    case (IIOSF_UNSIGNED_CHAR): {
        Image8 *img_8 = new Image8 ();
        img = img_8;
        break;
    }

    case (IIOSF_UNSIGNED_SHORT): {
        Image16 *img_16 = new Image16 ();
        img = img_16;
        break;
    }

    case (IIOSF_LOGLUV24):
    case (IIOSF_LOGLUV32):
    case (IIOSF_FLOAT): {
        Imagefloat *img_float = new Imagefloat ();
        img = img_float;
        break;
    }

    default:
        return IMIO_FILETYPENOTSUPPORTED;
    }

    img->setSampleFormat(sFormat);
    img->setSampleArrangement(sArrangement);

    if (plistener) {
        plistener->setProgressStr ("PROGRESSBAR_LOADING");
        plistener->setProgress (0.0);
        img->setProgressListener (plistener);
    }

    // And load the image!

    int error = img->load (fname);

    if (error) {
        delete img;
        img = NULL;
        return error;
    }

    embProfile = img->getEmbeddedProfile ();

    idata = new ImageData (fname);

    if (idata->hasExif()) {
        int deg = 0;

        if (idata->getOrientation() == "Rotate 90 CW") {
            deg = 90;
        } else if (idata->getOrientation() == "Rotate 180") {
            deg = 180;
        } else if (idata->getOrientation() == "Rotate 270 CW") {
            deg = 270;
        }

        if (deg) {
            img->rotate(deg);
        }
    }

    if (plistener) {
        plistener->setProgressStr ("PROGRESSBAR_READY");
        plistener->setProgress (1.0);
    }

    wb = ColorTemp (1.0, 1.0, 1.0, 1.0);
    //this is probably a mistake if embedded profile is not D65

    return 0;
}

void StdImageSource::getImage (const ColorTemp &ctemp, int tran, Imagefloat* image, const PreviewProps &pp, const ToneCurveParams &hrp, const ColorManagementParams &cmp, const RAWParams &raw)
{

    // the code will use OpenMP as of now.

    img->getStdImage(ctemp, tran, image, pp, true, hrp);

    // Hombre: we could have rotated the image here too, with just few line of code, but:
    // 1. it would require other modifications in the engine, so "do not touch that little plonker!"
    // 2. it's more optimized like this

    // Flip if needed
    if (tran & TR_HFLIP) {
        image->hflip();
    }

    if (tran & TR_VFLIP) {
        image->vflip();
    }
}

void StdImageSource::convertColorSpace(Imagefloat* image, const ColorManagementParams &cmp, const ColorTemp &wb)
{
    colorSpaceConversion (image, cmp, embProfile, img->getSampleFormat());
}

void StdImageSource::colorSpaceConversion (Imagefloat* im, const ColorManagementParams &cmp, cmsHPROFILE embedded, IIOSampleFormat sampleFormat)
{

    bool skipTransform = false;
    cmsHPROFILE in = NULL;
    cmsHPROFILE out = iccStore->workingSpace (cmp.working);

    if (cmp.input == "(embedded)" || cmp.input == "" || cmp.input == "(camera)" || cmp.input == "(cameraICC)") {
        if (embedded) {
            in = embedded;
        } else {
            if (sampleFormat & (IIOSF_LOGLUV24 | IIOSF_LOGLUV32 | IIOSF_FLOAT)) {
                skipTransform = true;
            } else {
                in = iccStore->getsRGBProfile ();
            }
        }
    } else {
        if (cmp.input != "(none)") {
            in = iccStore->getProfile (cmp.input);

            if (in == NULL && embedded) {
                in = embedded;
            } else if (in == NULL) {
                if (sampleFormat & (IIOSF_LOGLUV24 | IIOSF_LOGLUV32 | IIOSF_FLOAT)) {
                    skipTransform = true;
                } else {
                    in = iccStore->getsRGBProfile ();
                }
            }
        }
    }

    if (!skipTransform && in) {
        if(in == embedded && cmsGetColorSpace(in) != cmsSigRgbData) { // if embedded profile is not an RGB profile, use sRGB
            printf("embedded profile is not an RGB profile, using sRGB as input profile\n");
            in = iccStore->getsRGBProfile ();
        }

        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (in, TYPE_RGB_FLT, out, TYPE_RGB_FLT, INTENT_RELATIVE_COLORIMETRIC,
                                   cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
        lcmsMutex->unlock ();

        if(hTransform) {
            // Convert to the [0.0 ; 1.0] range
            im->normalizeFloatTo1();

            im->ExecCMSTransform(hTransform);

            // Converting back to the [0.0 ; 65535.0] range
            im->normalizeFloatTo65535();

            cmsDeleteTransform(hTransform);
        } else {
            printf("Could not convert from %s to %s\n", in == embedded ? "embedded profile" : cmp.input.data(), cmp.working.data());
        }
    }
}

void StdImageSource::getFullSize (int& w, int& h, int tr)
{

    w = img->width;
    h = img->height;

    if ((tr & TR_ROT) == TR_R90 || (tr & TR_ROT) == TR_R270) {
        w = img->height;
        h = img->width;
    }
}

void StdImageSource::getSize (int tran, PreviewProps pp, int& w, int& h)
{

    w = pp.w / pp.skip + (pp.w % pp.skip > 0);
    h = pp.h / pp.skip + (pp.h % pp.skip > 0);
}

void StdImageSource::getAutoExpHistogram (LUTu & histogram, int& histcompr)
{
    if (img->getType() == sImage8) {
        Image8 *img_ = static_cast<Image8*>(img);
        img_->computeAutoHistogram(histogram, histcompr);
    } else if (img->getType() == sImage16) {
        Image16 *img_ = static_cast<Image16*>(img);
        img_->computeAutoHistogram(histogram, histcompr);
    } else if (img->getType() == sImagefloat) {
        Imagefloat *img_ = static_cast<Imagefloat*>(img);
        img_->computeAutoHistogram(histogram, histcompr);
    }
}

void StdImageSource::getAutoWBMultipliers (double &rm, double &gm, double &bm)
{
    if (redAWBMul != -1.) {
        rm = redAWBMul;
        gm = greenAWBMul;
        bm = blueAWBMul;
        return;
    }

    img->getAutoWBMultipliers(rm, gm, bm);

    redAWBMul   = rm;
    greenAWBMul = gm;
    blueAWBMul  = bm;
}

ColorTemp StdImageSource::getSpotWB (std::vector<Coord2D> &red, std::vector<Coord2D> &green, std::vector<Coord2D>& blue, int tran, double equal)
{
    int rn, gn, bn;
    double reds, greens, blues;
    img->getSpotWBData(reds, greens, blues, rn, gn, bn, red, green, blue, tran);
    double img_r, img_g, img_b;
    wb.getMultipliers (img_r, img_g, img_b);

    if( settings->verbose ) {
        printf ("AVG: %g %g %g\n", reds / rn, greens / gn, blues / bn);
    }

    return ColorTemp (reds / rn * img_r, greens / gn * img_g, blues / bn * img_b, equal);
}

}

