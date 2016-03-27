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

#include <tiff.h>
#include <tiffio.h>

#include "image8.h"
#include "image16.h"
#include "imagefloat.h"

const char rtengine::sImage8[] =     "Image8";
const char rtengine::sImage16[] =    "Image16";
const char rtengine::sImagefloat[] = "Imagefloat";
int rtengine::getCoarseBitMask( const procparams::CoarseTransformParams &coarse)
{
    int tr = TR_NONE;

    if (coarse.rotate == 90) {
        tr |= TR_R90;
    } else if (coarse.rotate == 180) {
        tr |= TR_R180;
    } else if (coarse.rotate == 270) {
        tr |= TR_R270;
    }

    if (coarse.hflip) {
        tr |= TR_HFLIP;
    }

    if (coarse.vflip) {
        tr |= TR_VFLIP;
    }

    return tr;
}

namespace rtengine
{

namespace
{

bool readScanlines (TIFF* tiff, guint8* data, const int count, const int sample, const gsize rowStride)
{
    for (int row = 0; row < count; ++row) {
        if (TIFFReadScanline (tiff, data, row, sample) < 0) {
            return false;
        }
        data += rowStride;
    }
    return true;
}

bool writeScanlines (TIFF* tiff, guint8* data, const int count, const int sample, const gsize rowStride)
{
    for (int row = 0; row < count; ++row) {
        if (TIFFWriteScanline (tiff, data, row, sample) < 0) {
            return false;
        }
        data += rowStride;
    }
    return true;
}

Image8* readImage8 (TIFF* tiff, const uint32 length, const uint32 width)
{
    std::unique_ptr< Image8 > image (new Image8 (width, length));

    const auto data = reinterpret_cast< guint8* > (image->data);
    const auto rowStride = 3 * sizeof (unsigned char) * width;

    const auto ok = readScanlines (tiff, data, length, 0, rowStride);

    return ok ? image.release () : nullptr;
}

Image16* readImage16 (TIFF* tiff, const uint32 length, const uint32 width)
{
    std::unique_ptr< Image16 > image (new Image16 (width, length));

    auto data = reinterpret_cast< guint8* > (image->data);
    const auto rowStride = image->getRowStride ();
    const auto planeStride = image->getPlaneStride ();

    auto ok = readScanlines (tiff, data, length, 0, rowStride);
    data += planeStride;
    ok = ok && readScanlines (tiff, data, length, 1, rowStride);
    data += planeStride;
    ok = ok && readScanlines (tiff, data, length, 2, rowStride);

    return ok ? image.release () : nullptr;
}

Imagefloat* readImagefloat (TIFF* tiff, const uint32 length, const uint32 width)
{
    std::unique_ptr< Imagefloat > image (new Imagefloat (width, length));

    auto data = reinterpret_cast< guint8* > (image->data);
    const auto rowStride = image->getRowStride ();
    const auto planeStride = image->getPlaneStride ();

    auto ok = readScanlines (tiff, data, length, 0, rowStride);
    data += planeStride;
    ok = ok && readScanlines (tiff, data, length, 1, rowStride);
    data += planeStride;
    ok = ok && readScanlines (tiff, data, length, 2, rowStride);

    return ok ? image.release() : nullptr;
}

}

ImageIO* IImage::readData (const char* fname)
{
    TIFF* tiff = TIFFOpen (fname, "r");
    if (!tiff) {
        return nullptr;
    }

    uint32 length, width;
    uint16 planarconfig, sampleformat;
    if (TIFFGetField (tiff, TIFFTAG_IMAGELENGTH, &length) == 0
     || TIFFGetField (tiff, TIFFTAG_IMAGEWIDTH, &width) == 0
     || TIFFGetField (tiff, TIFFTAG_PLANARCONFIG, &planarconfig) == 0
     || TIFFGetField (tiff, TIFFTAG_SAMPLEFORMAT, &sampleformat) == 0) {
        return nullptr;
    }

    ImageIO* image = nullptr;

    if (planarconfig == PLANARCONFIG_CONTIG && sampleformat == SAMPLEFORMAT_UINT) {
        image = readImage8 (tiff, length, width);
    }
    else if (planarconfig == PLANARCONFIG_SEPARATE && sampleformat == SAMPLEFORMAT_UINT) {
        image = readImage16 (tiff, length, width);
    }
    else if (planarconfig == PLANARCONFIG_SEPARATE && sampleformat == SAMPLEFORMAT_IEEEFP) {
        image = readImagefloat (tiff, length, width);
    }

    TIFFClose (tiff);
    return image;
}

bool IImage8::writeData (const char* fname)
{
    TIFF* tiff = TIFFOpen (fname, "w");
    if (!tiff) {
        return false;
    }

    TIFFSetField (tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField (tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField (tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField (tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    TIFFSetField (tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField (tiff, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField (tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField (tiff, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    TIFFSetField (tiff, TIFFTAG_PREDICTOR, PREDICTOR_HORIZONTAL);

    auto data = reinterpret_cast< guint8* > (this->data);
    const auto rowStride = 3 * sizeof (unsigned char) * width;

    const auto ok = writeScanlines (tiff, data, height, 0, rowStride);

    TIFFClose (tiff);
    return ok;
}

bool IImage16::writeData (const char* fname)
{
    TIFF* tiff = TIFFOpen (fname, "w");
    if (!tiff) {
        return false;
    }

    TIFFSetField (tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField (tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField (tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField (tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_SEPARATE);
    TIFFSetField (tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    TIFFSetField (tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField (tiff, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField (tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField (tiff, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    TIFFSetField (tiff, TIFFTAG_PREDICTOR, PREDICTOR_HORIZONTAL);

    auto data = reinterpret_cast< guint8* > (this->data);
    const auto rowStride = getRowStride ();
    const auto planeStride = getPlaneStride ();

    auto ok = writeScanlines (tiff, data, height, 0, rowStride);
    data += planeStride;
    ok = ok && writeScanlines (tiff, data, height, 1, rowStride);
    data += planeStride;
    ok = ok && writeScanlines (tiff, data, height, 2, rowStride);

    TIFFClose (tiff);
    return ok;
}

bool IImagefloat::writeData (const char* fname)
{
    TIFF* tiff = TIFFOpen (fname, "w");
    if (!tiff) {
        return false;
    }

    TIFFSetField (tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField (tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField (tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField (tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_SEPARATE);
    TIFFSetField (tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField (tiff, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField (tiff, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    TIFFSetField (tiff, TIFFTAG_PREDICTOR, PREDICTOR_FLOATINGPOINT);

    auto data = reinterpret_cast< guint8* > (this->data);
    const auto rowStride = getRowStride ();
    const auto planeStride = getPlaneStride ();

    auto ok = writeScanlines (tiff, data, height, 0, rowStride);
    data += planeStride;
    ok = ok && writeScanlines (tiff, data, height, 1, rowStride);
    data += planeStride;
    ok = ok && writeScanlines (tiff, data, height, 2, rowStride);

    TIFFClose (tiff);
    return ok;
}

}
