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

#include <giomm/inputstream.h>
#include <giomm/outputstream.h>

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

void readScanlines (const Glib::RefPtr< Gio::InputStream >& stream, guint8* data, const int count, const gsize rowSize, const gsize rowStride)
{
    for (int index = 0; index < count; ++index)
    {
        gsize bytesRead;
        stream->read_all (data, rowSize, bytesRead);
        assert (bytesRead == rowSize);

        data += rowStride;
    }
}

void writeScanlines (const Glib::RefPtr< Gio::OutputStream >& stream, const guint8* data, const int count, const gsize rowSize, const gsize rowStride)
{
    for (int index = 0; index < count; ++index)
    {
        gsize bytesWritten;
        stream->write_all (data, rowSize, bytesWritten);
        assert (bytesWritten == rowSize);

        data += rowStride;
    }
}

}
