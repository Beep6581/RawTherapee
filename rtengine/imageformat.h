/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Jean-Christophe Frisch <natureh.510@gmail.com>
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

namespace rtengine
{

//NB: Update the associated strings in languages files when updating the following enum
//    Look for "SAMPLEFORMAT_"
typedef enum IIO_Sample_Format {
    IIOSF_UNKNOWN        = 0,       // Unknown or Unsupported file type; Has to remain 0
    //IIOSF_SIGNED_INT         ,    // Not yet supported
    IIOSF_UNSIGNED_CHAR  = 1 << 0,
    IIOSF_UNSIGNED_SHORT = 1 << 1,
    IIOSF_LOGLUV24       = 1 << 2,
    IIOSF_LOGLUV32       = 1 << 3,
    IIOSF_FLOAT16        = 1 << 4,
    IIOSF_FLOAT24        = 1 << 5,
    IIOSF_FLOAT32        = 1 << 6
} IIOSampleFormat;

typedef enum IIO_Sample_Arrangement {
    IIOSA_UNKNOWN,       // Unknown or Unsupported file type
    IIOSA_CHUNKY,
    IIOSA_PLANAR
} IIOSampleArrangement;

typedef enum SensorType {
    ST_NONE,   // use this value if the image is already demosaiced (i.e. not a raw file)
    ST_BAYER,
    ST_FUJI_XTRANS,
    ST_FOVEON,
    //ST_FUJI_EXR
} eSensorType;

}
