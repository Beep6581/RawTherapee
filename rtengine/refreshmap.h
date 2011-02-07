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
#ifndef __REFRESHMAP__
#define __REFRESHMAP__

#include <glibmm.h>

// Use M_VOID is you wish to update the proc params without updating the preview at all !
#define M_VOID       (1<<15)
// Use M_MINUPDATE if you you wish to update the preview without modifying the image (think about it like a "refreshPreview")
#define M_MINUPDATE  (1<<14)

// Elementary functions that can be done to
// the preview image when an event occurs
#define M_PREPROC    (1<<9)
#define M_RAW        (1<<8)
#define M_INIT       (1<<7)
#define M_TRANSFORM  (1<<6)
#define M_BLURMAP    (1<<5)
#define M_AUTOEXP    (1<<4)
#define M_RGBCURVE   (1<<3)
#define M_LUMACURVE  (1<<2)
#define M_LUMINANCE  (1<<1)
#define M_COLOR      (1<<0)

// Bitfield of functions to do to the preview image when an event occurs
// Use those or create new ones for your new events
#define FIRST           65535
#define ALL             65535
#define TRANSFORM       (M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define RETINEX         (M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define AUTOEXP         (M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define RGBCURVE        (M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define LUMINANCECURVE  (M_LUMACURVE|M_LUMINANCE)
#define SHARPENING       M_LUMINANCE
#define IMPULSEDENOISE   M_LUMINANCE
#define DEFRINGE         M_LUMINANCE
#define LUMADENOISE      M_LUMINANCE
#define WHITEBALANCE    (M_INIT|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define DEMOSAIC        (M_RAW|M_INIT|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define DARKFRAME       (M_PREPROC|M_RAW|M_INIT|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define FLATFIELD       (M_PREPROC|M_RAW|M_INIT|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define COLORBOOST       M_COLOR
#define COLORDENOISE     M_COLOR
#define DIRPYRDENOISE   (M_COLOR|M_LUMINANCE)
#define CROP             M_MINUPDATE
#define RESIZE           M_VOID
#define EXIF             M_VOID
#define IPTC             M_VOID
#define EQUALIZER       (M_COLOR|M_LUMINANCE)
#define DIRPYREQUALIZER	(M_COLOR|M_LUMINANCE)
#define NONE             0

extern int refreshmap[];
#endif    
