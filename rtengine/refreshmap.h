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

#include <unordered_map>
#include "procevents.h"

// Use M_MODE_RAWCROPADJUST to update the rendering for On Preview Adjustment of the RawCrop tool
#define M_MODE_RAWCROPADJUST  (1<<19)

// Use M_VOID if you wish to update the proc params without updating the preview at all !
#define M_VOID       (1<<18)
// Use M_MINUPDATE if you wish to update the preview without modifying the image (think about it like a "refreshPreview")
// Must NOT be used with other event (i.e. will be used for MINUPDATE only)
#define M_MINUPDATE  (1<<17)
// Force high quality
#define M_HIGHQUAL   (1<<16)

// Elementary functions that can be done to
// the preview image when an event occurs
#define M_RELOADRAW   (1<<15)
#define M_MONITOR     (1<<14)
#define M_RETINEX     (1<<13)
#define M_CROP        (1<<12)
#define M_PREPROC     (1<<11)
#define M_RAW         (1<<10)
#define M_INIT        (1<<9)
#define M_LINDENOISE  (1<<8)
#define M_HDR         (1<<7)
#define M_TRANSFORM   (1<<6)
#define M_BLURMAP     (1<<5)
#define M_AUTOEXP     (1<<4)
#define M_RGBCURVE    (1<<3)
#define M_LUMACURVE   (1<<2)
#define M_LUMINANCE   (1<<1)
#define M_COLOR       (1<<0)

// Bitfield of functions to do to the preview image when an event occurs
// Use those or create new ones for your new events
#define FIRST            (M_PREPROC|M_RAW|M_INIT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR|M_MONITOR)  // without HIGHQUAL
#define ALL              (M_PREPROC|M_RAW|M_INIT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)  // without HIGHQUAL
#define DARKFRAME        (M_PREPROC|M_RAW|M_INIT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define FLATFIELD        (M_PREPROC|M_RAW|M_INIT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define DEMOSAIC                   (M_RAW|M_INIT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define ALLNORAW                         (M_INIT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define HDR                                     (M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define TRANSFORM                                                  (M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define AUTOEXP                                                                          (M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define RGBCURVE                                                                                   (M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define LUMINANCECURVE                                                                                        (M_LUMACURVE|M_LUMINANCE|M_COLOR)
#define SHARPENING                                                                                                        (M_LUMINANCE|M_COLOR)
#define IMPULSEDENOISE                                                                                                    (M_LUMINANCE|M_COLOR)
#define DEFRINGE                                                                                                          (M_LUMINANCE|M_COLOR)
#define DIRPYRDENOISE                                                                                                     (M_LUMINANCE|M_COLOR)
#define DIRPYREQUALIZER                                                                                                   (M_LUMINANCE|M_COLOR)
#define GAMMA             M_MONITOR
#define RAWCROP           (M_RELOADRAW|ALL)
#define CROP              M_CROP
#define RESIZE            M_VOID
#define EXIF              M_VOID
#define IPTC              M_VOID
#define MINUPDATE         M_MINUPDATE
#define RETINEX          (M_RETINEX|ALLNORAW)
#define MONITORTRANSFORM  M_MONITOR
#define OUTPUTPROFILE     M_MONITOR

extern int refreshmap[];

namespace rtengine {

class RefreshMapper {
public:
    static RefreshMapper *getInstance();
    ProcEvent newEvent();
    void mapEvent(ProcEvent event, int action);
    int getAction(ProcEvent event) const;
    
private:
    RefreshMapper();

    int next_event_;
    std::unordered_map<int, int> actions_;
};

} // namespace rtengine

#endif
