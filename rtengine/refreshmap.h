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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <unordered_map>

#include "procevents.h"

// clang-format off

// Use M_VOID if you wish to update the proc params without updating the preview at all!
constexpr int M_VOID      = (1<<20);
// Use M_MINUPDATE if you wish to update the preview without modifying the image (think about it like a "refreshPreview")
// Must NOT be used with other event (i.e. will be used for MINUPDATE only)
constexpr int M_MINUPDATE = (1<<19);
// Force high quality
constexpr int M_HIGHQUAL  = (1<<18);

// Elementary functions that can be done to
// the preview image when an event occurs
constexpr int M_WB         = (1<<17);
constexpr int M_SPOT       = (1<<16);
constexpr int M_CSHARP     = (1<<15);
constexpr int M_MONITOR    = (1<<14);
constexpr int M_RETINEX    = (1<<13);
constexpr int M_CROP       = (1<<12);
constexpr int M_PREPROC    = (1<<11);
constexpr int M_RAW        = (1<<10);
constexpr int M_INIT       = (1<<9);
constexpr int M_LINDENOISE = (1<<8);
constexpr int M_HDR        = (1<<7);
constexpr int M_TRANSFORM  = (1<<6);
constexpr int M_BLURMAP    = (1<<5);
constexpr int M_AUTOEXP    = (1<<4);
constexpr int M_RGBCURVE   = (1<<3);
constexpr int M_LUMACURVE  = (1<<2);
constexpr int M_LUMINANCE  = (1<<1);
constexpr int M_COLOR      = (1<<0);

// Bitfield of functions to do to the preview image when an event occurs
// Use those or create new ones for your new events
constexpr int FIRST            = (M_PREPROC|M_RAW|M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR|M_MONITOR);  // without HIGHQUAL
constexpr int ALL              = (M_PREPROC|M_RAW|M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);  // without HIGHQUAL
constexpr int DARKFRAME        = (M_PREPROC|M_RAW|M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int FLATFIELD        = (M_PREPROC|M_RAW|M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int DEMOSAIC         =           (M_RAW|M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int WB               = (M_WB|           M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int COMPR            =                 (M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR|M_MONITOR);
constexpr int ALLNORAW         =                 (M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int CAPTURESHARPEN   =                 (M_INIT|M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR|M_CSHARP);
constexpr int HDR              =                        (M_SPOT|M_LINDENOISE|M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int SPOTADJUST       =                        (M_SPOT|             M_HDR|            M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int TRANSFORM        =                        (M_SPOT|             M_HDR|M_TRANSFORM|M_BLURMAP|M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int AUTOEXP          =                        (M_SPOT|             M_HDR|                      M_AUTOEXP|M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int RGBCURVE         =                                                                                  (M_RGBCURVE|M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int LUMINANCECURVE   =                                                                                             (M_LUMACURVE|M_LUMINANCE|M_COLOR);
constexpr int SHARPENING       =                                                                                                         (M_LUMINANCE|M_COLOR);
constexpr int IMPULSEDENOISE   =                                                                                                         (M_LUMINANCE|M_COLOR);
constexpr int DEFRINGE         =                                                                                                         (M_LUMINANCE|M_COLOR);
constexpr int DIRPYRDENOISE    =                                                                                                         (M_LUMINANCE|M_COLOR);
constexpr int DIRPYREQUALIZER  =                                                                                                         (M_LUMINANCE|M_COLOR);
constexpr int GAMMA            =  M_VOID;
constexpr int CROP             =  M_CROP;
constexpr int RESIZE           =  M_VOID;
constexpr int EXIF             =  M_VOID;
constexpr int IPTC             =  M_VOID;
constexpr int MINUPDATE        =  M_MINUPDATE;
constexpr int RETINEX          = (M_RETINEX|ALLNORAW);
constexpr int MONITORTRANSFORM =  M_MONITOR;
constexpr int OUTPUTPROFILE    =  M_MONITOR;

// clang-format on

extern int refreshmap[];

namespace rtengine
{

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
