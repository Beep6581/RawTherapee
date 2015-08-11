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
#ifndef _PREVIEWIMAGE_
#define _PREVIEWIMAGE_

#include <gtkmm.h>
#include "cairomm/cairomm.h"

namespace rtengine
{

/** @brief Get a quick preview image out of a raw or standard file
 *
 * This class reads the full size preview image (at least the biggest one available) from the raw file,
 * or the fast demosaiced version if no suitable embedded preview is found.
 *
 * For standard image, it simply read it with fast conversion for 32 bits images
 */
class PreviewImage
{

private:
    Cairo::RefPtr<Cairo::ImageSurface> previewImage;

public:
    typedef enum mode {
        PIM_EmbeddedPreviewOnly,  /// Get the embedded image only, fail if doesn't exist
        PIM_EmbeddedOrRaw,        /// Get the embedded image if it exist, or use the raw file otherwise
        PIM_ForceRaw              /// Get a preview of the raw file, even if an embedded image exist
    } PreviewImageMode;

    PreviewImage (const Glib::ustring &fname, const Glib::ustring &ext, const PreviewImageMode mode);

    Cairo::RefPtr<Cairo::ImageSurface> getImage();

};

}

#endif
