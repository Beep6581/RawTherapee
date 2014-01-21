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
#include <thumbbrowserentry.h>

FileBrowserEntry::FileBrowserEntry (Thumbnail* thm, const Glib::ustring& fname) 
    : ThumbBrowserEntryBase (fname), thumbnail(thm) {
    
    previewOwner = false;
    italicstyle = thumbnail->getType() != FT_Raw;
    datetimeline = thumbnail->getDateTimeString ();
    exifline = thumbnail->getExifString ();
}

void ThumbBrowserEntry::obtainThumbnailSize () {

    if (thumbnail)
        thumbnail->getThumbnailSize (prew, preh);
}
Glib::RefPtr<Gdk::Pixbuf> ThumbBrowserEntry::editedIcon;
Glib::RefPtr<Gdk::Pixbuf> ThumbBrowserEntry::recentlySavedIcon;
Glib::RefPtr<Gdk::Pixbuf> ThumbBrowserEntry::enqueuedIcon;
std::vector<Glib::RefPtr<Gdk::Pixbuf> > ThumbBrowserEntry::getIconsOnImageArea () {

    std::vector<Glib::RefPtr<Gdk::Pixbuf> > ret;
    
    if (!thumbnail)
        return ret;

    if (thumbnail->hasProcParams() && editedIcon)
        ret.push_back (editedIcon);
    if (thumbnail->isRecentlySaved() && recentlySavedIcon)
        ret.push_back (recentlySavedIcon);
    if (thumbnail->isEnqueued () && enqueuedIcon)
        ret.push_back (enqueuedIcon);

   return ret;
}

ThumbnailButtonSet* ThumbBrowserEntry::getThumbButtonSet () {

    return (ThumbnailButtonSet*)buttonSet;
}
