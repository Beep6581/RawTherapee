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
#ifndef _FILEBROWSERENTRY_
#define _FILEBROWSERENTRY_

#include <gtkmm.h>
#include "thumbbrowserentrybase.h"
#include "thumbnail.h"
#include "filethumbnailbuttonset.h"
#include "thumbnaillistener.h"
#include "thumbimageupdater.h"
#include "imageareatoollistener.h"
#include "editenums.h"
#include "../rtengine/rtengine.h"
#include "crophandler.h"


class FileBrowserEntry;
struct FileBrowserEntryIdleHelper {
    FileBrowserEntry* fbentry;
    bool destroyed;
    int pending;
};

class FileThumbnailButtonSet;
class FileBrowserEntry : public ThumbBrowserEntryBase,
    public ThumbnailListener,
    public ThumbImageUpdateListener
{

    double scale;
    static bool iconsLoaded;
    bool wasInside;
    ImageAreaToolListener* iatlistener;
    int press_x, press_y, action_x, action_y;
    double rot_deg;
    bool landscape;
    rtengine::procparams::CropParams cropParams;
    CropGUIListener* cropgl;
    FileBrowserEntryIdleHelper* feih;

    ImgEditState state;

    bool onArea (CursorArea a, int x, int y);
    void updateCursor (int x, int y);
    void drawStraightenGuide (Cairo::RefPtr<Cairo::Context> c);
    void customBackBufferUpdate (Cairo::RefPtr<Cairo::Context> c);

public:

    static Glib::RefPtr<Gdk::Pixbuf> editedIcon;
    static Glib::RefPtr<Gdk::Pixbuf> recentlySavedIcon;
    static Glib::RefPtr<Gdk::Pixbuf> enqueuedIcon;

    FileBrowserEntry (Thumbnail* thm, const Glib::ustring& fname);
    ~FileBrowserEntry ();
    void draw (Cairo::RefPtr<Cairo::Context> cc);

    void setImageAreaToolListener (ImageAreaToolListener* l)
    {
        iatlistener = l;
    }

    FileThumbnailButtonSet* getThumbButtonSet ();

    void refreshThumbnailImage ();
    void refreshQuickThumbnailImage ();
    void calcThumbnailSize ();

    virtual std::vector<Glib::RefPtr<Gdk::Pixbuf> > getIconsOnImageArea ();
    virtual void getIconSize (int& w, int& h);

    // thumbnaillistener interface
    void procParamsChanged (Thumbnail* thm, int whoChangedIt);
    // thumbimageupdatelistener interface
    void updateImage (rtengine::IImage8* img, double scale, rtengine::procparams::CropParams cropParams);
    void _updateImage (rtengine::IImage8* img, double scale, rtengine::procparams::CropParams cropParams); // inside gtk thread

    virtual bool    motionNotify  (int x, int y);
    virtual bool    pressNotify   (int button, int type, int bstate, int x, int y);
    virtual bool    releaseNotify (int button, int type, int bstate, int x, int y);
};

#endif
