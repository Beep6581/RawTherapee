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

#include <atomic>
#include <tuple>
#include <gtkmm.h>

#include "cursormanager.h"
#include "guiutils.h"
#include "lwbuttonset.h"
#include "threadutils.h"
#include "options.h"
#include "thumbnail.h"

#include "../rtengine/coord2d.h"

class Thumbnail;
class ThumbBrowserBase;
class ThumbBrowserEntryBase
{

public:
    enum eWithFilename {
        WFNAME_NONE,
        WFNAME_REDUCED,
        WFNAME_FULL
    };

protected:
    int fnlabw, fnlabh; // dimensions of the filename label
    int dtlabw, dtlabh; // dimensions of the date/time label
    int exlabw, exlabh; // dimensions of the exif label
    int prew;       // width of the thumbnail
    int preh;       // height of the thumbnail
    int prex;
    int prey;

    int upperMargin;
    int borderWidth;
    int textGap;
    int sideMargin;
    int lowerMargin;


    MyRWMutex lockRW;  // Locks access to all image thumb changing actions

    std::vector<guint8> preview;  // holds the preview image. used in updateBackBuffer.

    Glib::ustring dispname;

    LWButtonSet* buttonSet;

    int width;      // minimal width
    int height;     // minimal height
// set by arrangeFiles() of thumbbrowser
    int exp_width;  // arranged width (backbuffer dimensions)
    int exp_height; // arranged height
    int startx;     // x coord. in the widget
    int starty;     // y coord. in the widget

    int ofsX, ofsY; // offset due to the scrolling of the parent

    std::atomic<int> redrawRequests;

    ThumbBrowserBase* parent;
    ThumbBrowserEntryBase* original;

    Glib::RefPtr<BackBuffer> backBuffer;
    bool bbSelected, bbFramed;
    guint8* bbPreview;
    std::vector<Glib::RefPtr<Gdk::Pixbuf>> bbIcons;
    std::vector<Glib::RefPtr<Gdk::Pixbuf>> bbSpecificityIcons;
    CursorShape cursor_type;

    void drawFrame (Cairo::RefPtr<Cairo::Context> cr, const Gdk::RGBA& bg, const Gdk::RGBA& fg);
    void getTextSizes (int& w, int& h);

    // called during updateBackBuffer for custom overlays
    virtual void customBackBufferUpdate (Cairo::RefPtr<Cairo::Context> c) {}

private:
    const std::string collate_name;
    const std::string collate_exif;

public:

    Thumbnail* thumbnail;

// thumbnail preview properties:
    Glib::ustring filename;
    Glib::ustring exifline;
    Glib::ustring datetimeline;

// misc attributes
    bool selected;
    bool drawable;
    bool filtered;
    bool framed;
    bool processing;
    bool italicstyle;
    bool edited;
    bool recentlysaved;
    bool updatepriority;
    eWithFilename withFilename;

    explicit ThumbBrowserEntryBase (const Glib::ustring& fname, Thumbnail *thm);
    virtual ~ThumbBrowserEntryBase ();

    void setParent (ThumbBrowserBase* l)
    {
        parent = l;
    }

    void updateBackBuffer ();
    void resize (int h);
    virtual void draw (Cairo::RefPtr<Cairo::Context> cc);

    void addButtonSet (LWButtonSet* bs);
    int getMinimalHeight () const
    {
        return height;
    }
    int getMinimalWidth () const
    {
        return width;
    }

    int getEffectiveWidth () const
    {
        return exp_width;
    }
    int getEffectiveHeight () const
    {
        return exp_height;
    }
    int getPreviewHeight () const
    {
        return preh;
    }
    int getStartX () const
    {
        return startx;
    }
    int getStartY () const
    {
        return starty;
    }
    int getX () const
    {
        return ofsX + startx;
    }
    int getY () const
    {
        return ofsY + starty;
    }

    bool inside (int x, int y) const;
    rtengine::Coord2D getPosInImgSpace (int x, int y) const;
    bool insideWindow (int x, int y, int w, int h) const;
    void setPosition (int x, int y, int w, int h);
    void setOffset (int x, int y);

    bool compare (const ThumbBrowserEntryBase& other, Options::SortMethod method) const
    {
        int cmp = 0;
        switch (method){
        case Options::SORT_BY_NAME:
            return collate_name < other.collate_name;
        case Options::SORT_BY_DATE:
            cmp = thumbnail->getDateTime().compare(other.thumbnail->getDateTime());
            break;
        case Options::SORT_BY_EXIF:
            cmp = collate_exif.compare(other.collate_exif);
            break;
        case Options::SORT_BY_RANK:
            cmp = thumbnail->getRank() - other.thumbnail->getRank();
            break;
        case Options::SORT_BY_LABEL:
            cmp = thumbnail->getColorLabel() - other.thumbnail->getColorLabel();
            break;
        case Options::SORT_METHOD_COUNT: abort();
        }

        // Always fall back to sorting by name
        if (!cmp)
            cmp = collate_name.compare(other.collate_name);

        return cmp < 0;
    }

    virtual void refreshThumbnailImage () = 0;
    virtual void refreshQuickThumbnailImage () {}
    virtual void calcThumbnailSize () = 0;

    virtual void drawProgressBar (Glib::RefPtr<Gdk::Window> win, const Gdk::RGBA& foregr, const Gdk::RGBA& backgr, int x, int w, int y, int h) {}

    virtual std::vector<Glib::RefPtr<Gdk::Pixbuf>> getIconsOnImageArea ();
    virtual std::vector<Glib::RefPtr<Gdk::Pixbuf>> getSpecificityIconsOnImageArea ();
    virtual void getIconSize (int& w, int& h) const = 0;

    virtual bool motionNotify (int x, int y);
    virtual bool pressNotify (int button, int type, int bstate, int x, int y);
    virtual bool releaseNotify (int button, int type, int bstate, int x, int y);
    virtual std::tuple<Glib::ustring, bool> getToolTip (int x, int y) const;

    inline ThumbBrowserEntryBase* getOriginal() const
    {
        return original;
    }

    inline void setOriginal(ThumbBrowserEntryBase* original)
    {
        this->original = original;
    }

};
