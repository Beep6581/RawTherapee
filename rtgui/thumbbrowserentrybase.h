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
#ifndef _THUMBNAILBROWSERENTRYBASE_
#define _THUMBNAILBROWSERENTRYBASE_

#include <gtkmm.h>
#include "lwbuttonset.h"
#include "thumbnail.h"
#include "threadutils.h"

class ThumbBrowserBase;
class ThumbBrowserEntryBase {

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

    guint8* preview;  // holds the preview image. used in updateBackBuffer. TODO Olli: Make a cache to reduce mem significantly

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
    
    int redrawRequests;
    
    ThumbBrowserBase* parent;
    
    Glib::RefPtr<Gdk::Pixmap> backBuffer;
    bool bbSelected, bbFramed;
    guint8* bbPreview;
    std::vector<Glib::RefPtr<Gdk::Pixbuf> > bbIcons;
    
    void drawFrame (Cairo::RefPtr<Cairo::Context> cr, const Gdk::Color& bg, const Gdk::Color& fg);
    void getTextSizes (int& w, int& h);
    
    // called during updateBackBuffer for custom overlays
    virtual void customBackBufferUpdate (Cairo::RefPtr<Cairo::Context> c) {}
    
  public:
  
    Thumbnail* thumbnail;

// thumbnail preview properties:
    Glib::ustring filename;
    Glib::ustring shortname;
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
  
    ThumbBrowserEntryBase   (const Glib::ustring& fname);   
    virtual ~ThumbBrowserEntryBase  ();
    
    void setParent (ThumbBrowserBase* l) { parent = l; }

    void updateBackBuffer   ();
    void resize             (int h);  
    virtual void draw       ();
    
    void addButtonSet       (LWButtonSet* bs);
    int getMinimalHeight    () { return height; }
    int getMinimalWidth     () { return width; }

    int getEffectiveWidth   () const { return exp_width; }
    int getEffectiveHeight  () const { return exp_height; }
    int getPreviewHeight    () const { return preh; }
    int getStartX           () const { return startx; }
    int getStartY           () const { return starty; }
    int getX                () const { return ofsX+startx; }
    int getY                () const { return ofsY+starty; }

    bool inside             (int x, int y);
    bool insideWindow       (int x, int y, int w, int h);
    void setPosition        (int x, int y, int w, int h);
    void setOffset (int x, int y);

    bool operator< (ThumbBrowserEntryBase& other) { return shortname.casefold()>other.shortname.casefold(); } 
    
    virtual void refreshThumbnailImage () {}
    virtual void refreshQuickThumbnailImage () {}
    virtual void calcThumbnailSize () {}

    virtual void drawProgressBar (Glib::RefPtr<Gdk::Window> win, Glib::RefPtr<Gdk::GC> gc, const Gdk::Color& foregr, const Gdk::Color& backgr, int x, int w, int y, int h) {}

    virtual std::vector<Glib::RefPtr<Gdk::Pixbuf> > getIconsOnImageArea () { std::vector<Glib::RefPtr<Gdk::Pixbuf> >  r; return r; }
    virtual void getIconSize (int& w, int& h) { w=0; h=0; }

    virtual bool    motionNotify  (int x, int y);
    virtual bool    pressNotify   (int button, int type, int bstate, int x, int y);
    virtual bool    releaseNotify (int button, int type, int bstate, int x, int y);
    virtual Glib::ustring getToolTip (int x, int y);
};

#endif
