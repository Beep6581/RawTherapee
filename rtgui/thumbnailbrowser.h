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
#ifndef _THUMBNAILBROWSER_
#define _THUMBNAILBROWSER_

#include <gtkmm.h>
#include "thumbnail.h"
#include "filecatalog.h"

class ThumbBrowserEntry
{

public:
// set by arrangeFiles():
    int width;      // minimal width
    int height;     // minimal height
    int exp_width;  // ararnged width
    int startx;     // x coord. in the widget
    int starty;     // y coord. in the widget
// thumbnail preview properties:
    int prew;       // width of the thumbnail
    int preh;       // height of the thumbnail
    guint8* preview;
// file and directory attributes:
    Glib::ustring filename;
    Glib::ustring shortname;
    Glib::ustring dirname;
// the associated thumbnail instance:
    Thumbnail* thumbnail;

    ThumbBrowserEntry (Thumbnail* thm, Glib::ustring fname, Glib::ustring sname, Glib::ustring dname, int h)
        : thumbnail(thm), filename(fname), shortname(sname), dirname(dname), preh(h)
    {
        preview = thumbnail ? (guint8*) thumbnail->getThumbnailImage (prew, preh) : NULL;
    }

    bool operator< (FileDescr& other)
    {
        return shortname > other.shortname;
    }
};

class ThumbBrowser  : public Gtk::DrawingArea
{

protected:
    int dx, dy, w, h;

    Glib::RefPtr<Gdk::GC> gc_;
    Gdk::Color black;
    Gdk::Color white;
    Gdk::Color blue;
    Gdk::Color bluew;

    std::vector<ThumbBrowserEntry*> fd;
    std::vector<ThumbBrowserEntry*> selected;

    int rowHeight;
    int numOfRows;

    ThumbBrowserListener* tbl;

    void arrangeFiles (int rows);

public:

    ThumbBrowser ();

    void addEntry (ThumbBrowserEntry* entry);
    void setThumbBrowserListener (ThumbBrowserListener* l)
    {
        tbl = l;
    }

    virtual void on_realize();
    virtual bool on_expose_event(GdkEventExpose* event);
    virtual bool on_button_press_event (GdkEventButton* event);
    virtual bool on_button_release_event (GdkEventButton* event);
    virtual void previewReady (FileDescr* fdn);

    void resized (Gtk::Allocation& req);
    void redraw ();
    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
};

#endif
