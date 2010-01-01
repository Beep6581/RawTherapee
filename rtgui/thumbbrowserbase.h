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
#ifndef _THUMBNAILBROWSERBASE_
#define _THUMBNAILBROWSERBASE_

#include <gtkmm.h>
#include <thumbbrowserentrybase.h>
#include <set>

class ThumbBrowserBase  :  public Gtk::VBox {  

    class Internal : public Gtk::DrawingArea {
        
            Glib::RefPtr<Gdk::GC> gc_;
            int ofsX, ofsY;
            ThumbBrowserBase* parent;
            bool dirty;
        public:
            Internal ();
            void setParent (ThumbBrowserBase* p);
            void on_realize();
            bool on_expose_event(GdkEventExpose* event);
            bool on_button_press_event (GdkEventButton* event);
            bool on_button_release_event (GdkEventButton* event);
            bool on_motion_notify_event (GdkEventMotion* event);
            bool on_scroll_event (GdkEventScroll* event);
            bool on_key_press_event (GdkEventKey* event);
            bool on_query_tooltip (int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip);
            void setPosition (int x, int y);
            
            void setDirty () { dirty = true; }
            bool isDirty  () { return dirty; }
    };

  protected:

    Internal internal;
    Gtk::HScrollbar hscroll;
    Gtk::VScrollbar vscroll;
    
    int inW, inH;
    
    void resizeThumbnailArea (int w, int h);
    void internalAreaResized (Gtk::Allocation& req);
    void buttonPressed (int x, int y, int button, GdkEventType type, int state, int clx, int cly, int clw, int clh);
    
  public:

    enum Arrangement {TB_Horizontal, TB_Vertical};
    void configScrollBars ();
    void scrollChanged ();
    void scroll (int direction);

  protected:

    int eventTime;

    std::vector<ThumbBrowserEntryBase*> fd;
    std::vector<ThumbBrowserEntryBase*> selected;
    ThumbBrowserEntryBase* lastClicked;
    
    int previewHeight;

    Arrangement arrangement;

    std::set<Glib::ustring> editedFiles;

    void arrangeFiles ();
    void zoomChanged (bool zoomIn);

  public:
   
    ThumbBrowserBase ();

    void zoomIn ()  { zoomChanged (true); }
    void zoomOut () { zoomChanged (false); }
    
    const std::vector<ThumbBrowserEntryBase*>& getEntries () { return fd; } 
    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
    void redraw ();   // arrange files and draw area
    void refreshThumbImages (); // refresh thumbnail sizes, re-generate thumbnail images, arrange and draw
    void refreshEditedState (const std::set<Glib::ustring>& efiles);
     
    void initEntry (ThumbBrowserEntryBase* entry);

    void getScrollPosition (double& h, double& v);
    void setScrollPosition (double h, double v);
     
    void setArrangement (Arrangement a);
    virtual bool checkFilter (ThumbBrowserEntryBase* entry) { return true; }
    virtual void rightClicked (ThumbBrowserEntryBase* entry) {}
    virtual void doubleClicked (ThumbBrowserEntryBase* entry) {}
    virtual bool keyPressed (GdkEventKey* event) {}
    virtual void selectionChanged () {}
    
    virtual void redrawNeeded (ThumbBrowserEntryBase* entry) {}
    virtual void thumbRearrangementNeeded () {}

    Gtk::Widget* getDrawingArea () { return &internal; }
};

#endif
