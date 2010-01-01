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

class ThumbBrowserBase  : public Gtk::DrawingArea {  

  public:

    enum Arrangement {TB_Horizontal, TB_Vertical};

  protected:

    Glib::RefPtr<Gdk::GC> gc_;
  
    std::vector<ThumbBrowserEntryBase*> fd;
    std::vector<ThumbBrowserEntryBase*> selected;
    ThumbBrowserEntryBase* lastClicked;
    
    int previewHeight;

    Glib::ustring fileInEditor;
    Arrangement arrangement;

    void arrangeFiles ();

  public:
   
    ThumbBrowserBase ();

    void setPreviewHeight (int h);
    const std::vector<ThumbBrowserEntryBase*>& getEntries () { return fd; } 
    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
    void redraw ();   // arrange files and draw area
    void refreshAll (); // refresh thumbnail sizes, re-generate thumbnail images, arrange and draw
    void resized (Gtk::Allocation& req);
    void setOpenedFileName (const Glib::ustring& fname);
    
    virtual void on_realize();
    virtual bool on_expose_event(GdkEventExpose* event);
    virtual bool on_button_press_event (GdkEventButton* event);
    virtual bool on_button_release_event (GdkEventButton* event);
    virtual bool on_motion_notify_event (GdkEventMotion* event);
  
    void setArrangement (Arrangement a);
    virtual bool checkFilter (ThumbBrowserEntryBase* entry) { return true; }
    virtual void rightClicked (ThumbBrowserEntryBase* entry) {}
    virtual void doubleClicked (ThumbBrowserEntryBase* entry) {}
    
};

#endif
