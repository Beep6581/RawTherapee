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
#ifndef __IMAGEAREA_H__
#define __IMAGEAREA_H__

#include <gtkmm.h>
#include <cropguilistener.h>
#include <imageareapanel.h>
#include <editenums.h>
#include <toolbar.h>
#include <previewhandler.h>
#include <imageareatoollistener.h>
#include <cropwindow.h>
#include <zoompanel.h>
#include <indclippedpanel.h>

class ImageAreaPanel;
class ImageArea : public Gtk::DrawingArea, public CropWindowListener {

    friend class ZoomPanel;

  protected:

    bool          showInfo;
    Glib::ustring infotext;
    Glib::RefPtr<Pango::Layout> ilayout;
    Glib::RefPtr<Pango::Layout> deglayout;
    Glib::RefPtr<Gdk::Pixbuf>   ipixbuf;
	bool showClippedH, showClippedS;

    ImageAreaPanel* parent;   
    CropWindow* mainCropWindow;
    std::list<CropWindow*> cropWins;
    PreviewHandler* previewHandler;
    rtengine::StagedImageProcessor* ipc;

    int lastClosedX, lastClosedY, lastClosedW, lastClosedH;

    bool        dirty;
    CropWindow* focusGrabber;
    CropGUIListener* cropgl;
	PointerMotionListener* pmlistener;
    ImageAreaToolListener* listener;

    CropWindow* getCropWindow (int x, int y);

  public:
        
    ZoomPanel* zoomPanel;
	IndicateClippedPanel* indClippedPanel;

    ImageArea (ImageAreaPanel* p);
    ~ImageArea ();
    
    void setImProcCoordinator (rtengine::StagedImageProcessor* ipc_);
    
    void getScrollImageSize (int& w, int& h);
    void getScrollPosition  (int& x, int& y);
    void setScrollPosition  (int x, int y);     // called by the imageareapanel when the scrollbars have been changed
    
    // enabling and setting text of info area
    void setInfoText (Glib::ustring text);
    void infoEnabled (bool e);

    // widget base events
    void on_realize ();
    bool on_expose_event        (GdkEventExpose* event);
    bool on_motion_notify_event (GdkEventMotion* event);
    bool on_button_press_event  (GdkEventButton* event);
    bool on_button_release_event (GdkEventButton* event);
    bool on_scroll_event        (GdkEventScroll* event);
    void on_resized             (Gtk::Allocation& req);
    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
    void updateScrollbars       ();

    void            setCropGUIListener       (CropGUIListener* l);
	void 			setPointerMotionListener (PointerMotionListener* pml);
    void            setImageAreaToolListener (ImageAreaToolListener* l) { listener = l; }
    void            setPreviewHandler        (PreviewHandler* ph);
    PreviewHandler* getPreviewHandler        ()                         { return previewHandler; }
    
    void grabFocus          (CropWindow* cw);
    void unGrabFocus        ();
    void addCropWindow      ();
    void cropWindowSelected (CropWindow* cw);
    void cropWindowClosed   (CropWindow* cw);
    ToolMode getToolMode    ();
    void setToolHand        ();
    void straightenReady    (double rotDeg);
    void spotWBSelected     (int x, int y);
    int  getSpotWBRectSize  ();
    void redraw             ();
    
    void zoomFit     ();
    double getZoom   ();  
    void   setZoom   (double zoom);  
    
    // cropwindowlistener interface
    void cropPositionChanged   (CropWindow* cw);
    void cropWindowSizeChanged (CropWindow* cw);  
    void cropZoomChanged       (CropWindow* cw);
    void initialImageArrived   (CropWindow* cw) ;

    CropWindow* getMainCropWindow () { return mainCropWindow; }
};



#endif
