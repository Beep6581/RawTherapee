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
#ifndef _PREVIEWWINDOW_
#define _PREVIEWWINDOW_

#include <gtkmm.h>
#include <previewhandler.h>
#include <cropwindow.h>

class PreviewWindow : public Gtk::DrawingArea, public PreviewListener, public CropWindowListener {

    private:      
        Glib::RefPtr<Gdk::Pixmap> backBuffer;
        int oldwidth, oldheight;
        PreviewHandler* previewHandler;
        sigc::connection rconn;
        CropWindow* mainCropWin;
        ImageArea* imageArea;
        int imgX, imgY, imgW, imgH;
        double zoom;
		Gdk::Cursor* cCropMoving;
		Gdk::Cursor* cNormal;
		int press_x, press_y;
		bool isMoving;
		
        void updatePreviewImage     ();
        void getObservedFrameArea   (int& x, int& y, int& w, int& h);
        
    public:
        PreviewWindow ();
		~PreviewWindow ();
        
        void setPreviewHandler  (PreviewHandler* ph);
        void setImageArea       (ImageArea* ia);

		void on_realize 			();
        void on_resized         	(Gtk::Allocation& req);
        bool on_expose_event    	(GdkEventExpose* event);
		bool on_motion_notify_event (GdkEventMotion* event);
		bool on_button_press_event  (GdkEventButton* event);
		bool on_button_release_event(GdkEventButton* event);
        
        // PreviewListener interface
        void previewImageChanged ();
        
        // CropWindowListener interface
        void cropPositionChanged   (CropWindow* w);
        void cropWindowSizeChanged (CropWindow* w);
        void cropZoomChanged       (CropWindow* w);
};

#endif
