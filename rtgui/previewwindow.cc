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
#include <previewwindow.h>
#include <guiutils.h>
#include <imagearea.h>

PreviewWindow::PreviewWindow () : previewHandler(NULL), mainCropWin(NULL), isMoving(false) {

    rconn = signal_size_allocate().connect( sigc::mem_fun(*this, &PreviewWindow::on_resized) );
}

PreviewWindow::~PreviewWindow () {

    delete cCropMoving;
#ifndef WIN32
    delete cNormal;
#endif
}

void PreviewWindow::on_realize () {

	Gtk::DrawingArea::on_realize ();
	add_events(Gdk::EXPOSURE_MASK | Gdk::POINTER_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::SCROLL_MASK);
	cCropMoving = new Gdk::Cursor (Gdk::FLEUR);
#ifdef _WIN32
	cNormal = new Gdk::Cursor (Gdk::LAST_CURSOR);
#else
	cNormal = new Gdk::Cursor (Gdk::ARROW);
#endif
}

void PreviewWindow::getObservedFrameArea (int& x, int& y, int& w, int& h) {

    if (mainCropWin) {
        int cropX, cropY, cropW, cropH;
        mainCropWin->getCropRectangle (cropX, cropY, cropW, cropH);
        // translate it to screen coordinates
        x = imgX + cropX*zoom;
        y = imgY + cropY*zoom;
        w = cropW * zoom;
        h = cropH * zoom;
    }
}

void PreviewWindow::updatePreviewImage () {

    int W = get_width(), H = get_height();
    backBuffer = Gdk::Pixmap::create (get_window(), W, H, -1);
    backBuffer->draw_rectangle (get_style()->get_base_gc(Gtk::STATE_NORMAL), true, 0, 0, W, H);
    if (previewHandler) {
        Glib::RefPtr<Gdk::Pixbuf> resPixbuf = previewHandler->getRoughImage (W, H, zoom);
        if (resPixbuf) {
            imgW = resPixbuf->get_width();
            imgH = resPixbuf->get_height();
            imgX = (W-imgW)/2;
            imgY = (H-imgH)/2;
            backBuffer->draw_pixbuf (get_style()->get_base_gc(Gtk::STATE_NORMAL), resPixbuf, 0, 0, imgX, imgY, -1, -1, Gdk::RGB_DITHER_NONE, 0, 0);
            Cairo::RefPtr<Cairo::Context> cr = backBuffer->create_cairo_context();
            if (previewHandler->getCropParams().enabled)
                drawCrop (cr, imgX, imgY, imgW, imgH, 0, 0, zoom, previewHandler->getCropParams());
        }
    }
}

void PreviewWindow::setPreviewHandler (PreviewHandler* ph) {

    previewHandler = ph;
    previewHandler->addPreviewImageListener (this);
}

void PreviewWindow::on_resized (Gtk::Allocation& req) {

    updatePreviewImage ();
    queue_draw ();
}

bool PreviewWindow::on_expose_event (GdkEventExpose* event) {

    if (backBuffer) {
        Glib::RefPtr<Gdk::Window> window = get_window();

        int bufferW, bufferH;
        backBuffer->get_size (bufferW, bufferH);

        if (!mainCropWin) {
            mainCropWin = imageArea->getMainCropWindow ();
            if (mainCropWin)
                mainCropWin->addCropWindowListener (this);
        }

        if (get_width()!=bufferW && get_height()!=bufferH)
            updatePreviewImage ();
            
        window->draw_drawable (get_style()->get_base_gc(Gtk::STATE_NORMAL), backBuffer, 0, 0, 0, 0, -1, -1);

        if (mainCropWin) {
            Cairo::RefPtr<Cairo::Context> cr = get_window()->create_cairo_context();
            int x, y, w, h;
            getObservedFrameArea (x, y, w, h);
            cr->set_source_rgb (1.0, 1.0, 1.0);
            cr->set_line_width (3);
            cr->rectangle (x-1.5, y-1.5, w+2, h+2);
            cr->stroke ();
            cr->set_source_rgb (1.0, 0.0, 0.0);
            cr->set_line_width (1);
            cr->rectangle (x-1.5, y-1.5, w+2, h+2);
            cr->stroke ();
        }
    }
    return true;
}
        
void PreviewWindow::previewImageChanged () {

    updatePreviewImage ();
    queue_draw ();
}

void PreviewWindow::setImageArea (ImageArea* ia) {

    imageArea = ia;
    mainCropWin = ia->getMainCropWindow ();
    if (mainCropWin)
        mainCropWin->addCropWindowListener (this);
}

void PreviewWindow::cropPositionChanged (CropWindow* w) {

    queue_draw ();
}

void PreviewWindow::cropWindowSizeChanged (CropWindow* w)  {

    queue_draw ();
}

void PreviewWindow::cropZoomChanged (CropWindow* w) {

    queue_draw ();
}

bool PreviewWindow::on_motion_notify_event (GdkEventMotion* event) {

    if (!mainCropWin)
		return true;

	int x, y, w, h;
	getObservedFrameArea (x, y, w, h);
	bool inside = event->x > x-6 && event->x < x+w-1+6 && event->y > y-6 && event->y < y+h-1+6;
	bool moreInside = event->x > x+6 && event->x < x+w-1-6 && event->y > y+6 && event->y < y+h-1-6;    

	if (isMoving) 
        mainCropWin->remoteMove ((event->x - press_x)/zoom, (event->y - press_y)/zoom);
	else if (inside && !moreInside)
        get_window()->set_cursor (*cCropMoving);
	else
        get_window()->set_cursor (*cNormal);

}

bool PreviewWindow::on_button_press_event (GdkEventButton* event) {

    if (!mainCropWin)
		return true;

	int x, y, w, h;
	getObservedFrameArea (x, y, w, h);
	bool inside = event->x > x-6 && event->x < x+w-1+6 && event->y > y-6 && event->y < y+h-1+6;
	bool moreInside = event->x > x+6 && event->x < x+w-1-6 && event->y > y+6 && event->y < y+h-1-6;    

	if (!isMoving) {
		isMoving = true;
        if (!inside || moreInside) {
            mainCropWin->remoteMove ((event->x - (x+w/2))/zoom, (event->y - (y+h/2))/zoom);
    		press_x = x+w/2;
    		press_y = y+h/2;
    	}
    	else {
    		press_x = event->x;
    		press_y = event->y;
    	}
        get_window()->set_cursor (*cCropMoving);
	}
}

bool PreviewWindow::on_button_release_event (GdkEventButton* event) {

    if (!mainCropWin)
		return true;

	if (isMoving) {
		isMoving = false;
		get_window()->set_cursor (*cNormal);
		mainCropWin->remoteMoveReady ();
	}
}
