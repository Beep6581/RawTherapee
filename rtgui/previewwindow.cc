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
#include "previewwindow.h"
#include "guiutils.h"
#include "imagearea.h"
#include "cursormanager.h"

PreviewWindow::PreviewWindow () : previewHandler(NULL), mainCropWin(NULL), imageArea(NULL), imgX(0), imgY(0), imgW(0), imgH(0),
    zoom(0.0), isMoving(false), needsUpdate(false), cursor_type(CSUndefined)
{
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    style->add_class(GTK_STYLE_CLASS_BACKGROUND);
    style->add_class(GTK_STYLE_CLASS_FLAT);
    rconn = signal_size_allocate().connect( sigc::mem_fun(*this, &PreviewWindow::on_resized) );
}

void PreviewWindow::on_realize ()
{

    Gtk::DrawingArea::on_realize ();
    add_events(Gdk::POINTER_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::SCROLL_MASK);
}

void PreviewWindow::getObservedFrameArea (int& x, int& y, int& w, int& h)
{

    if (mainCropWin) {
        int cropX, cropY, cropW, cropH;
        mainCropWin->getCropRectangle (cropX, cropY, cropW, cropH);
        // translate it to screen coordinates
        x = imgX + round(cropX * zoom);
        y = imgY + round(cropY * zoom);
        w = round(cropW * zoom);
        h = round(cropH * zoom);
    }
}

void PreviewWindow::updatePreviewImage ()
{

    int W = get_width(), H = get_height();
    Glib::RefPtr<Gdk::Window> wind = get_window();

    if( ! wind ) {
        needsUpdate = true;
        return;
    }

    backBuffer = Cairo::RefPtr<BackBuffer> ( new BackBuffer(W, H, wind) );
    Cairo::RefPtr<Cairo::ImageSurface> surface = backBuffer->getSurface();
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    Cairo::RefPtr<Cairo::Context> cc = Cairo::Context::create(surface);
    style->render_background(cc, 0, 0, W, H);
    Gdk::RGBA c = style->get_background_color(Gtk::STATE_FLAG_NORMAL);
    cc->set_antialias(Cairo::ANTIALIAS_NONE);
    cc->set_line_join(Cairo::LINE_JOIN_MITER);

    if (previewHandler) {
        Glib::RefPtr<Gdk::Pixbuf> resPixbuf = previewHandler->getRoughImage (W, H, zoom);

        if (resPixbuf) {
            imgW = resPixbuf->get_width();
            imgH = resPixbuf->get_height();
            imgX = (W - imgW) / 2;
            imgY = (H - imgH) / 2;
            Gdk::Cairo::set_source_pixbuf(cc, resPixbuf, imgX, imgY);
            cc->rectangle(imgX, imgY, imgW, imgH);
            cc->fill();

            if (previewHandler->getCropParams().enabled) {
                drawCrop (cc, imgX, imgY, imgW, imgH, 0, 0, zoom, previewHandler->getCropParams(), true, false);
            }
        }
    }
}

void PreviewWindow::setPreviewHandler (PreviewHandler* ph)
{

    previewHandler = ph;

    if (previewHandler) {
        previewHandler->addPreviewImageListener (this);
    }
}

void PreviewWindow::on_resized (Gtk::Allocation& req)
{

    updatePreviewImage ();
    queue_draw ();
}

bool PreviewWindow::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    if (backBuffer) {
        int bufferW, bufferH;
        bufferW = backBuffer->getWidth();
        bufferH = backBuffer->getHeight();

        if (!mainCropWin && imageArea) {
            mainCropWin = imageArea->getMainCropWindow ();

            if (mainCropWin) {
                mainCropWin->addCropWindowListener (this);
            }
        }

        if ((get_width() != bufferW && get_height() != bufferH) || needsUpdate) {
            needsUpdate = false;
            updatePreviewImage ();
        }

        backBuffer->copySurface(cr, NULL);

        if (mainCropWin && zoom > 0.0) {
            if(mainCropWin->getZoom() > mainCropWin->cropHandler.getFitZoom()) {
                int x, y, w, h;
                getObservedFrameArea (x, y, w, h);
                double rectX = x + 0.5;
                double rectY = y + 0.5;
                double rectW = std::min(w, (int)(imgW - (x - imgX) - 1));
                double rectH = std::min(h, (int)(imgH - (y - imgY) - 1));

                // draw a black "shadow" line
                cr->set_source_rgba (0.0, 0.0, 0.0, 0.65);
                cr->set_line_width (1.);
                cr->set_line_join(Cairo::LINE_JOIN_MITER);
                cr->rectangle (rectX + 1., rectY + 1, rectW, rectH);
                cr->stroke ();

                // draw a "frame" line. Color of frame line can be set in preferences
                cr->set_source_rgba(options.navGuideBrush[0], options.navGuideBrush[1], options.navGuideBrush[2], options.navGuideBrush[3]); //( 1.0, 1.0, 1.0, 1.0);
                cr->rectangle (rectX, rectY, rectW, rectH);
                cr->stroke ();
            }
        }
    }

    return true;
}

void PreviewWindow::previewImageChanged ()
{

    updatePreviewImage ();
    queue_draw ();
}

void PreviewWindow::setImageArea (ImageArea* ia)
{

    imageArea = ia;
    mainCropWin = ia->getMainCropWindow ();

    if (mainCropWin) {
        mainCropWin->addCropWindowListener (this);
    }
}

void PreviewWindow::cropPositionChanged (CropWindow* w)
{

    queue_draw ();
}

void PreviewWindow::cropWindowSizeChanged (CropWindow* w)
{

    queue_draw ();
}

void PreviewWindow::cropZoomChanged (CropWindow* w)
{

    queue_draw ();
}

bool PreviewWindow::on_motion_notify_event (GdkEventMotion* event)
{

    if (!mainCropWin) {
        return true;
    }


    if(mainCropWin->getZoom() > mainCropWin->cropHandler.getFitZoom()) {
        int x, y, w, h;
        getObservedFrameArea (x, y, w, h);
        bool inside = event->x > x - 6 && event->x < x + w - 1 + 6 && event->y > y - 6 && event->y < y + h - 1 + 6;
        bool moreInside = event->x > x + 6 && event->x < x + w - 1 - 6 && event->y > y + 6 && event->y < y + h - 1 - 6;

        CursorShape newType = cursor_type;

        if (isMoving) {
            mainCropWin->remoteMove ((event->x - press_x) / zoom, (event->y - press_y) / zoom);
        } else if (inside && !moreInside) {
            newType = CSClosedHand;
        } else {
            newType = CSArrow;
        }

        if (newType != cursor_type) {
            cursor_type = newType;
            CursorManager::setWidgetCursor(get_window(), cursor_type);
        }
    }

    return true;
}

bool PreviewWindow::on_button_press_event (GdkEventButton* event)
{

    if (!mainCropWin) {
        return true;
    }

    if(mainCropWin->getZoom() > mainCropWin->cropHandler.getFitZoom()) {
        int x, y, w, h;
        getObservedFrameArea (x, y, w, h);
        bool inside = event->x > x - 6 && event->x < x + w - 1 + 6 && event->y > y - 6 && event->y < y + h - 1 + 6;
        bool moreInside = event->x > x + 6 && event->x < x + w - 1 - 6 && event->y > y + 6 && event->y < y + h - 1 - 6;

        if (!isMoving) {
            isMoving = true;

            if (!inside || moreInside) {
                mainCropWin->remoteMove ((event->x - (x + w / 2)) / zoom, (event->y - (y + h / 2)) / zoom);
                press_x = x + w / 2;
                press_y = y + h / 2;
            } else {
                press_x = event->x;
                press_y = event->y;
            }

            if (cursor_type != CSClosedHand) {
                cursor_type = CSClosedHand;
                CursorManager::setWidgetCursor(get_window(), cursor_type);
            }
        }
    }

    return true;
}

bool PreviewWindow::on_button_release_event (GdkEventButton* event)
{

    if (!mainCropWin) {
        return true;
    }

    if (isMoving) {
        isMoving = false;

        if (cursor_type != CSArrow) {
            cursor_type = CSArrow;
            CursorManager::setWidgetCursor(get_window(), cursor_type);
        }

        mainCropWin->remoteMoveReady ();
    }

    return true;
}
