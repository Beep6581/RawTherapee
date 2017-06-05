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
#include "previewhandler.h"
#include "cropwindow.h"
#include "guiutils.h"
#include "cursormanager.h"

class PreviewWindow : public Gtk::DrawingArea, public PreviewListener, public CropWindowListener
{

private:
    Cairo::RefPtr<BackBuffer> backBuffer;
    PreviewHandler* previewHandler;
    sigc::connection rconn;
    CropWindow* mainCropWin;
    ImageArea* imageArea;
    int imgX, imgY, imgW, imgH;
    double zoom;
    int press_x, press_y;
    bool isMoving;
    bool needsUpdate;
    CursorShape cursor_type;

    void updatePreviewImage     ();
    void getObservedFrameArea   (int& x, int& y, int& w, int& h);

public:
    PreviewWindow ();

    void setPreviewHandler  (PreviewHandler* ph);
    void setImageArea       (ImageArea* ia);

    void on_realize             ();
    void on_resized             (Gtk::Allocation& req);
    bool on_draw                (const ::Cairo::RefPtr< Cairo::Context> &cr);
    bool on_motion_notify_event (GdkEventMotion* event);
    bool on_button_press_event  (GdkEventButton* event);
    bool on_button_release_event(GdkEventButton* event);
    Gtk::SizeRequestMode get_request_mode_vfunc () const;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const;
    void get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const;

    // PreviewListener interface
    void previewImageChanged ();

    // CropWindowListener interface
    void cropPositionChanged   (CropWindow* w);
    void cropWindowSizeChanged (CropWindow* w);
    void cropZoomChanged       (CropWindow* w);
};

#endif
