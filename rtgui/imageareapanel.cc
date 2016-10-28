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
#include "imageareapanel.h"

ImageAreaPanel::ImageAreaPanel () : before(nullptr), after(nullptr)
{

    imageArea = new ImageArea (this);

    Gtk::HBox*  hb1   = Gtk::manage (new Gtk::HBox ());
    Gtk::Frame* frame = Gtk::manage (new Gtk::Frame ());

    frame->add (*imageArea);
    frame->set_shadow_type (Gtk::SHADOW_IN );
    hb1->pack_start (*frame, Gtk::PACK_EXPAND_WIDGET);

    pack_start (*hb1);
    frame->show ();
    imageArea->show ();
    hb1->show ();

}

ImageAreaPanel::~ImageAreaPanel ()
{

    delete imageArea;
}

void ImageAreaPanel::syncBeforeAfterViews ()
{

    if (before && this == after) {
        before->synchronize ();
    } else if (after && this == before) {
        after->synchronize ();
    }

    queue_draw ();
}

void ImageAreaPanel::setBeforeAfterViews (ImageAreaPanel* bef, ImageAreaPanel* aft)
{

    before = bef;
    after = aft;
    syncBeforeAfterViews ();
}

void ImageAreaPanel::zoomChanged ()
{

    if (after && this == before) {
        after->imageArea->setZoom (imageArea->getZoom ());
    } else if (before && this == after) {
        before->imageArea->setZoom (imageArea->getZoom ());
    }
}

void ImageAreaPanel::synchronize ()
{

    if (after && this == before) {
        int imgw, imgh, x, y;
        after->imageArea->getScrollImageSize (imgw, imgh);
        after->imageArea->getScrollPosition (x, y);

        if (imgw > 0 && imgh > 0) {
            imageArea->setScrollPosition (x, y);
            imageArea->queue_draw ();
        }
    } else if (before && this == after) {
        int imgw, imgh, x, y;
        before->imageArea->getScrollImageSize (imgw, imgh);
        before->imageArea->getScrollPosition (x, y);

        if (imgw > 0 && imgh > 0) {
            imageArea->setScrollPosition (x, y);
            imageArea->queue_draw ();
        }
    }

}

