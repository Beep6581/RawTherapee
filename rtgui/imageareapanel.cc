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
#include <imageareapanel.h>

ImageAreaPanel::ImageAreaPanel () : before(NULL), after(NULL) {

    set_border_width (2);

    Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* hb2 = Gtk::manage (new Gtk::HBox ());
    hscroll = Gtk::manage (new Gtk::HScrollbar ());
    vscroll = Gtk::manage (new Gtk::VScrollbar ());
    imageArea = new ImageArea (this);
    Gtk::Frame* frame = Gtk::manage (new Gtk::Frame ());
    frame->add (*imageArea);
    frame->set_shadow_type (Gtk::SHADOW_IN );
    hb1->pack_start (*frame);
    hb1->pack_end (*vscroll, Gtk::PACK_SHRINK, 0);

    pack_start (*hb1);
    vscroll->show ();
    frame->show ();
    imageArea->show ();
    hb1->show ();

    Gtk::HBox* tmp = Gtk::manage (new Gtk::HBox ());
    hb2->pack_start (*hscroll);
    Gtk::Requisition vcr = vscroll->size_request ();
    tmp->set_size_request (vcr.width, vcr.width);
    hb2->pack_end (*tmp, Gtk::PACK_SHRINK, 0);

    pack_start (*hb2,Gtk::PACK_SHRINK, 0);
    hscroll->show ();
    tmp->show ();
    hb2->show ();

    hscroll->set_update_policy (Gtk::UPDATE_CONTINUOUS);
    vscroll->set_update_policy (Gtk::UPDATE_CONTINUOUS);

    vscrollconn = vscroll->signal_value_changed().connect( sigc::mem_fun(*this, &ImageAreaPanel::scrollChanged) );
    hscrollconn = hscroll->signal_value_changed().connect( sigc::mem_fun(*this, &ImageAreaPanel::scrollChanged) );

    imageArea->signal_size_allocate().connect( sigc::mem_fun(*this, &ImageAreaPanel::imageAreaResized) );
}

ImageAreaPanel::~ImageAreaPanel () {

    delete imageArea;
}

void ImageAreaPanel::configScrollBars () {

    int imgw, imgh;
    imageArea->getScrollImageSize (imgw, imgh);

    if (imgw>0 && imgh>0) {
  
        int iw = imageArea->get_width ();
        int ih = imageArea->get_height ();

        hscrollconn.block (true);
        vscrollconn.block (true);

        hscroll->get_adjustment()->set_upper (imgw);
        vscroll->get_adjustment()->set_upper (imgh);
        hscroll->get_adjustment()->set_lower (0);
        vscroll->get_adjustment()->set_lower (0);
        hscroll->get_adjustment()->set_step_increment (imgw/100);
        vscroll->get_adjustment()->set_step_increment (imgh/100);
        hscroll->get_adjustment()->set_page_increment (imgw/5);
        vscroll->get_adjustment()->set_page_increment (imgh/5);
        hscroll->get_adjustment()->set_page_size (iw);
        vscroll->get_adjustment()->set_page_size (ih);
        
        int x, y;
        imageArea->getScrollPosition (x, y);
        hscroll->set_value (x);
        vscroll->set_value (y);

        if (before && this==after)
            before->synchronize ();
        else if (after && this==before)
            after->synchronize ();

        hscrollconn.block (false);
        vscrollconn.block (false);
    }
}

void ImageAreaPanel::refreshScrollBars () {
    
    configScrollBars ();
    queue_draw ();
}

void ImageAreaPanel::imageAreaResized (Gtk::Allocation& req) {

    configScrollBars ();
    queue_draw ();
}

void ImageAreaPanel::scrollChanged () {

    imageArea->setScrollPosition ((int)(hscroll->get_value()), (int)(vscroll->get_value()));
    imageArea->queue_draw ();
#ifdef _WIN32
    gdk_window_process_updates (get_window()->gobj(), true);
#endif
    if (before && this==after) 
        before->synchronize ();
    else if (after && this==before)
        after->synchronize ();
}

void ImageAreaPanel::setBeforeAfterViews (ImageAreaPanel* bef, ImageAreaPanel* aft) {

    before = bef;
    after = aft;
    configScrollBars ();
}

void ImageAreaPanel::zoomChanged () {

    if (after && this==before)
        after->imageArea->setZoom (imageArea->getZoom ());
    else if (before && this==after)
        before->imageArea->setZoom (imageArea->getZoom ());
}

void ImageAreaPanel::synchronize () {

    hscrollconn.block (true);
    vscrollconn.block (true);

    if (after && this==before) {
        int imgw, imgh, x, y;
        after->imageArea->getScrollImageSize (imgw, imgh);
        after->imageArea->getScrollPosition (x, y);
        int bimgw, bimgh;
        imageArea->getScrollImageSize (bimgw, bimgh);
        imageArea->setScrollPosition (x*bimgw/imgw, y*bimgh/imgh);
        imageArea->queue_draw ();
    }
    else if (before && this==after) {
        int imgw, imgh, x, y;
        before->imageArea->getScrollImageSize (imgw, imgh);
        before->imageArea->getScrollPosition (x, y);
        int bimgw, bimgh;
        imageArea->getScrollImageSize (bimgw, bimgh);
        imageArea->setScrollPosition (x*bimgw/imgw, y*bimgh/imgh);
        imageArea->queue_draw ();
    }

    hscrollconn.block (false);
    vscrollconn.block (false);
}

