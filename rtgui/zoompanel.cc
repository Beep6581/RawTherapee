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
#include <zoompanel.h>
#include <multilangmgr.h>
#include <imagearea.h>

ZoomPanel::ZoomPanel (ImageArea* iarea) : iarea(iarea) {

    set_border_width (0);

    Gtk::Label* label = Gtk::manage (new Gtk::Label (Glib::ustring("<b>") + "Zoom" + ":</b>  "));
    label->set_use_markup (true);
    
    pack_start (*label, Gtk::PACK_SHRINK, 4);
    
    Gtk::Image* imageOut = Gtk::manage (new Gtk::Image (Gtk::StockID ("gtk-zoom-out"), Gtk::ICON_SIZE_SMALL_TOOLBAR));
    imageOut->set_padding(0,0);
    Gtk::Image* imageIn = Gtk::manage (new Gtk::Image (Gtk::StockID ("gtk-zoom-in"), Gtk::ICON_SIZE_SMALL_TOOLBAR));
    imageIn->set_padding(0,0);
    Gtk::Image* image11 =Gtk::manage ( new Gtk::Image (Gtk::StockID ("gtk-zoom-100"), Gtk::ICON_SIZE_SMALL_TOOLBAR));
    image11->set_padding(0,0);
    Gtk::Image* imageFit = Gtk::manage (new Gtk::Image (Gtk::StockID ("gtk-zoom-fit"), Gtk::ICON_SIZE_SMALL_TOOLBAR));
    imageFit->set_padding(0,0);

    zoomOut = Gtk::manage (new Gtk::Button());
    zoomOut->add (*imageOut);
    zoomOut->set_relief(Gtk::RELIEF_NONE);
    zoomIn = Gtk::manage (new Gtk::Button());
    zoomIn->add (*imageIn);
    zoomIn->set_relief(Gtk::RELIEF_NONE);
    zoomFit = Gtk::manage (new Gtk::Button());
    zoomFit->add (*imageFit);
    zoomFit->set_relief(Gtk::RELIEF_NONE);
    zoom11 = Gtk::manage (new Gtk::Button());
    zoom11->add (*image11);
    zoom11->set_relief(Gtk::RELIEF_NONE);

    pack_start (*zoomOut, Gtk::PACK_SHRINK, 0);
    pack_start (*zoomIn, Gtk::PACK_SHRINK, 0);
    pack_start (*zoomFit, Gtk::PACK_SHRINK, 0);
    pack_start (*zoom11, Gtk::PACK_SHRINK, 0);

    zoomLabel = Gtk::manage (new Gtk::Label ());
    pack_start (*zoomLabel, Gtk::PACK_SHRINK, 4);

    Gtk::Image* imageCrop = Gtk::manage (new Gtk::Image (Gtk::StockID ("gtk-add"), Gtk::ICON_SIZE_SMALL_TOOLBAR));
    imageCrop->set_padding(0,0);
    newCrop = Gtk::manage (new Gtk::Button());
    newCrop->add (*imageCrop);
    newCrop->set_relief(Gtk::RELIEF_NONE);
    pack_start (*newCrop, Gtk::PACK_SHRINK, 4);

    show_all_children ();

    zoomIn->signal_clicked().connect ( sigc::mem_fun(*this, &ZoomPanel::zoomInClicked) );
    zoomOut->signal_clicked().connect( sigc::mem_fun(*this, &ZoomPanel::zoomOutClicked) );
    zoomFit->signal_clicked().connect( sigc::mem_fun(*this, &ZoomPanel::zoomFitClicked) );
    zoom11->signal_clicked().connect ( sigc::mem_fun(*this, &ZoomPanel::zoom11Clicked) );
    newCrop->signal_clicked().connect ( sigc::mem_fun(*this, &ZoomPanel::newCropClicked) );

    zoomIn->set_tooltip_text ("Zoom In");
    zoomOut->set_tooltip_text ("Zoom Out");
    zoom11->set_tooltip_text ("Zoom to 100%");
    zoomFit->set_tooltip_text ("Fit to screen");
    newCrop->set_tooltip_text ("Add new crop window");

    zoomLabel->set_text ("(100%)");
}

void ZoomPanel::zoomInClicked () {

    if (iarea->mainCropWindow)
        iarea->mainCropWindow->zoomIn ();
}

void ZoomPanel::zoomOutClicked () {

    if (iarea->mainCropWindow) 
        iarea->mainCropWindow->zoomOut ();
}

void ZoomPanel::zoomFitClicked () {

    if (iarea->mainCropWindow) 
        iarea->mainCropWindow->zoomFit ();
}

void ZoomPanel::zoom11Clicked () {

    if (iarea->mainCropWindow) 
        iarea->mainCropWindow->zoom11 ();
}

void ZoomPanel::refreshZoomLabel () {

    if (iarea->mainCropWindow) {
        int z = (int)(iarea->mainCropWindow->getZoom () * 100);
        zoomLabel->set_text (Glib::ustring::compose("%1%%", z));
    }
}

void ZoomPanel::newCropClicked () {
    
    iarea->addCropWindow ();
}
