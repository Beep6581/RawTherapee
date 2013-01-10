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
#include "zoompanel.h"
#include "multilangmgr.h"
#include "imagearea.h"
#include "rtimage.h"

ZoomPanel::ZoomPanel (ImageArea* iarea) : iarea(iarea) {

    set_border_width (0);

    Gtk::Image* imageOut = Gtk::manage (new RTImage ("gtk-zoom-out.png"));
    imageOut->set_padding(0,0);
    Gtk::Image* imageIn = Gtk::manage (new RTImage ("gtk-zoom-in.png"));
    imageIn->set_padding(0,0);
    Gtk::Image* image11 =Gtk::manage ( new RTImage ("gtk-zoom-100.png"));
    image11->set_padding(0,0);
    Gtk::Image* imageFit = Gtk::manage (new RTImage ("gtk-zoom-fit.png"));
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

    Gtk::Image* imageCrop = Gtk::manage (new RTImage ("new-detail-window.png"));
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

    zoomIn->set_tooltip_markup (M("ZOOMPANEL_ZOOMIN"));
    zoomOut->set_tooltip_markup (M("ZOOMPANEL_ZOOMOUT"));
    zoom11->set_tooltip_markup (M("ZOOMPANEL_ZOOM100"));
    zoomFit->set_tooltip_markup (M("ZOOMPANEL_ZOOMFITSCREEN"));
    newCrop->set_tooltip_markup (M("ZOOMPANEL_NEWCROPWINDOW"));

    zoomLabel->set_text (M("ZOOMPANEL_100"));
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
		if (z<100) {
			zoomLabel->set_text (Glib::ustring::compose(" %1%%", z));
		} else {
			zoomLabel->set_text (Glib::ustring::compose("%1%%", z));
		}
    }
}

void ZoomPanel::newCropClicked () {
    
    iarea->addCropWindow ();
}
