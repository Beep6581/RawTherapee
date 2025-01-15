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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "zoompanel.h"
#include "multilangmgr.h"
#include "imagearea.h"
#include "rtimage.h"

ZoomPanel::ZoomPanel (ImageArea* iarea) : iarea(iarea)
{
    set_name ("EditorZoomPanel");

    Gtk::Image* imageOut = Gtk::manage (new RTImage ("magnifier-minus", Gtk::ICON_SIZE_LARGE_TOOLBAR));
    imageOut->set_padding(0, 0);
    Gtk::Image* imageIn = Gtk::manage (new RTImage ("magnifier-plus", Gtk::ICON_SIZE_LARGE_TOOLBAR));
    imageIn->set_padding(0, 0);
    Gtk::Image* image11 = Gtk::manage ( new RTImage ("magnifier-1to1", Gtk::ICON_SIZE_LARGE_TOOLBAR));
    image11->set_padding(0, 0);
    Gtk::Image* imageFit = Gtk::manage (new RTImage ("magnifier-fit", Gtk::ICON_SIZE_LARGE_TOOLBAR));
    imageFit->set_padding(0, 0);
    Gtk::Image* imageFitCrop = Gtk::manage (new RTImage ("magnifier-crop", Gtk::ICON_SIZE_LARGE_TOOLBAR));
    imageFit->set_padding(0, 0);

    zoomOut = Gtk::manage (new Gtk::Button());
    zoomOut->add (*imageOut);
    zoomOut->set_relief(Gtk::RELIEF_NONE);
    setExpandAlignProperties(zoomOut, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    zoomIn = Gtk::manage (new Gtk::Button());
    zoomIn->add (*imageIn);
    zoomIn->set_relief(Gtk::RELIEF_NONE);
    setExpandAlignProperties(zoomIn, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    zoomFit = Gtk::manage (new Gtk::Button());
    zoomFit->add (*imageFit);
    zoomFit->set_relief(Gtk::RELIEF_NONE);
    setExpandAlignProperties(zoomFit, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    zoomFitCrop = Gtk::manage (new Gtk::Button());
    zoomFitCrop->add (*imageFitCrop);
    zoomFitCrop->set_relief(Gtk::RELIEF_NONE);
    setExpandAlignProperties(zoomFitCrop, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    zoom11 = Gtk::manage (new Gtk::Button());
    zoom11->add (*image11);
    zoom11->set_relief(Gtk::RELIEF_NONE);
    setExpandAlignProperties(zoom11, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);

    attach_next_to (*zoomOut, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*zoomIn, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*zoomFit, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*zoomFitCrop, Gtk::POS_RIGHT, 1, 1);
    attach_next_to (*zoom11, Gtk::POS_RIGHT, 1, 1);

    zoomLabel = Gtk::manage (new Gtk::Label ());
    setExpandAlignProperties(zoomLabel, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    attach_next_to (*zoomLabel, Gtk::POS_RIGHT, 1, 1);

    Gtk::Image* imageCrop = Gtk::manage (new RTImage ("window-add", Gtk::ICON_SIZE_LARGE_TOOLBAR));
    imageCrop->set_padding(0, 0);
    newCrop = Gtk::manage (new Gtk::Button());
    newCrop->add (*imageCrop);
    newCrop->set_relief(Gtk::RELIEF_NONE);
    setExpandAlignProperties(newCrop, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    attach_next_to (*newCrop, Gtk::POS_RIGHT, 1, 1);

    show_all_children ();

    zoomIn->signal_clicked().connect ( sigc::mem_fun(*this, &ZoomPanel::zoomInClicked) );
    zoomOut->signal_clicked().connect( sigc::mem_fun(*this, &ZoomPanel::zoomOutClicked) );
    zoomFit->signal_clicked().connect( sigc::mem_fun(*this, &ZoomPanel::zoomFitClicked) );
    zoomFitCrop->signal_clicked().connect( sigc::mem_fun(*this, &ZoomPanel::zoomFitCropClicked) );
    zoom11->signal_clicked().connect ( sigc::mem_fun(*this, &ZoomPanel::zoom11Clicked) );
    newCrop->signal_clicked().connect ( sigc::mem_fun(*this, &ZoomPanel::newCropClicked) );

    zoomIn->set_tooltip_markup (M("ZOOMPANEL_ZOOMIN"));
    zoomOut->set_tooltip_markup (M("ZOOMPANEL_ZOOMOUT"));
    zoom11->set_tooltip_markup (M("ZOOMPANEL_ZOOM100"));
    zoomFit->set_tooltip_markup (M("ZOOMPANEL_ZOOMFITSCREEN"));
    zoomFitCrop->set_tooltip_markup (M("ZOOMPANEL_ZOOMFITCROPSCREEN"));
    newCrop->set_tooltip_markup (M("ZOOMPANEL_NEWCROPWINDOW"));

    zoomLabel->set_text (M("ZOOMPANEL_100"));
}

void ZoomPanel::zoomInClicked ()
{

    if (iarea->mainCropWindow) {
        iarea->mainCropWindow->zoomIn ();
    }
}

void ZoomPanel::zoomOutClicked ()
{

    if (iarea->mainCropWindow) {
        iarea->mainCropWindow->zoomOut ();
    }
}

void ZoomPanel::zoomFitClicked ()
{

    if (iarea->mainCropWindow) {
        iarea->mainCropWindow->zoomFit ();
    }
}

void ZoomPanel::zoomFitCropClicked ()
{

    if (iarea->mainCropWindow) {
        iarea->mainCropWindow->zoomFitCrop ();
    }
}

void ZoomPanel::zoom11Clicked ()
{

    if (iarea->mainCropWindow) {
        iarea->mainCropWindow->zoom11 ();
    }
}

void ZoomPanel::refreshZoomLabel ()
{

    if (iarea->mainCropWindow) {
        int z = (int)(iarea->mainCropWindow->getZoom () * 100);

        if (z < 100) {
            zoomLabel->set_text (Glib::ustring::compose(" %1%%", z));
        } else {
            zoomLabel->set_text (Glib::ustring::compose("%1%%", z));
        }
    }
}

void ZoomPanel::newCropClicked ()
{

    iarea->addCropWindow ();
}
