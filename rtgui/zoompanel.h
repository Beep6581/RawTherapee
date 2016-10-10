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
 *l
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _ZOOMPANEL_
#define _ZOOMPANEL_

#include <gtkmm.h>

class ImageArea;
class ZoomPanel : public Gtk::Grid
{

protected:

    Gtk::Button*    zoomOut;
    Gtk::Button*    zoomIn;
    Gtk::Button*    zoomFit;
    Gtk::Button*    zoomFitCrop;
    Gtk::Button*    zoom11;
    Gtk::Button*    newCrop;
    Gtk::Label*     zoomLabel;
    ImageArea*      iarea;

public:

    explicit ZoomPanel (ImageArea* iarea);

    void zoomInClicked      ();
    void zoomOutClicked     ();
    void zoomFitClicked     ();
    void zoomFitCropClicked ();
    void zoom11Clicked      ();
    void newCropClicked     ();
    void refreshZoomLabel   ();
};

#endif

