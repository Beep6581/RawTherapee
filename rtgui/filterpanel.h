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
#ifndef _FILTERPANEL_
#define _FILTERPANEL_

#include <gtkmm.h>
#include "exiffiltersettings.h"

class FilterPanelListener
{

public:
    virtual void exifFilterChanged () {}
};

class FilterPanel : public Gtk::VBox
{

protected:
    Gtk::ListViewText*      filetype;
    Gtk::ListViewText*      camera;
    Gtk::ListViewText*      lens;
    Gtk::ListViewText*      expcomp;
    Gtk::Entry* fnumberFrom;
    Gtk::Entry* fnumberTo;
    Gtk::Entry* shutterFrom;
    Gtk::Entry* shutterTo;
    Gtk::Entry* focalFrom;
    Gtk::Entry* focalTo;
    Gtk::Entry* isoFrom;
    Gtk::Entry* isoTo;
    Gtk::CheckButton* enabled;
    Gtk::CheckButton* enaFNumber;
    Gtk::CheckButton* enaShutter;
    Gtk::CheckButton* enaFocalLen;
    Gtk::CheckButton* enaISO;
    Gtk::CheckButton* enaExpComp;
    Gtk::CheckButton* enaCamera;
    Gtk::CheckButton* enaLens;
    Gtk::CheckButton* enaFiletype;

    int conns;
    sigc::connection sChange[22];

    ExifFilterSettings curefs;
    FilterPanelListener* listener;

public:
    FilterPanel ();

    void setFilterPanelListener (FilterPanelListener* l)
    {
        listener = l;
    }

    void setFilter (ExifFilterSettings& defefs, bool updateLists);
    ExifFilterSettings getFilter ();
    bool isEnabled ();

    void valueChanged ();
    void setEnabled(bool enabledState)
    {
        enabled->set_active(enabledState);
    }
};

#endif
