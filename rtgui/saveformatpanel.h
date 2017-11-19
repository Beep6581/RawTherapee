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
#ifndef __SAVEFORMATPANEL_H__
#define __SAVEFORMATPANEL_H__

#include <gtkmm.h>
#include "adjuster.h"
#include "guiutils.h"
#include "options.h"

class FormatChangeListener
{

public:
    virtual ~FormatChangeListener () {}
    virtual void formatChanged (Glib::ustring f) {}
};

class SaveFormatPanel : public Gtk::Grid, public AdjusterListener
{

protected:
    Adjuster*           jpegQual;
    Gtk::CheckButton*   tiffUncompressed;
    MyComboBoxText*     format;
    MyComboBoxText*     jpegSubSamp;
    Gtk::Grid*          formatOpts;
    Gtk::Grid*          jpegOpts;
    Gtk::Label*         jpegSubSampLabel;
    FormatChangeListener* listener;
    Glib::ustring       fstr[6];
    Gtk::CheckButton*   savesPP;


public:

    SaveFormatPanel ();
    ~SaveFormatPanel ();
    void        setListener     (FormatChangeListener* l)
    {
        listener = l;
    }

    void        init            (SaveFormat& sf);
    SaveFormat  getFormat       ();

    void        formatChanged   ();
    void        adjusterChanged (Adjuster* a, double newval);
};

#endif
