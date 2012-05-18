/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2011 Oliver Duis <oduis@oliverduis.de>
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
#ifndef _LENSPROFILE_H_
#define _LENSPROFILE_H_

#include <gtkmm.h>
#include "toolpanel.h"
#include "guiutils.h"

class LensProfilePanel : public Gtk::VBox, public FoldableToolPanel {

protected:

    MyFileChooserButton *fcbLCPFile;
    Gtk::HBox *hbLCPFile;
    Gtk::Button *btnReset;
    Gtk::Label *lLCPFileHead;
    bool lcpFileChanged;
    sigc::connection conLCPFile;

public:

    LensProfilePanel ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);

    void onLCPFileChanged ();
    void onLCPFileReset   ();
};

#endif
