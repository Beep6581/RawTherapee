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
#ifndef _FLATFIELD_H_
#define _FLATFIELD_H_

#include <auto_ptr.h>
#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "../rtengine/rawimage.h"
#include "guiutils.h"

class FFProvider {
  public:
    virtual rtengine::RawImage* getFF() = 0;
    // add other info here
};

class FlatField : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel {

protected:

	MyFileChooserButton *flatFieldFile;
	std::auto_ptr<FileChooserLastFolderPersister> flatFieldFilePersister;
	Gtk::Label *ffLabel;
	Gtk::Label *ffInfo;
	Gtk::Button *flatFieldFileReset;
	Gtk::CheckButton* flatFieldAutoSelect;
	Adjuster* flatFieldBlurRadius;
	MyComboBoxText* flatFieldBlurType;
	Gtk::HBox *hbff;
	bool ffChanged;
	bool lastFFAutoSelect;
	FFProvider *ffp;
	sigc::connection flatFieldFileconn, flatFieldAutoSelectconn, flatFieldBlurTypeconn;
        Glib::ustring lastShortcutPath;

public:

	FlatField ();

	void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);

    void adjusterChanged            (Adjuster* a, double newval);
    void flatFieldFileChanged       ();
    void flatFieldFile_Reset        ();
    void flatFieldAutoSelectChanged ();
    void flatFieldBlurTypeChanged   ();
    void setShortcutPath(Glib::ustring path);
    void setFFProvider              (FFProvider* p) { ffp = p; };
};

#endif
