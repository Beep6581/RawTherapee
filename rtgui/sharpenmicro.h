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
 *
 *
 *  Manuel Llorens' algorithm of micro-contrast sharpening
 *
 *
 */
#ifndef _SHARPENMICRO_H_
#define _SHARPENMICRO_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"

class SharpenMicro : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel {

protected:

	Gtk::CheckButton* enabled;
	Gtk::CheckButton* matrix;
	Adjuster* amount;
	Adjuster* uniformity;

	bool lastEnabled;
	sigc::connection enaconn;
	sigc::connection matrixconn;
	bool lastmatrix;

public:

	SharpenMicro           ();

	void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL);
	void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
	void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
	void setBatchMode        (bool batchMode);
	void trimValues          (rtengine::procparams::ProcParams* pp);
	void setAdjusterBehavior (bool amountadd, bool uniformityadd );
	void adjusterChanged     (Adjuster* a, double newval);

	void enabled_toggled     ();
	void matrix_toggled      ();
	

};

#endif
