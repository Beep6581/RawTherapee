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
#ifndef _VIBRANCE_
#define _VIBRANCE_

#include <gtkmm.h>
#include "adjuster.h"
#include "thresholdadjuster.h"
#include "toolpanel.h"

class Vibrance : public Gtk::VBox, public AdjusterListener, public ThresholdAdjusterListener, public FoldableToolPanel, public ThresholdCurveProvider {

protected:
	Gtk::CheckButton* enabled;
	Adjuster* pastels;
	Adjuster* saturated;
    ThresholdAdjuster* psThreshold;
	Gtk::CheckButton* protectSkins;
	Gtk::CheckButton* avoidColorShift;
	Gtk::CheckButton* pastSatTog;
    bool lastEnabled;
	bool lastProtectSkins;
	bool lastAvoidColorShift;
	bool lastPastSatTog;

	sigc::connection enaconn;
	sigc::connection pskinsconn;
	sigc::connection ashiftconn;
	sigc::connection pastsattogconn;

public:

	Vibrance                 ();

	void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL);
	void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
	void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
	void setBatchMode        (bool batchMode);
	void trimValues          (rtengine::procparams::ProcParams* pp);
	void setAdjusterBehavior (bool amountadd, bool passadd, bool psthreshdadd);
	void adjusterChanged     (Adjuster* a, double newval);
	void adjusterChanged     (ThresholdAdjuster* a, int newBottom, int newTop);

	void enabled_toggled         ();
	void protectskins_toggled    ();
	void avoidcolorshift_toggled ();
	void pastsattog_toggled      ();
	std::vector<double> getCurvePoints(ThresholdSelector* tAdjuster) const;
};


#endif
