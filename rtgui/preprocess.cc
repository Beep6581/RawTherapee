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
#include "preprocess.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>

using namespace rtengine;
using namespace rtengine::procparams;

PreProcess::PreProcess () : Gtk::VBox(), FoldableToolPanel(this)
{
	set_border_width(4);

	hotDeadPixel = Gtk::manage(new Gtk::CheckButton((M("TP_PREPROCESS_HOTDEADPIXFILT"))));
	hotDeadPixel->set_tooltip_markup (M("TP_PREPROCESS_HOTDEADPIXFILT_TOOLTIP"));

	lineDenoise = Gtk::manage(new Adjuster (M("TP_PREPROCESS_LINEDENOISE"),0,1000,1,0));
	lineDenoise->setAdjusterListener (this);
	if (lineDenoise->delay < 1000) lineDenoise->delay = 1000;
	lineDenoise->show();

	greenEqThreshold = Gtk::manage(new Adjuster (M("TP_PREPROCESS_GREENEQUIL"),0,100,1,0));
	greenEqThreshold->setAdjusterListener (this);
	if (greenEqThreshold->delay < 1000) greenEqThreshold->delay = 1000;
	greenEqThreshold->show();

	pack_start( *lineDenoise, Gtk::PACK_SHRINK, 4);

	pack_start( *Gtk::manage (new  Gtk::HSeparator()));

	pack_start( *greenEqThreshold, Gtk::PACK_SHRINK, 4);

	pack_start( *Gtk::manage (new  Gtk::HSeparator()));

	pack_start( *hotDeadPixel, Gtk::PACK_SHRINK, 4);

	hdpixelconn = hotDeadPixel->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::hotDeadPixelChanged), true);
}

void PreProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();
	hdpixelconn.block (true);

	if(pedited ){
		hotDeadPixel->set_inconsistent (!pedited->raw.hotDeadPixelFilter);
		lineDenoise->setEditedState( pedited->raw.linenoise ? Edited : UnEdited );
		greenEqThreshold->setEditedState( pedited->raw.greenEq ? Edited : UnEdited );
	}

	lastHot = pp->raw.hotdeadpix_filt;

	hotDeadPixel->set_active (pp->raw.hotdeadpix_filt);
	lineDenoise->setValue (pp->raw.linenoise);
	greenEqThreshold->setValue (pp->raw.greenthresh);

	hdpixelconn.block (false);
	enableListener ();
}

void PreProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.hotdeadpix_filt = hotDeadPixel->get_active();
	pp->raw.linenoise = lineDenoise->getIntValue();
	pp->raw.greenthresh = greenEqThreshold->getIntValue();

	if (pedited) {
		pedited->raw.linenoise = lineDenoise->getEditedState ();
		pedited->raw.greenEq= greenEqThreshold->getEditedState ();
		pedited->raw.hotDeadPixelFilter = !hotDeadPixel->get_inconsistent();
	}
}

void PreProcess::adjusterChanged (Adjuster* a, double newval)
{
	if (listener) {

		Glib::ustring value = a->getTextValue();

		if (a == greenEqThreshold)
			listener->panelChanged (EvPreProcessGEquilThresh,  value );
		else if (a == lineDenoise)
			listener->panelChanged (EvPreProcessLineDenoise,  value );
	}
}

void PreProcess::setBatchMode(bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	lineDenoise->showEditedCB ();
	greenEqThreshold->showEditedCB ();
}

void PreProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	lineDenoise->setDefault( defParams->raw.linenoise);
	greenEqThreshold->setDefault (defParams->raw.greenthresh);

	if (pedited) {
		lineDenoise->setDefaultEditedState( pedited->raw.linenoise ? Edited : UnEdited);
		greenEqThreshold->setDefaultEditedState(pedited->raw.greenEq ? Edited : UnEdited);
	} else {
		lineDenoise->setDefaultEditedState( Irrelevant );
		greenEqThreshold->setDefaultEditedState(Irrelevant );
	}
}

void PreProcess::hotDeadPixelChanged ()
{
    if (batchMode) {
        if (hotDeadPixel->get_inconsistent()) {
        	hotDeadPixel->set_inconsistent (false);
        	hdpixelconn.block (true);
        	hotDeadPixel->set_active (false);
        	hdpixelconn.block (false);
        }
        else if (lastHot)
        	hotDeadPixel->set_inconsistent (true);

        lastHot = hotDeadPixel->get_active ();
    }
    if (listener)
        listener->panelChanged (EvPreProcessHotDeadPixel, hotDeadPixel->get_active()?M("GENERAL_ENABLED"):M("GENERAL_DISABLED"));
}

void PreProcess::setAdjusterBehavior (bool linedenoiseadd, bool greenequiladd) {

   	lineDenoise->setAddMode(linedenoiseadd);
   	greenEqThreshold->setAddMode(greenequiladd);
}

void PreProcess::trimValues (rtengine::procparams::ProcParams* pp) {

	lineDenoise->trimValue(pp->raw.linenoise);
	greenEqThreshold->trimValue(pp->raw.greenthresh);
}
