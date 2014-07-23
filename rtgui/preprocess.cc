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

PreProcess::PreProcess () : FoldableToolPanel(this)
{
	hotDeadPixel = Gtk::manage(new Gtk::CheckButton((M("TP_PREPROCESS_HOTDEADPIXFILT"))));
	hotDeadPixel->set_tooltip_markup (M("TP_PREPROCESS_HOTDEADPIXFILT_TOOLTIP"));

	pack_start( *hotDeadPixel, Gtk::PACK_SHRINK, 4);

	hdpixelconn = hotDeadPixel->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::hotDeadPixelChanged), true);
}

void PreProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();
	hdpixelconn.block (true);

	if(pedited ){
		hotDeadPixel->set_inconsistent (!pedited->raw.hotDeadPixelFilter);
	}

	lastHot = pp->raw.hotdeadpix_filt;

	hotDeadPixel->set_active (pp->raw.hotdeadpix_filt);

	hdpixelconn.block (false);
	enableListener ();
}

void PreProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.hotdeadpix_filt = hotDeadPixel->get_active();

	if (pedited) {
		pedited->raw.hotDeadPixelFilter = !hotDeadPixel->get_inconsistent();
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
