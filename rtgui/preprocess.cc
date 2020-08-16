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
#include <sstream>

#include "preprocess.h"

#include "guiutils.h"
#include "options.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

PreProcess::PreProcess () : FoldableToolPanel(this, "preprocess", M("TP_PREPROCESS_LABEL"), options.prevdemo != PD_Sidecar)
{

    Gtk::HBox* hotdeadPixel = Gtk::manage( new Gtk::HBox () );
    hotdeadPixel->set_spacing(4);
    hotPixel = Gtk::manage(new Gtk::CheckButton((M("TP_PREPROCESS_HOTPIXFILT"))));
    deadPixel = Gtk::manage(new Gtk::CheckButton((M("TP_PREPROCESS_DEADPIXFILT"))));

    hotPixel->set_tooltip_markup (M("TP_PREPROCESS_HOTPIXFILT_TOOLTIP"));
    deadPixel->set_tooltip_markup (M("TP_PREPROCESS_DEADPIXFILT_TOOLTIP"));

    hotdeadPixel->pack_start( *hotPixel, Gtk::PACK_SHRINK);
    hotdeadPixel->pack_start( *deadPixel, Gtk::PACK_SHRINK, 0);
    pack_start(*hotdeadPixel, Gtk::PACK_SHRINK, 0);
    hdThreshold = Gtk::manage (new Adjuster (M("TP_RAW_HD"), 20, 200, 2, 100));
    hdThreshold->set_tooltip_markup (M("TP_RAW_HD_TOOLTIP"));
    hdThreshold->setAdjusterListener (this);

    hdThreshold->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    hdThreshold->show();
    pack_start( *hdThreshold, Gtk::PACK_SHRINK, 4);

//  hotdeadPixel->show();
    hpixelconn = hotPixel->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::hotPixelChanged), true);
    dpixelconn = deadPixel->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::deadPixelChanged), true);
}

void PreProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    hpixelconn.block (true);
    dpixelconn.block (true);

    if(pedited ) {
        hotPixel->set_inconsistent (!pedited->raw.hotPixelFilter);
        deadPixel->set_inconsistent (!pedited->raw.deadPixelFilter);
    }

    lastHot = pp->raw.hotPixelFilter;
    lastDead = pp->raw.deadPixelFilter;
    hotPixel->set_active (pp->raw.hotPixelFilter);
    deadPixel->set_active (pp->raw.deadPixelFilter);
    hdThreshold->setValue (pp->raw.hotdeadpix_thresh);
    hpixelconn.block (false);
    dpixelconn.block (false);
    enableListener ();
}

void PreProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.hotPixelFilter = hotPixel->get_active();
    pp->raw.deadPixelFilter = deadPixel->get_active();
    pp->raw.hotdeadpix_thresh = hdThreshold->getIntValue();

    if (pedited) {
        pedited->raw.hotdeadpix_thresh = hdThreshold->getEditedState ();
        pedited->raw.hotPixelFilter = !hotPixel->get_inconsistent();
        pedited->raw.deadPixelFilter = !deadPixel->get_inconsistent();
    }
}

void PreProcess::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        if (a == hdThreshold) {
            listener->panelChanged (EvPreProcessHotDeadThresh, a->getTextValue() );
        }
    }
}

void PreProcess::hotPixelChanged ()
{
    if (batchMode) {
        if (hotPixel->get_inconsistent()) {
            hotPixel->set_inconsistent (false);
            hpixelconn.block (true);
            hotPixel->set_active (false);
            hpixelconn.block (false);
        } else if (lastHot) {
            hotPixel->set_inconsistent (true);
        }

        lastHot = hotPixel->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPreProcessHotPixel, hotPixel->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void PreProcess::deadPixelChanged ()
{
    if (batchMode) {
        if (deadPixel->get_inconsistent()) {
            deadPixel->set_inconsistent (false);
            dpixelconn.block (true);
            deadPixel->set_active (false);
            dpixelconn.block (false);
        } else if (lastDead) {
            deadPixel->set_inconsistent (true);
        }

        lastDead = deadPixel->get_active ();
    }

    if (listener) {
        listener->panelChanged (EvPreProcessDeadPixel, deadPixel->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
