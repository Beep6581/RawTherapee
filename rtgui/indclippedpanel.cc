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
#include <indclippedpanel.h>
#include <options.h>
#include <multilangmgr.h>
#include <imagearea.h>

IndicateClippedPanel::IndicateClippedPanel (ImageArea* ia) : imageArea(ia) {

    indclippedh = Gtk::manage (new Gtk::ToggleButton ());
    indclippedh->set_relief(Gtk::RELIEF_NONE);
    indclippedh->add (*Gtk::manage (new Gtk::Image (argv0+"/images/warnhl.png")));   
    indclippedh->set_tooltip_text (M("MAIN_TOOLTIP_INDCLIPPEDH"));

    indclippeds = Gtk::manage (new Gtk::ToggleButton ());
    indclippeds->set_relief(Gtk::RELIEF_NONE);
    indclippeds->add (*Gtk::manage (new Gtk::Image (argv0+"/images/warnsh.png")));   
    indclippeds->set_tooltip_text (M("MAIN_TOOLTIP_INDCLIPPEDS"));

    indclippedh->set_active (options.showClippedHighlights);
    indclippeds->set_active (options.showClippedShadows);
	
    pack_start (*indclippedh, Gtk::PACK_SHRINK, 0);
    pack_start (*indclippeds, Gtk::PACK_SHRINK, 0);

    indclippedh->signal_toggled().connect( sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled) );
    indclippeds->signal_toggled().connect( sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled) );
	
	show_all ();
}

void IndicateClippedPanel::buttonToggled () {

	imageArea->queue_draw ();
}
