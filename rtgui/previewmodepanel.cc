/*
 *  This file is part of RawTherapee.
 *
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
#include "previewmodepanel.h"
#include "options.h"
#include "multilangmgr.h"
#include "imagearea.h"
#include "rtimage.h"

PreviewModePanel::PreviewModePanel (ImageArea* ia) : imageArea(ia) {

    previewR = Gtk::manage (new Gtk::ToggleButton ("R"));
    previewR->set_relief(Gtk::RELIEF_NONE);
    previewR->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWR"));
    
    previewG = Gtk::manage (new Gtk::ToggleButton ("G"));
    previewG->set_relief(Gtk::RELIEF_NONE);
    previewG->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWG"));

    previewB = Gtk::manage (new Gtk::ToggleButton ("B"));
    previewB->set_relief(Gtk::RELIEF_NONE);
    previewB->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWB"));

    previewL = Gtk::manage (new Gtk::ToggleButton ("L"));
    previewL->set_relief(Gtk::RELIEF_NONE);
    previewL->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWL"));

    previewFocusMask = Gtk::manage (new Gtk::ToggleButton ("F"));
    previewFocusMask->set_relief(Gtk::RELIEF_NONE);
    previewFocusMask->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWFOCUSMASK"));
    previewFocusMask->hide();//TODO re-enable when Focus Mask is developed

    previewR->set_active (false);
    previewG->set_active (false);
    previewB->set_active (false);
    previewL->set_active (false);
    previewFocusMask->set_active (false);

    pack_start (*previewR, Gtk::PACK_SHRINK, 0);
    pack_start (*previewG, Gtk::PACK_SHRINK, 0);
    pack_start (*previewB, Gtk::PACK_SHRINK, 0);
    pack_start (*previewL, Gtk::PACK_SHRINK, 0);
    //pack_start (*previewFocusMask, Gtk::PACK_SHRINK, 0); //TODO re-enable when Focus Mask is developed

    connR = previewR->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewR) );
    connG = previewG->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewG) );
    connB = previewB->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewB) );
    connL = previewL->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewL) );
    connFocusMask = previewFocusMask->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewFocusMask) );

	//show_all ();
}

//TODO: use functions below for shortcuts
void PreviewModePanel::toggleR () {

}
void PreviewModePanel::toggleG () {

}
void PreviewModePanel::toggleB () {

}
void PreviewModePanel::toggleL () {

}
void PreviewModePanel::toggleFocusMask () {

}

void PreviewModePanel::buttonToggled (Gtk::ToggleButton* tbpreview) {
	// only 0 or 1 button at a time can remain pressed
	
	connR.block(true);
	connG.block(true);
	connB.block(true);
	connL.block(true);
	connFocusMask.block(true);

	if (tbpreview==previewR){
		//
		previewG->set_active(false);
		previewB->set_active(false);
		previewL->set_active(false);
		previewFocusMask->set_active(false);
	}

	if (tbpreview==previewG){
		previewR->set_active(false);
		//
		previewB->set_active(false);
		previewL->set_active(false);
		previewFocusMask->set_active(false);
	}

	if (tbpreview==previewB){
		previewR->set_active(false);
		previewG->set_active(false);
		//
		previewL->set_active(false);
		previewFocusMask->set_active(false);
	}

	if (tbpreview==previewL){
		previewR->set_active(false);
		previewG->set_active(false);
		previewB->set_active(false);
		//
		previewFocusMask->set_active(false);
	}

	if (tbpreview==previewFocusMask){
		previewR->set_active(false);
		previewG->set_active(false);
		previewB->set_active(false);
		previewL->set_active(false);
		//
	}

	connR.block(false);
	connG.block(false);
	connB.block(false);
	connL.block(false);
	connFocusMask.block(false);

	imageArea->queue_draw ();
}
