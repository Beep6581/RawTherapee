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

    iR  = new RTImage ("previewmodeR-on.png");
    iG  = new RTImage ("previewmodeG-on.png");
    iB  = new RTImage ("previewmodeB-on.png");
    iL  = new RTImage ("previewmodeL-on.png");
    iF  = new RTImage ("previewmodeF-on.png");
    
    igR = new RTImage ("previewmodeR-off.png");
    igG = new RTImage ("previewmodeG-off.png");
    igB = new RTImage ("previewmodeB-off.png");
    igL = new RTImage ("previewmodeL-off.png");
    igF = new RTImage ("previewmodeF-off.png");

    previewR = Gtk::manage (new Gtk::ToggleButton ());
    previewR->set_relief(Gtk::RELIEF_NONE);
    previewR->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWR"));
    previewR->set_image(*igR);
    
    previewG = Gtk::manage (new Gtk::ToggleButton ());
    previewG->set_relief(Gtk::RELIEF_NONE);
    previewG->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWG"));
    previewG->set_image(*igG);

    previewB = Gtk::manage (new Gtk::ToggleButton ());
    previewB->set_relief(Gtk::RELIEF_NONE);
    previewB->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWB"));
    previewB->set_image(*igB);

    previewL = Gtk::manage (new Gtk::ToggleButton ());
    previewL->set_relief(Gtk::RELIEF_NONE);
    previewL->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWL"));
    previewL->set_image(*igL);

    previewFocusMask = Gtk::manage (new Gtk::ToggleButton ());
    previewFocusMask->set_relief(Gtk::RELIEF_NONE);
    previewFocusMask->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWFOCUSMASK"));
    previewFocusMask->set_image(*igF);

    previewR->set_active (false);
    previewG->set_active (false);
    previewB->set_active (false);
    previewL->set_active (false);
    previewFocusMask->set_active (false);

    pack_start (*previewR, Gtk::PACK_SHRINK, 0);
    pack_start (*previewG, Gtk::PACK_SHRINK, 0);
    pack_start (*previewB, Gtk::PACK_SHRINK, 0);
    pack_start (*previewL, Gtk::PACK_SHRINK, 0);
    pack_start (*previewFocusMask, Gtk::PACK_SHRINK, 0);

    connR = previewR->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewR) );
    connG = previewG->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewG) );
    connB = previewB->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewB) );
    connL = previewL->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewL) );
    connFocusMask = previewFocusMask->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled),previewFocusMask) );

	//show_all ();
}

PreviewModePanel::~PreviewModePanel (){
	delete iR;
	delete iG;
	delete iB;
	delete iL;
	delete iF;
	delete igR;
	delete igG;
	delete igB;
	delete igL;
	delete igF;
}
//TODO: use functions below for shortcuts
void PreviewModePanel::toggleR () {
    previewR->set_active(!previewR->get_active());
}
void PreviewModePanel::toggleG () {
    previewG->set_active(!previewG->get_active());
}
void PreviewModePanel::toggleB () {
    previewB->set_active(!previewB->get_active());
}
void PreviewModePanel::toggleL () {
    previewL->set_active(!previewL->get_active());
}
void PreviewModePanel::toggleFocusMask () {
    previewFocusMask->set_active(!previewFocusMask->get_active());
}

void PreviewModePanel::buttonToggled (Gtk::ToggleButton* tbpreview) {
	
	connR.block(true);
	connG.block(true);
	connB.block(true);
	connL.block(true);
	connFocusMask.block(true);

	// control state of the buttons
	// only 0 or 1 button at a time can remain pressed
	if (tbpreview!=previewR) previewR->set_active(false);
	if (tbpreview!=previewG) previewG->set_active(false);
	if (tbpreview!=previewB) previewB->set_active(false);
	if (tbpreview!=previewL) previewL->set_active(false);
	if (tbpreview!=previewFocusMask) previewFocusMask->set_active(false);

	// set image based on button's state
	previewR->set_image(previewR->get_active()?*iR:*igR);
	previewG->set_image(previewG->get_active()?*iG:*igG);
	previewB->set_image(previewB->get_active()?*iB:*igB);
	previewL->set_image(previewL->get_active()?*iL:*igL);
	previewFocusMask->set_image(previewFocusMask->get_active()?*iF:*igF);

	connR.block(false);
	connG.block(false);
	connB.block(false);
	connL.block(false);
	connFocusMask.block(false);

	imageArea->queue_draw ();
}
