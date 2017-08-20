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

PreviewModePanel::PreviewModePanel (ImageArea* ia) : imageArea(ia)
{

    iR  = new RTImage ("previewmodeR-on.png");
    iG  = new RTImage ("previewmodeG-on.png");
    iB  = new RTImage ("previewmodeB-on.png");
    iL  = new RTImage ("previewmodeL-on.png");
    iF  = new RTImage ("previewmodeF-on.png");
    iBC0 = new RTImage ("previewmodeBC0-on.png");
    iBC1 = new RTImage ("previewmodeBC1-on.png");
    iBC2 = new RTImage ("previewmodeBC2-on.png");
    iBC3 = new RTImage ("previewmodeBC3-on.png");

    igR = new RTImage ("previewmodeR-off.png");
    igG = new RTImage ("previewmodeG-off.png");
    igB = new RTImage ("previewmodeB-off.png");
    igL = new RTImage ("previewmodeL-off.png");
    igF = new RTImage ("previewmodeF-off.png");
    igBC0 = new RTImage ("previewmodeBC0-off.png");
    igBC1 = new RTImage ("previewmodeBC1-off.png");
    igBC2 = new RTImage ("previewmodeBC2-off.png");
    igBC3 = new RTImage ("previewmodeBC3-off.png");

    backColor0 = Gtk::manage (new Gtk::ToggleButton ());
    backColor0->set_relief(Gtk::RELIEF_NONE);
    backColor0->set_tooltip_markup (M("MAIN_TOOLTIP_BACKCOLOR0"));
    backColor0->set_image(options.bgcolor == 0 ? *iBC0 : *igBC0);

    backColor1 = Gtk::manage (new Gtk::ToggleButton ());
    backColor1->set_relief(Gtk::RELIEF_NONE);
    backColor1->set_tooltip_markup (M("MAIN_TOOLTIP_BACKCOLOR1"));
    backColor1->set_image(options.bgcolor == 1 ? *iBC1 : *igBC1);

    backColor3 = Gtk::manage (new Gtk::ToggleButton ());
    backColor3->set_relief(Gtk::RELIEF_NONE);
    backColor3->set_tooltip_markup (M("MAIN_TOOLTIP_BACKCOLOR3"));
    backColor3->set_image(options.bgcolor == 2 ? *iBC3 : *igBC3);

    backColor2 = Gtk::manage (new Gtk::ToggleButton ());
    backColor2->set_relief(Gtk::RELIEF_NONE);
    backColor2->set_tooltip_markup (M("MAIN_TOOLTIP_BACKCOLOR2"));
    backColor2->set_image(options.bgcolor == 2 ? *iBC2 : *igBC2);

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

    backColor0->set_active (options.bgcolor == 0);
    backColor1->set_active (options.bgcolor == 1);
    backColor2->set_active (options.bgcolor == 2);
    backColor3->set_active (options.bgcolor == 3);

    vbbackColor = Gtk::manage (new Gtk::VBox ());
    vbbackColor->pack_start (*backColor0, Gtk::PACK_SHRINK, 0);
    vbbackColor->pack_start (*backColor1, Gtk::PACK_SHRINK, 0);
    vbbackColor->pack_start (*backColor2, Gtk::PACK_SHRINK, 0);
    vbbackColor->pack_start (*backColor3, Gtk::PACK_SHRINK, 0);
    pack_start (*vbbackColor, Gtk::PACK_SHRINK, 0);

    pack_start (*previewR, Gtk::PACK_SHRINK, 0);
    pack_start (*previewG, Gtk::PACK_SHRINK, 0);
    pack_start (*previewB, Gtk::PACK_SHRINK, 0);
    pack_start (*previewL, Gtk::PACK_SHRINK, 0);
    pack_start (*previewFocusMask, Gtk::PACK_SHRINK, 0);

    connR = previewR->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewR) );
    connG = previewG->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewG) );
    connB = previewB->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewB) );
    connL = previewL->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewL) );
    connFocusMask = previewFocusMask->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewFocusMask) );

    connbackColor0 = backColor0->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled_backColor), backColor0) );
    connbackColor1 = backColor1->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled_backColor), backColor1) );
    connbackColor2 = backColor2->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled_backColor), backColor2) );
    connbackColor3 = backColor3->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled_backColor), backColor3) );

    //show_all ();
}

PreviewModePanel::~PreviewModePanel ()
{
    delete iR;
    delete iG;
    delete iB;
    delete iL;
    delete iF;
    delete iBC0;
    delete iBC1;
    delete iBC2;
    delete iBC3;
    delete igR;
    delete igG;
    delete igB;
    delete igL;
    delete igF;
    delete igBC0;
    delete igBC1;
    delete igBC2;
    delete igBC3;
}
//toggle Functions below are for shortcuts
void PreviewModePanel::toggleR ()
{
    previewR->set_active(!previewR->get_active());
}
void PreviewModePanel::toggleG ()
{
    previewG->set_active(!previewG->get_active());
}
void PreviewModePanel::toggleB ()
{
    previewB->set_active(!previewB->get_active());
}
void PreviewModePanel::toggleL ()
{
    previewL->set_active(!previewL->get_active());
}
void PreviewModePanel::toggleFocusMask ()
{
    previewFocusMask->set_active(!previewFocusMask->get_active());
}

void PreviewModePanel::togglebackColor0 ()
{
    backColor0->set_active(!backColor0->get_active());
}
void PreviewModePanel::togglebackColor1 ()
{
    backColor1->set_active(!backColor1->get_active());
}
void PreviewModePanel::togglebackColor2 ()
{
    backColor2->set_active(!backColor2->get_active());
}
void PreviewModePanel::togglebackColor3 ()
{
    backColor3->set_active(!backColor3->get_active());
}

void PreviewModePanel::buttonToggled (Gtk::ToggleButton* tbpreview)
{

    connR.block(true);
    connG.block(true);
    connB.block(true);
    connL.block(true);
    connFocusMask.block(true);

    // control state of the buttons
    // only 0 or 1 button at a time can remain pressed
    if (tbpreview != previewR) {
        previewR->set_active(false);
    }

    if (tbpreview != previewG) {
        previewG->set_active(false);
    }

    if (tbpreview != previewB) {
        previewB->set_active(false);
    }

    if (tbpreview != previewL) {
        previewL->set_active(false);
    }

    if (tbpreview != previewFocusMask) {
        previewFocusMask->set_active(false);
    }

    // set image based on button's state
    previewR->set_image(previewR->get_active() ? *iR : *igR);
    previewG->set_image(previewG->get_active() ? *iG : *igG);
    previewB->set_image(previewB->get_active() ? *iB : *igB);
    previewL->set_image(previewL->get_active() ? *iL : *igL);
    previewFocusMask->set_image(previewFocusMask->get_active() ? *iF : *igF);

    connR.block(false);
    connG.block(false);
    connB.block(false);
    connL.block(false);
    connFocusMask.block(false);

    imageArea->queue_draw ();

    // this will redraw the linked Before image area
    // which is set when before/after view is enabled
    if (imageArea->iLinkedImageArea != nullptr) {
        imageArea->iLinkedImageArea->queue_draw ();
    }
}

int PreviewModePanel::GetbackColor()
{
    int backColor = 0;

    if (backColor0->get_active ()) {
        backColor = 0;
    }

    if (backColor1->get_active ()) {
        backColor = 1;
    }

    if (backColor2->get_active ()) {
        backColor = 2;
    }

    if (backColor3->get_active ()) {
        backColor = 3;
    }

    return backColor;
}

void PreviewModePanel::togglebackColor()
{
    int backColor = GetbackColor();

    if(backColor == 0) {
        togglebackColor1();
    } else if(backColor == 1) {
        togglebackColor3();
    } else if(backColor == 3) {
        togglebackColor2();
    } else {
        togglebackColor0();
    }
}

void PreviewModePanel::buttonToggled_backColor (Gtk::ToggleButton* tbbackColor)
{

    connbackColor0.block(true);
    connbackColor1.block(true);
    connbackColor2.block(true);
    connbackColor3.block(true);

    // control the state of the buttons
    // Exactly 1 button at a time must remain pressed
    if (tbbackColor == backColor0 && !backColor0->get_active()) {
        backColor0->set_active(true);
    }

    if (tbbackColor == backColor1 && !backColor1->get_active()) {
        backColor1->set_active(true);
    }

    if (tbbackColor == backColor2 && !backColor2->get_active()) {
        backColor2->set_active(true);
    }

    if (tbbackColor == backColor3 && !backColor3->get_active()) {
        backColor3->set_active(true);
    }

    if (tbbackColor != backColor0) {
        backColor0->set_active(false);
    }

    if (tbbackColor != backColor1) {
        backColor1->set_active(false);
    }

    if (tbbackColor != backColor2) {
        backColor2->set_active(false);
    }

    if (tbbackColor != backColor3) {
        backColor3->set_active(false);
    }

    // set image based on button's state
    backColor0->set_image(backColor0->get_active() ? *iBC0 : *igBC0);
    backColor1->set_image(backColor1->get_active() ? *iBC1 : *igBC1);
    backColor2->set_image(backColor2->get_active() ? *iBC2 : *igBC2);
    backColor3->set_image(backColor3->get_active() ? *iBC3 : *igBC3);

    connbackColor0.block(false);
    connbackColor1.block(false);
    connbackColor2.block(false);
    connbackColor3.block(false);

    //TODO not sure if queue_draw is necessary, but will need to reach to backColor of the Before view
    imageArea->queue_draw ();

    // this will redraw the linked Before image area
    // which is set when before/after view is enabled
    if (imageArea->iLinkedImageArea != nullptr) {
        imageArea->iLinkedImageArea->queue_draw ();
    }
}
