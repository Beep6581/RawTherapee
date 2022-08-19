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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "previewmodepanel.h"
#include "options.h"
#include "multilangmgr.h"
#include "imagearea.h"
#include "rtimage.h"

PreviewModePanel::PreviewModePanel (ImageArea* ia) :
    imageArea(ia),
    nR("square-toggle-red-on-narrow"), ngR("square-toggle-red-off-narrow"),
    nG("square-toggle-green-on-narrow"), ngG("square-toggle-green-off-narrow"),
    nB("square-toggle-blue-on-narrow"), ngB("square-toggle-blue-off-narrow"),
    nL("square-toggle-luminosity-on-narrow"), ngL("square-toggle-luminosity-off-narrow"),
    nBC0("square-toggle-theme-on-narrow"), ngBC0("square-toggle-theme-off-narrow"),
    nBC1("square-toggle-black-on-narrow"), ngBC1("square-toggle-black-off-narrow"),
    nBC2("square-toggle-white-on-narrow"), ngBC2("square-toggle-white-off-narrow"),
    nBC3("square-toggle-gray-on-narrow"), ngBC3("square-toggle-gray-off-narrow"),
    iR(Gtk::manage (new RTImage(ngR, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    iG(Gtk::manage (new RTImage(ngG, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    iB(Gtk::manage (new RTImage(ngB, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    iL(Gtk::manage (new RTImage(ngL, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    iBC0(Gtk::manage (new RTImage(options.bgcolor == 0 ? nBC0 : ngBC0, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    iBC1(Gtk::manage (new RTImage(options.bgcolor == 1 ? nBC1 : ngBC1, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    iBC2(Gtk::manage (new RTImage(options.bgcolor == 2 ? nBC2 : ngBC2, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    iBC3(Gtk::manage (new RTImage(options.bgcolor == 3 ? nBC3 : ngBC3, Gtk::ICON_SIZE_LARGE_TOOLBAR)))
{
    backColor0 = Gtk::manage (new Gtk::ToggleButton ());
    backColor0->get_style_context()->add_class("narrowbutton");
    backColor0->set_relief(Gtk::RELIEF_NONE);
    backColor0->set_tooltip_markup (M("MAIN_TOOLTIP_BACKCOLOR0"));
    backColor0->set_image(*iBC0);

    backColor1 = Gtk::manage (new Gtk::ToggleButton ());
    backColor1->get_style_context()->add_class("narrowbutton");
    backColor1->set_relief(Gtk::RELIEF_NONE);
    backColor1->set_tooltip_markup (M("MAIN_TOOLTIP_BACKCOLOR1"));
    backColor1->set_image(*iBC1);

    backColor3 = Gtk::manage (new Gtk::ToggleButton ());
    backColor3->get_style_context()->add_class("narrowbutton");
    backColor3->set_relief(Gtk::RELIEF_NONE);
    backColor3->set_tooltip_markup (M("MAIN_TOOLTIP_BACKCOLOR3"));
    backColor3->set_image(*iBC3);

    backColor2 = Gtk::manage (new Gtk::ToggleButton ());
    backColor2->get_style_context()->add_class("narrowbutton");
    backColor2->set_relief(Gtk::RELIEF_NONE);
    backColor2->set_tooltip_markup (M("MAIN_TOOLTIP_BACKCOLOR2"));
    backColor2->set_image(*iBC2);

    previewR = Gtk::manage (new Gtk::ToggleButton ());
    previewR->get_style_context()->add_class("narrowbutton");
    previewR->set_relief(Gtk::RELIEF_NONE);
    previewR->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWR"));
    previewR->set_image(*iR);

    previewG = Gtk::manage (new Gtk::ToggleButton ());
    previewG->get_style_context()->add_class("narrowbutton");
    previewG->set_relief(Gtk::RELIEF_NONE);
    previewG->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWG"));
    previewG->set_image(*iG);

    previewB = Gtk::manage (new Gtk::ToggleButton ());
    previewB->get_style_context()->add_class("narrowbutton");
    previewB->set_relief(Gtk::RELIEF_NONE);
    previewB->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWB"));
    previewB->set_image(*iB);

    previewL = Gtk::manage (new Gtk::ToggleButton ());
    previewL->get_style_context()->add_class("narrowbutton");
    previewL->set_relief(Gtk::RELIEF_NONE);
    previewL->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWL"));
    previewL->set_image(*iL);

    previewR->set_active (false);
    previewG->set_active (false);
    previewB->set_active (false);
    previewL->set_active (false);

    backColor0->set_active (options.bgcolor == 0);
    backColor1->set_active (options.bgcolor == 1);
    backColor2->set_active (options.bgcolor == 2);
    backColor3->set_active (options.bgcolor == 3);

    pack_start (*backColor0, Gtk::PACK_SHRINK, 0);
    pack_start (*backColor1, Gtk::PACK_SHRINK, 0);
    pack_start (*backColor3, Gtk::PACK_SHRINK, 0);
    pack_start (*backColor2, Gtk::PACK_SHRINK, 0);

    pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_VERTICAL)), Gtk::PACK_SHRINK, 2);

    pack_start (*previewR, Gtk::PACK_SHRINK, 0);
    pack_start (*previewG, Gtk::PACK_SHRINK, 0);
    pack_start (*previewB, Gtk::PACK_SHRINK, 0);
    pack_start (*previewL, Gtk::PACK_SHRINK, 0);

    connR = previewR->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewR) );
    connG = previewG->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewG) );
    connB = previewB->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewB) );
    connL = previewL->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled), previewL) );

    connbackColor0 = backColor0->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled_backColor), backColor0) );
    connbackColor1 = backColor1->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled_backColor), backColor1) );
    connbackColor2 = backColor2->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled_backColor), backColor2) );
    connbackColor3 = backColor3->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &PreviewModePanel::buttonToggled_backColor), backColor3) );

    //show_all ();
}

PreviewModePanel::~PreviewModePanel () {}

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

    // set image based on button's state
    iR->set_from_icon_name(previewR->get_active() ? nR : ngR);
    iG->set_from_icon_name(previewG->get_active() ? nG : ngG);
    iB->set_from_icon_name(previewB->get_active() ? nB : ngB);
    iL->set_from_icon_name(previewL->get_active() ? nL : ngL);

    connR.block(false);
    connG.block(false);
    connB.block(false);
    connL.block(false);

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
    iBC0->set_from_icon_name(backColor0->get_active() ? nBC0 : ngBC0);
    iBC1->set_from_icon_name(backColor1->get_active() ? nBC1 : ngBC1);
    iBC2->set_from_icon_name(backColor2->get_active() ? nBC2 : ngBC2);
    iBC3->set_from_icon_name(backColor3->get_active() ? nBC3 : ngBC3);

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
