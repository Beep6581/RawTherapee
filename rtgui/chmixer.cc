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
#include "chmixer.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

ChMixer::ChMixer (): FoldableToolPanel(this, "chmixer", M("TP_CHMIXER_LABEL"), false, true)
{

    imgIcon[0] = Gtk::manage (new RTImage ("Chanmixer-RR.png"));
    imgIcon[1] = Gtk::manage (new RTImage ("Chanmixer-RG.png"));
    imgIcon[2] = Gtk::manage (new RTImage ("Chanmixer-RB.png"));
    imgIcon[3] = Gtk::manage (new RTImage ("Chanmixer-GR.png"));
    imgIcon[4] = Gtk::manage (new RTImage ("Chanmixer-GG.png"));
    imgIcon[5] = Gtk::manage (new RTImage ("Chanmixer-GB.png"));
    imgIcon[6] = Gtk::manage (new RTImage ("Chanmixer-BR.png"));
    imgIcon[7] = Gtk::manage (new RTImage ("Chanmixer-BG.png"));
    imgIcon[8] = Gtk::manage (new RTImage ("Chanmixer-BB.png"));

    Gtk::Label* rlabel = Gtk::manage (new Gtk::Label ());
    rlabel->set_markup (Glib::ustring("\t<span foreground=\"#b00000\"><b>") + M("TP_CHMIXER_RED") + Glib::ustring(":</b></span>"));
    rlabel->set_alignment(Gtk::ALIGN_START);

    red[0] = Gtk::manage (new Adjuster ("",   -200, 200, 1, 100, imgIcon[0]));
    red[1] = Gtk::manage (new Adjuster ("", -200, 200, 1, 0, imgIcon[1]));
    red[2] = Gtk::manage (new Adjuster ("",  -200, 200, 1, 0, imgIcon[2]));

    Gtk::HSeparator* rsep = Gtk::manage (new Gtk::HSeparator ());

    pack_start (*rlabel);

    for (int i = 0; i < 3; i++) {
        pack_start (*red[i]);
    }

    pack_start (*rsep, Gtk::PACK_EXPAND_WIDGET, 4);

    Gtk::Label* glabel = Gtk::manage (new Gtk::Label ());
    glabel->set_markup (Glib::ustring("\t<span foreground=\"#0b8c21\"><b>") + M("TP_CHMIXER_GREEN") + Glib::ustring(":</b></span>"));
    glabel->set_alignment(Gtk::ALIGN_START);


    green[0] = Gtk::manage (new Adjuster ("",   -200, 200, 1, 0, imgIcon[3]));
    green[1] = Gtk::manage (new Adjuster ("", -200, 200, 1, 100, imgIcon[4]));
    green[2] = Gtk::manage (new Adjuster ("",  -200, 200, 1, 0, imgIcon[5]));

    Gtk::HSeparator* gsep = Gtk::manage (new Gtk::HSeparator ());

    pack_start (*glabel);

    for (int i = 0; i < 3; i++) {
        pack_start (*green[i]);
    }

    pack_start (*gsep, Gtk::PACK_EXPAND_WIDGET, 4);

    Gtk::Label* blabel = Gtk::manage (new Gtk::Label ());
    blabel->set_markup (Glib::ustring("\t<span foreground=\"#1377d7\"><b>") + M("TP_CHMIXER_BLUE") + Glib::ustring(":</b></span>"));
    blabel->set_alignment(Gtk::ALIGN_START);
    blue[0] = Gtk::manage (new Adjuster ("",   -200, 200, 1, 0, imgIcon[6]));
    blue[1] = Gtk::manage (new Adjuster ("", -200, 200, 1, 0, imgIcon[7]));
    blue[2] = Gtk::manage (new Adjuster ("",  -200, 200, 1, 100, imgIcon[8]));

    for (int i = 0; i < 3; i++) {
        red[i]->setAdjusterListener (this);
        green[i]->setAdjusterListener (this);
        blue[i]->setAdjusterListener (this);
    }

    pack_start (*blabel);

    for (int i = 0; i < 3; i++) {
        pack_start (*blue[i]);
    }

    show_all();
}

void ChMixer::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    setEnabled(pp->chmixer.enabled);
    
    if (pedited) {
        for (int i = 0; i < 3; i++) {
            red[i]->setEditedState (pedited->chmixer.red[i] ? Edited : UnEdited);
            green[i]->setEditedState (pedited->chmixer.green[i] ? Edited : UnEdited);
            blue[i]->setEditedState (pedited->chmixer.blue[i] ? Edited : UnEdited);
        }
        set_inconsistent(multiImage && !pedited->chmixer.enabled);
    }

    for (int i = 0; i < 3; i++) {
        red[i]->setValue (pp->chmixer.red[i]);
        green[i]->setValue (pp->chmixer.green[i]);
        blue[i]->setValue (pp->chmixer.blue[i]);
    }

    enableListener ();
}

void ChMixer::write (ProcParams* pp, ParamsEdited* pedited)
{

    for (int i = 0; i < 3; i++) {
        pp->chmixer.red[i] = (int) red[i]->getValue ();
        pp->chmixer.green[i] = (int) green[i]->getValue ();
        pp->chmixer.blue[i] = (int) blue[i]->getValue ();
    }
    pp->chmixer.enabled = getEnabled();

    if (pedited) {
        for (int i = 0; i < 3; i++) {
            pedited->chmixer.red[i] = red[i]->getEditedState ();
            pedited->chmixer.green[i] = green[i]->getEditedState ();
            pedited->chmixer.blue[i] = blue[i]->getEditedState ();
        }
        pedited->chmixer.enabled = !get_inconsistent();
    }
}

void ChMixer::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    for (int i = 0; i < 3; i++) {
        red[i]->setDefault (defParams->chmixer.red[i]);
        green[i]->setDefault (defParams->chmixer.green[i]);
        blue[i]->setDefault (defParams->chmixer.blue[i]);
    }

    if (pedited)
        for (int i = 0; i < 3; i++) {
            red[i]->setDefaultEditedState (pedited->chmixer.red[i] ? Edited : UnEdited);
            green[i]->setDefaultEditedState (pedited->chmixer.green[i] ? Edited : UnEdited);
            blue[i]->setDefaultEditedState (pedited->chmixer.blue[i] ? Edited : UnEdited);
        }
    else
        for (int i = 0; i < 3; i++) {
            red[i]->setDefaultEditedState (Irrelevant);
            green[i]->setDefaultEditedState (Irrelevant);
            blue[i]->setDefaultEditedState (Irrelevant);
        }
}

void ChMixer::adjusterChanged (Adjuster* a, double newval)
{

    if (listener && getEnabled()) {
        Glib::ustring descr = Glib::ustring::compose ("R=%1,%2,%3\nG=%4,%5,%6\nB=%7,%8,%9",
                              (int)red[0]->getValue(), (int)red[1]->getValue(), (int)red[2]->getValue(),
                              (int)green[0]->getValue(), (int)green[1]->getValue(), (int)green[2]->getValue(),
                              (int)blue[0]->getValue(), (int)blue[1]->getValue(), (int)blue[2]->getValue());
        listener->panelChanged (EvChMixer, descr);
    }
}


void ChMixer::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvChMixer, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvChMixer, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvChMixer, M("GENERAL_DISABLED"));
        }
    }
}


void ChMixer::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);

    for (int i = 0; i < 3; i++) {
        red[i]->showEditedCB ();
        green[i]->showEditedCB ();
        blue[i]->showEditedCB ();
    }
}

void ChMixer::setAdjusterBehavior (bool rgbadd)
{

    for (int i = 0; i < 3; i++) {
        red[i]->setAddMode(rgbadd);
        green[i]->setAddMode(rgbadd);
        blue[i]->setAddMode(rgbadd);
    }
}

void ChMixer::trimValues (rtengine::procparams::ProcParams* pp)
{

    for (int i = 0; i < 3; i++) {
        red[i]->trimValue(pp->chmixer.red[i]);
        green[i]->trimValue(pp->chmixer.green[i]);
        blue[i]->trimValue(pp->chmixer.blue[i]);
    }
}
