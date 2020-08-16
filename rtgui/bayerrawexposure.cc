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
#include "bayerrawexposure.h"

#include "guiutils.h"
#include "options.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

BayerRAWExposure::BayerRAWExposure () : FoldableToolPanel(this, "bayerrawexposure", M("TP_EXPOS_BLACKPOINT_LABEL"), options.prevdemo != PD_Sidecar)
{
    PexBlack1 = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACK_1"), -2048, 2048, 1.0, 0)); //black level
    PexBlack1->setAdjusterListener (this);

    PexBlack1->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    PexBlack1->show();
    PexBlack2 = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACK_2"), -2048, 2048, 1.0, 0)); //black level
    PexBlack2->setAdjusterListener (this);

    PexBlack2->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    PexBlack2->show();
    PexBlack3 = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACK_3"), -2048, 2048, 1.0, 0)); //black level
    PexBlack3->setAdjusterListener (this);

    PexBlack3->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    PexBlack3->show();
    PexBlack0 = Gtk::manage(new Adjuster (M("TP_RAWEXPOS_BLACK_0"), -2048, 2048, 1.0, 0)); //black level
    PexBlack0->setAdjusterListener (this);

    PexBlack0->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    PexBlack0->show();
    PextwoGreen = Gtk::manage(new CheckBox(M("TP_RAWEXPOS_TWOGREEN"), multiImage));// two green
    PextwoGreen->set_active (true);
    PextwoGreen->setCheckBoxListener (this);

    pack_start( *PexBlack1, Gtk::PACK_SHRINK, 0);//black R
    pack_start( *PexBlack0, Gtk::PACK_SHRINK, 0);//black G1
    pack_start( *PexBlack3, Gtk::PACK_SHRINK, 0);//black G2
    pack_start( *PexBlack2, Gtk::PACK_SHRINK, 0);//black B
    pack_start( *PextwoGreen, Gtk::PACK_SHRINK, 0);//black 2 green

    PexBlack0->setLogScale(100, 0);
    PexBlack1->setLogScale(100, 0);
    PexBlack2->setLogScale(100, 0);
    PexBlack3->setLogScale(100, 0);
}

void BayerRAWExposure::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if(pedited ) {
        PexBlack0->setEditedState( pedited->raw.bayersensor.exBlack0 ? Edited : UnEdited );
        PexBlack1->setEditedState( pedited->raw.bayersensor.exBlack1 ? Edited : UnEdited );
        PexBlack2->setEditedState( pedited->raw.bayersensor.exBlack2 ? Edited : UnEdited );
        PexBlack3->setEditedState( pedited->raw.bayersensor.exBlack3 ? Edited : UnEdited );
    }

    PextwoGreen->setValue (pp->raw.bayersensor.twogreen);

    PexBlack0->setValue (pp->raw.bayersensor.black0);//black
    PexBlack1->setValue (pp->raw.bayersensor.black1);//black
    PexBlack2->setValue (pp->raw.bayersensor.black2);//black

    if(!PextwoGreen->getLastActive()) {
        PexBlack3->setValue (pp->raw.bayersensor.black3);
    } else {
        PexBlack3->setValue (PexBlack0->getValue());
    }

    enableListener ();
}

void BayerRAWExposure::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->raw.bayersensor.black0 = PexBlack0->getValue();// black
    pp->raw.bayersensor.black1 = PexBlack1->getValue();// black
    pp->raw.bayersensor.black2 = PexBlack2->getValue();// black
    pp->raw.bayersensor.twogreen = PextwoGreen->getLastActive();

    if(PextwoGreen->getLastActive()) {
        pp->raw.bayersensor.black3 = pp->raw.bayersensor.black0;   // active or desactive 2 green together
    } else {
        pp->raw.bayersensor.black3 = PexBlack3->getValue();
    }

    if (pedited) {
        pedited->raw.bayersensor.exBlack0 = PexBlack0->getEditedState ();//black
        pedited->raw.bayersensor.exBlack1 = PexBlack1->getEditedState ();//black
        pedited->raw.bayersensor.exBlack2 = PexBlack2->getEditedState ();//black
        pedited->raw.bayersensor.exBlack3 = PexBlack3->getEditedState ();//black
        pedited->raw.bayersensor.exTwoGreen = !PextwoGreen->get_inconsistent();
    }

}

void BayerRAWExposure::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        Glib::ustring value = a->getTextValue();

        if (a == PexBlack0) {
            if(!PextwoGreen->getLastActive()) {
                listener->panelChanged (EvPreProcessExpBlackzero,  value );
            } else {
                listener->panelChanged (EvPreProcessExpBlackzero,  value );
                PexBlack3->setValue (PexBlack0->getValue());
            }
        } else if (a == PexBlack1) {
            listener->panelChanged (EvPreProcessExpBlackone,  value );
        } else if (a == PexBlack2) {
            listener->panelChanged (EvPreProcessExpBlacktwo,  value );
        } else if (a == PexBlack3)    {
            if(!PextwoGreen->getLastActive()) {
                listener->panelChanged (EvPreProcessExpBlackthree,  value );
            } else {
                listener->panelChanged (EvPreProcessExpBlackthree,  value );
                PexBlack0->setValue (PexBlack3->getValue());
            }
        }
    }
}

void BayerRAWExposure::checkBoxToggled (CheckBox* c, CheckValue newval)
{
    if (c == PextwoGreen) {
        if (listener) {
            listener->panelChanged (EvPreProcessExptwoGreen, PextwoGreen->getLastActive() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
            if (PextwoGreen->getLastActive()) {
                PexBlack3->setValue (PexBlack0->getValue());//two green together
            }
        }
    }
}

void BayerRAWExposure::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    PexBlack0->showEditedCB ();//black
    PexBlack1->showEditedCB ();//black
    PexBlack2->showEditedCB ();//black
    PexBlack3->showEditedCB ();//black

}

void BayerRAWExposure::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    PexBlack0->setDefault( defParams->raw.bayersensor.black0);
    PexBlack1->setDefault( defParams->raw.bayersensor.black1);
    PexBlack2->setDefault( defParams->raw.bayersensor.black2);
    PexBlack3->setDefault( defParams->raw.bayersensor.black3);

    if (pedited) {
        PexBlack0->setDefaultEditedState( pedited->raw.bayersensor.exBlack0 ? Edited : UnEdited);
        PexBlack1->setDefaultEditedState( pedited->raw.bayersensor.exBlack1 ? Edited : UnEdited);
        PexBlack2->setDefaultEditedState( pedited->raw.bayersensor.exBlack2 ? Edited : UnEdited);
        PexBlack3->setDefaultEditedState( pedited->raw.bayersensor.exBlack3 ? Edited : UnEdited);

    } else {
        PexBlack0->setDefaultEditedState( Irrelevant );
        PexBlack1->setDefaultEditedState( Irrelevant );
        PexBlack2->setDefaultEditedState( Irrelevant );
        PexBlack3->setDefaultEditedState( Irrelevant );

    }
}

void BayerRAWExposure::setAdjusterBehavior (bool pexblackadd)
{

    PexBlack0->setAddMode(pexblackadd);
    PexBlack1->setAddMode(pexblackadd);
    PexBlack2->setAddMode(pexblackadd);
    PexBlack3->setAddMode(pexblackadd);
}

void BayerRAWExposure::trimValues (rtengine::procparams::ProcParams* pp)
{

    PexBlack0->trimValue(pp->raw.bayersensor.black0);
    PexBlack1->trimValue(pp->raw.bayersensor.black1);
    PexBlack2->trimValue(pp->raw.bayersensor.black2);
    PexBlack3->trimValue(pp->raw.bayersensor.black3);
}
