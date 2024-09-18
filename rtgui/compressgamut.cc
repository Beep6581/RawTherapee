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
#include "compressgamut.h"

#include "eventmapper.h"

#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring Compressgamut::TOOL_NAME = "compressgamut";

Compressgamut::Compressgamut () : FoldableToolPanel(this, TOOL_NAME, M("TP_COMPRESSGAMUT_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvcgColorspace = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_CG_COLORSPACE");
    
    Gtk::Box* hb = Gtk::manage (new Gtk::Box ());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_COMPRESSGAMUT_MAIN_COLORSPACE") + ": ")), Gtk::PACK_SHRINK);
    colorspace = Gtk::manage(new MyComboBoxText());
    colorspace->append(M("TP_COMPRESSGAMUT_REC2020"));
    colorspace->append(M("TP_COMPRESSGAMUT_PROPHOTO"));
    colorspace->append(M("TP_COMPRESSGAMUT_SRGB"));
    colorspace->append(M("TP_COMPRESSGAMUT_DCIP3"));
    colorspace->append(M("TP_COMPRESSGAMUT_ACESP1"));
    hb->pack_start(*colorspace);
    pack_start(*hb);

    pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL)));

    th_c = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_CYANTH"), 0., 0.999, 0.001, 0.815));
    th_m = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_MAGENTATH"), 0., 0.999, 0.001, 0.803));
    th_y = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_YELLOWTH"), 0., 0.999, 0.001, 0.880));

    Gtk::Frame *thFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_THRESHOLD")));
    thFrame->set_label_align(0.025, 0.5);
    Gtk::Box *thVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    thVBox->pack_start (*th_c);
    thVBox->pack_start (*th_m);
    thVBox->pack_start (*th_y);
    thFrame->add(*thVBox);
    pack_start(*thFrame, Gtk::PACK_EXPAND_WIDGET);

    d_c = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_CYANLIM"), 1.001, 2.0, 0.001, 1.147));
    d_m = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_MAGENTALIM"), 1.001, 2.0, 0.001, 1.264));
    d_y = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_YELLOWLIM"), 1.001, 2.0, 0.001, 1.312));
    
    Gtk::Frame *limFrame = Gtk::manage(new Gtk::Frame(M("TP_COMPRESSGAMUT_LIMIT")));
    limFrame->set_label_align(0.025, 0.5);
    Gtk::Box *limVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    limVBox->pack_start (*d_c);
    limVBox->pack_start (*d_m);
    limVBox->pack_start (*d_y);
    limFrame->add(*limVBox);
    pack_start(*limFrame, Gtk::PACK_EXPAND_WIDGET);


    rolloff = Gtk::manage(new Gtk::CheckButton(M("TP_COMPRESSGAMUT_ROLLOFF")));
    pwr = Gtk::manage (new Adjuster (M("TP_COMPRESSGAMUT_PWR"), 0.5, 2.0, 0.01, 1.2));
    rolloffconn = rolloff->signal_pressed().connect ( sigc::mem_fun (*this, &Compressgamut::rolloff_change) );

    Gtk::Frame *rollFrame = Gtk::manage(new Gtk::Frame());
    rollFrame->set_label_align(0.025, 0.5);
    rollFrame->set_label_widget(*rolloff);
    ToolParamBlock* const rollvBox = Gtk::manage(new ToolParamBlock());
    rollvBox->pack_start (*pwr);
    rollFrame->add(*rollvBox);
    pack_start(*rollFrame, Gtk::PACK_EXPAND_WIDGET);
    
/*

    colorspace->signal_changed().connect(sigc::mem_fun(*this, &ShadowsHighlights::colorspaceChanged));
    
    show_all_children ();
    */
}

void Compressgamut::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();
    /*
    if (pedited) {
        radius->setEditedState       (pedited->sh.radius ? Edited : UnEdited);
        highlights->setEditedState   (pedited->sh.highlights ? Edited : UnEdited);
        h_tonalwidth->setEditedState (pedited->sh.htonalwidth ? Edited : UnEdited);
        shadows->setEditedState      (pedited->sh.shadows ? Edited : UnEdited);
        s_tonalwidth->setEditedState (pedited->sh.stonalwidth ? Edited : UnEdited);
        set_inconsistent             (multiImage && !pedited->sh.enabled);

    }

    setEnabled (pp->sh.enabled);

    radius->setValue        (pp->sh.radius);
    highlights->setValue    (pp->sh.highlights);
    h_tonalwidth->setValue  (pp->sh.htonalwidth);
    shadows->setValue       (pp->sh.shadows);
    s_tonalwidth->setValue  (pp->sh.stonalwidth);

    if (pedited && !pedited->sh.lab) {
        colorspace->set_active(2);
    } else if (pp->sh.lab) {
        colorspace->set_active(1);
    } else {
        colorspace->set_active(0);
    }
    */
    enableListener ();
}

void Compressgamut::write (ProcParams* pp, ParamsEdited* pedited)
{
    /*
    pp->sh.radius        = (int)radius->getValue ();
    pp->sh.highlights    = (int)highlights->getValue ();
    pp->sh.htonalwidth   = (int)h_tonalwidth->getValue ();
    pp->sh.shadows       = (int)shadows->getValue ();
    pp->sh.stonalwidth   = (int)s_tonalwidth->getValue ();
    pp->sh.enabled       = getEnabled();

    if (colorspace->get_active_row_number() == 0) {
        pp->sh.lab = false;
    } else if (colorspace->get_active_row_number() == 1) {
        pp->sh.lab = true;
    }

    if (pedited) {
        pedited->sh.radius          = radius->getEditedState ();
        pedited->sh.highlights      = highlights->getEditedState ();
        pedited->sh.htonalwidth     = h_tonalwidth->getEditedState ();
        pedited->sh.shadows         = shadows->getEditedState ();
        pedited->sh.stonalwidth     = s_tonalwidth->getEditedState ();
        pedited->sh.enabled         = !get_inconsistent();
        pedited->sh.lab = colorspace->get_active_row_number() != 2;
    }
    */
}

void Compressgamut::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
/*
    radius->setDefault (defParams->sh.radius);
    highlights->setDefault (defParams->sh.highlights);
    h_tonalwidth->setDefault (defParams->sh.htonalwidth);
    shadows->setDefault (defParams->sh.shadows);
    s_tonalwidth->setDefault (defParams->sh.stonalwidth);

    if (pedited) {
        radius->setDefaultEditedState       (pedited->sh.radius ? Edited : UnEdited);
        highlights->setDefaultEditedState   (pedited->sh.highlights ? Edited : UnEdited);
        h_tonalwidth->setDefaultEditedState (pedited->sh.htonalwidth ? Edited : UnEdited);
        shadows->setDefaultEditedState      (pedited->sh.shadows ? Edited : UnEdited);
        s_tonalwidth->setDefaultEditedState (pedited->sh.stonalwidth ? Edited : UnEdited);
    } else {
        radius->setDefaultEditedState       (Irrelevant);
        highlights->setDefaultEditedState   (Irrelevant);
        h_tonalwidth->setDefaultEditedState (Irrelevant);
        shadows->setDefaultEditedState      (Irrelevant);
        s_tonalwidth->setDefaultEditedState (Irrelevant);
    }
    */
}

void Compressgamut::adjusterChanged (Adjuster* a, double newval)
{
    /*
    if (listener && getEnabled()) {
        const Glib::ustring costr = Glib::ustring::format ((int)a->getValue());

        if (a == highlights) {
            listener->panelChanged (EvSHHighlights, costr);
        } else if (a == h_tonalwidth) {
            listener->panelChanged (EvSHHLTonalW, costr);
        } else if (a == shadows) {
            listener->panelChanged (EvSHShadows, costr);
        } else if (a == s_tonalwidth) {
            listener->panelChanged (EvSHSHTonalW, costr);
        } else if (a == radius) {
            listener->panelChanged (EvSHRadius, costr);
        }
    }
    */
}

void Compressgamut::rolloff_change ()
{
}


void Compressgamut::enabledChanged ()
{
    /*
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvSHEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvSHEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvSHEnabled, M("GENERAL_DISABLED"));
        }
    }
    */
}
/*
void ShadowsHighlights::colorspaceChanged()
{
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged(EvSHColorspace, colorspace->get_active_text());
    }
}
*/
void Compressgamut::setBatchMode (bool batchMode)
{
    /*
    ToolPanel::setBatchMode (batchMode);
    radius->showEditedCB ();
    highlights->showEditedCB ();
    h_tonalwidth->showEditedCB ();
    shadows->showEditedCB ();
    s_tonalwidth->showEditedCB ();
    colorspace->append(M("GENERAL_UNCHANGED"));    
    */
}
/*
void ShadowsHighlights::setAdjusterBehavior (bool hadd, bool sadd)
{

    highlights->setAddMode(hadd);
    shadows->setAddMode(sadd);
}
*/
void Compressgamut::trimValues (rtengine::procparams::ProcParams* pp)
{
    /*
    highlights->trimValue(pp->sh.highlights);
    shadows->trimValue(pp->sh.shadows);
    */
}
