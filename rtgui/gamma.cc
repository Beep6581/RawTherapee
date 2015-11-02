/*
 *  This file is part of RawTherapee.
 */
#include "gamma.h"
#include "mycurve.h"

using namespace rtengine;
using namespace rtengine::procparams;

Gamma::Gamma () : FoldableToolPanel(this, "gamma", M("TP_GAMMADIF_LABEL"), false, true)
{

    Gtk::VBox * gammaVBox = Gtk::manage ( new Gtk::VBox());
    gammaVBox->set_border_width(4);
    gammaVBox->set_spacing(2);


    gabox = Gtk::manage (new Gtk::HBox ());
    labmga = Gtk::manage (new Gtk::Label (M("TP_GAMMADIF_METHOD") + ":"));
    gabox->pack_start (*labmga, Gtk::PACK_SHRINK, 1);

    gammaMethod = Gtk::manage (new MyComboBoxText ());
    gammaMethod->append_text (M("TP_GAMMADIF_ONE"));
    gammaMethod->append_text (M("TP_GAMMADIF_ONEABS"));
    gammaMethod->append_text (M("TP_GAMMADIF_TWO"));
    gammaMethod->append_text (M("TP_GAMMADIF_THR"));
    gammaMethod->set_active(0);
    gammaMethodConn = gammaMethod->signal_changed().connect ( sigc::mem_fun(*this, &Gamma::gammaMethodChanged) );
    gammaMethod->set_tooltip_markup (M("TP_GAMMADIF_METHOD_TOOLTIP"));

    gabox->pack_start (*gammaMethod);
    gammaVBox->pack_start(*gabox);

    gamm = Gtk::manage (new Adjuster (M("TP_GAMMADIF_GAMMA"), 0.5, 3., 0.01, 1.));
    slop = Gtk::manage (new Adjuster (M("TP_GAMMADIF_SLOPE"), 0., 20., 0.01, 2.));

    outp = Gtk::manage(new Gtk::CheckButton((M("TP_GAMMADIF_OUTP"))));
    outp->set_active (false);
    gammaVBox->pack_start(*outp, Gtk::PACK_SHRINK, 0);
    outpconn = outp->signal_toggled().connect ( sigc::mem_fun(*this, &Gamma::outpChanged));

    gammaVBox->pack_start (*gamm);
    gamm->show ();

    gammaVBox->pack_start (*slop);
    slop->show ();


    gamm->setAdjusterListener (this);
    slop->setAdjusterListener (this);


    pack_start (*gammaVBox);

    disableListener();
    enableListener();

}

Gamma::~Gamma()
{

}




void Gamma::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    gammaMethodConn.block(true);


    if (pedited) {
        gamm->setEditedState (pedited->gamma.gamm ? Edited : UnEdited);
        slop->setEditedState (pedited->gamma.slop ? Edited : UnEdited);


        if (!pedited->gamma.gammaMethod) {
            gammaMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

    }

    gamm->setValue    (pp->gamma.gamm);
    slop->setValue      (pp->gamma.slop);
    outp->set_active (pp->gamma.outp);
    lastoutp = pp->gamma.outp;

    setEnabled (pp->gamma.enabled);


    if (pp->gamma.gammaMethod == "one") {
        gammaMethod->set_active (0);
    } else if (pp->gamma.gammaMethod == "oneabs") {
        gammaMethod->set_active (1);
    } else if (pp->gamma.gammaMethod == "two") {
        gammaMethod->set_active (2);
    } else if (pp->gamma.gammaMethod == "thr") {
        gammaMethod->set_active (3);
    }


    gammaMethodChanged ();


    gammaMethodConn.block(false);


    enableListener ();
}



void Gamma::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->gamma.gamm    = gamm->getValue ();
    pp->gamma.slop    = slop->getValue ();
    pp->gamma.enabled      = getEnabled();
    pp->gamma.outp = outp->get_active();

    if (pedited) {
        pedited->gamma.gammaMethod    = gammaMethod->get_active_text() != M("GENERAL_UNCHANGED");

        //%%%%%%%%%%%%%%%%%%%%%%
        pedited->gamma.gamm     = gamm->getEditedState ();
        pedited->gamma.slop     = slop->getEditedState ();
        pedited->gamma.outp = !outp->get_inconsistent();

    }

    if (gammaMethod->get_active_row_number() == 0) {
        pp->gamma.gammaMethod = "one";
    } else if (gammaMethod->get_active_row_number() == 1) {
        pp->gamma.gammaMethod = "oneabs";
    } else if (gammaMethod->get_active_row_number() == 2) {
        pp->gamma.gammaMethod = "two";
    } else if (gammaMethod->get_active_row_number() == 3) {
        pp->gamma.gammaMethod = "thr";
    }



}

void Gamma::outpChanged()
{
    if (batchMode) {
        if (outp->get_inconsistent()) {
            outp->set_inconsistent (false);
            outpconn.block (true);
            outp->set_active (false);
            outpconn.block (false);
        } else if (lastoutp) {
            outp->set_inconsistent (true);
        }

        lastoutp = outp->get_active ();
    }
    if (outp->get_active()) {
        gamm->set_sensitive(false);
        slop->set_sensitive(false);
    }
    else {
        gamm->set_sensitive(true);
        slop->set_sensitive(true);
    }

    if (listener) {
        if (outp->get_active()) {
            listener->panelChanged (Evoutp, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evoutp, M("GENERAL_DISABLED"));
        }
    }
}

void Gamma::gammaMethodChanged()
{

    if (gammaMethod->get_active_row_number() == 1) outp->set_sensitive(true);
    else outp->set_sensitive(false);

    if (!listener || !getEnabled()) {
        return;
    }

    if (listener) {
        listener->panelChanged (EvgammaMethod, gammaMethod->get_active_text ());
    }
}




void Gamma::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    gamm->setDefault (defParams->gamma.gamm);
    slop->setDefault (defParams->gamma.slop);

    if (pedited) {
        gamm->setDefaultEditedState (pedited->gamma.gamm ? Edited : UnEdited);
        slop->setDefaultEditedState (pedited->gamma.slop ? Edited : UnEdited);

    } else {
        gamm->setDefaultEditedState (Irrelevant);
        slop->setDefaultEditedState (Irrelevant);
    }
}

void Gamma::setAdjusterBehavior (bool gammAdd, bool slopeAdd)
{
    gamm->setAddMode(gammAdd);
    slop->setAddMode(slopeAdd);
}


void Gamma::adjusterChanged (Adjuster* a, double newval)
{

    if (!listener || !getEnabled()) {
        return;
    }

    if (a == gamm) {
        listener->panelChanged (EvGgamm, gamm->getTextValue());
    } else if (a == slop) {
        listener->panelChanged (EvGslop, slop->getTextValue());
    }

}


void Gamma::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvGammaEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvGammaEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvGammaEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void Gamma::trimValues (rtengine::procparams::ProcParams* pp)
{
    gamm->trimValue(pp->gamma.gamm);
    slop->trimValue(pp->gamma.slop);

}

void Gamma::setBatchMode (bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    gamm->showEditedCB ();
    slop->showEditedCB ();
}
