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
#include "sharpening.h"
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

Sharpening::Sharpening () : FoldableToolPanel(this, "sharpening", M("TP_SHARPENING_LABEL"), true, true)
{

    setEnabledTooltipMarkup(M("TP_SHARPENING_TOOLTIP"));

    Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
    hb->set_border_width (4);
    hb->show ();
    Gtk::Label* ml = Gtk::manage (new Gtk::Label (M("TP_SHARPENING_METHOD") + ":"));
    ml->show ();
    method = Gtk::manage (new MyComboBoxText ());
    method->append_text (M("TP_SHARPENING_USM"));
    method->append_text (M("TP_SHARPENING_RLD"));
    method->show ();
    hb->pack_start(*ml, Gtk::PACK_SHRINK, 4);
    hb->pack_start(*method);
    pack_start (*hb);

    rld = new Gtk::VBox ();
    dradius = Gtk::manage (new Adjuster (M("TP_SHARPENING_EDRADIUS"), 0.5, 2.5, 0.01, 0.75));
    damount = Gtk::manage (new Adjuster (M("TP_SHARPENING_RLD_AMOUNT"), 0.0, 100, 1, 75));
    ddamping = Gtk::manage (new Adjuster (M("TP_SHARPENING_RLD_DAMPING"), 0, 100, 1, 20));
    diter = Gtk::manage (new Adjuster (M("TP_SHARPENING_RLD_ITERATIONS"), 5, 100, 1, 30));
    rld->pack_start (*dradius);
    rld->pack_start (*damount);
    rld->pack_start (*ddamping);
    rld->pack_start (*diter);
    dradius->show ();
    damount->show ();
    ddamping->show ();
    diter->show ();
    rld->show ();

    usm = new Gtk::VBox ();
    usm->show ();


    Gtk::HSeparator *hsep6a = Gtk::manage (new  Gtk::HSeparator());
    amount = Gtk::manage (new Adjuster (M("TP_SHARPENING_AMOUNT"), 1, 1000, 1, 200));
    radius = Gtk::manage (new Adjuster (M("TP_SHARPENING_RADIUS"), 0.3, 3, 0.01, 0.5));
    threshold = Gtk::manage (new ThresholdAdjuster (M("TP_SHARPENING_THRESHOLD"), 0., 2000., 20., 80., 2000., 1200., 0, false));
    threshold->setAdjusterListener (this);
    pack_start(*hsep6a, Gtk::PACK_SHRINK, 2);

    pack_start (*usm);

    usm->pack_start(*radius);
    usm->pack_start(*threshold);
    usm->pack_start(*amount);
    hsep6a->show ();
    radius->show ();
    threshold->show ();
    amount->show ();

    Gtk::HSeparator *hsep6 = Gtk::manage (new  Gtk::HSeparator());
    edgesonly = Gtk::manage (new Gtk::CheckButton (M("TP_SHARPENING_ONLYEDGES")));
    edgesonly->set_active (false);
    edgebox = new Gtk::VBox ();
    eradius = Gtk::manage (new Adjuster (M("TP_SHARPENING_EDRADIUS"), 0.5, 2.5, 0.1, 1.9));
    etolerance = Gtk::manage (new Adjuster (M("TP_SHARPENING_EDTOLERANCE"), 10, 10000, 100, 1000));
    usm->pack_start(*hsep6, Gtk::PACK_SHRINK, 2);
    usm->pack_start(*edgesonly);
    edgebox->pack_start(*eradius);
    edgebox->pack_start(*etolerance);
    edgebox->show ();
    edgebin = Gtk::manage (new Gtk::VBox ());
    usm->pack_start (*edgebin);
    edgebin->show ();
    hsep6->show();
    edgesonly->show();
    eradius->show();
    etolerance->show();

    Gtk::HSeparator *hsep6b = Gtk::manage (new  Gtk::HSeparator());
    halocontrol = Gtk::manage (new Gtk::CheckButton (M("TP_SHARPENING_HALOCONTROL")));
    halocontrol->set_active (false);
    hcbox = new Gtk::VBox ();
    hcamount = Gtk::manage (new Adjuster (M("TP_SHARPENING_HCAMOUNT"), 1, 100, 1, 75));
    usm->pack_start(*hsep6b, Gtk::PACK_SHRINK, 2);
    usm->pack_start(*halocontrol);
    hcbox->pack_start(*hcamount);
    hcbox->show ();
    hcbin = Gtk::manage (new Gtk::VBox ());
    usm->pack_start (*hcbin);
    hcbin->show ();
    hsep6b->show ();
    halocontrol->show ();
    hcamount->show ();

    dradius->setAdjusterListener (this);
    damount->setAdjusterListener (this);
    ddamping->setAdjusterListener (this);
    diter->setAdjusterListener (this);
    radius->setAdjusterListener (this);
    amount->setAdjusterListener (this);
    eradius->setAdjusterListener (this);
    etolerance->setAdjusterListener (this);
    hcamount->setAdjusterListener (this);

    edgebox->reference ();
    hcbox->reference ();
    usm->reference ();
    rld->reference ();

    eonlyConn = edgesonly->signal_toggled().connect( sigc::mem_fun(*this, &Sharpening::edgesonly_toggled) );
    hcConn    = halocontrol->signal_toggled().connect( sigc::mem_fun(*this, &Sharpening::halocontrol_toggled) );
    method->signal_changed().connect( sigc::mem_fun(*this, &Sharpening::method_changed) );
}

Sharpening::~Sharpening ()
{

    delete usm;
    delete rld;
    delete edgebox;
    delete hcbox;
}


void Sharpening::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        amount->setEditedState      (pedited->sharpening.amount ? Edited : UnEdited);
        radius->setEditedState      (pedited->sharpening.radius ? Edited : UnEdited);
        threshold->setEditedState   (pedited->sharpening.threshold ? Edited : UnEdited);
        eradius->setEditedState     (pedited->sharpening.edges_radius ? Edited : UnEdited);
        etolerance->setEditedState  (pedited->sharpening.edges_tolerance ? Edited : UnEdited);
        hcamount->setEditedState    (pedited->sharpening.halocontrol_amount ? Edited : UnEdited);
        damount->setEditedState     (pedited->sharpening.deconvamount ? Edited : UnEdited);
        dradius->setEditedState     (pedited->sharpening.deconvradius ? Edited : UnEdited);
        diter->setEditedState       (pedited->sharpening.deconviter ? Edited : UnEdited);
        ddamping->setEditedState    (pedited->sharpening.deconvdamping ? Edited : UnEdited);

        halocontrol->set_inconsistent   (multiImage && !pedited->sharpening.halocontrol);
        edgesonly->set_inconsistent     (multiImage && !pedited->sharpening.edgesonly);
        set_inconsistent                (multiImage && !pedited->sharpening.enabled);
    }

    setEnabled (pp->sharpening.enabled);

    eonlyConn.block (true);
    edgesonly->set_active (pp->sharpening.edgesonly);
    eonlyConn.block (false);
    lastEdgesOnly = pp->sharpening.edgesonly;

    hcConn.block (true);
    halocontrol->set_active (pp->sharpening.halocontrol);
    hcConn.block (false);
    lastHaloControl = pp->sharpening.halocontrol;

    amount->setValue        (pp->sharpening.amount);
    radius->setValue        (pp->sharpening.radius);
    threshold->setValue<int>(pp->sharpening.threshold);
    eradius->setValue       (pp->sharpening.edges_radius);
    etolerance->setValue    (pp->sharpening.edges_tolerance);
    hcamount->setValue      (pp->sharpening.halocontrol_amount);

    dradius->setValue       (pp->sharpening.deconvradius);
    damount->setValue       (pp->sharpening.deconvamount);
    diter->setValue         (pp->sharpening.deconviter);
    ddamping->setValue      (pp->sharpening.deconvdamping);

    if (!batchMode) {
        removeIfThere (edgebin, edgebox, false);

        if (edgesonly->get_active ()) {
            edgebin->pack_start (*edgebox);
        }

        removeIfThere (hcbin, hcbox, false);

        if (halocontrol->get_active ()) {
            hcbin->pack_start (*hcbox);
        }

    }

    if (pedited && !pedited->sharpening.method) {
        method->set_active (2);
    } else if (pp->sharpening.method == "usm") {
        method->set_active (0);
    } else if (pp->sharpening.method == "rld") {
        method->set_active (1);
    }

    enableListener ();
}

void Sharpening::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->sharpening.amount           = (int)amount->getValue();
    pp->sharpening.enabled          = getEnabled ();
    pp->sharpening.radius           = radius->getValue ();
    pp->sharpening.threshold        = threshold->getValue<int> ();
    pp->sharpening.edgesonly        = edgesonly->get_active ();
    pp->sharpening.edges_radius     = eradius->getValue ();
    pp->sharpening.edges_tolerance  = (int)etolerance->getValue ();
    pp->sharpening.halocontrol      = halocontrol->get_active ();
    pp->sharpening.halocontrol_amount = (int)hcamount->getValue ();
    pp->sharpening.deconvradius     = dradius->getValue ();
    pp->sharpening.deconviter       = (int)diter->getValue ();
    pp->sharpening.deconvamount     = (int)damount->getValue ();
    pp->sharpening.deconvdamping    = (int)ddamping->getValue ();

    if (method->get_active_row_number() == 0) {
        pp->sharpening.method = "usm";
    } else if (method->get_active_row_number() == 1) {
        pp->sharpening.method = "rld";
    }

    if (pedited) {
        pedited->sharpening.amount          = amount->getEditedState ();
        pedited->sharpening.radius          = radius->getEditedState ();
        pedited->sharpening.threshold       = threshold->getEditedState ();
        pedited->sharpening.edges_radius    = eradius->getEditedState ();
        pedited->sharpening.edges_tolerance = etolerance->getEditedState ();
        pedited->sharpening.halocontrol_amount = hcamount->getEditedState ();
        pedited->sharpening.deconvamount    = damount->getEditedState ();
        pedited->sharpening.deconvradius    = dradius->getEditedState ();
        pedited->sharpening.deconviter      = diter->getEditedState ();
        pedited->sharpening.deconvdamping   = ddamping->getEditedState ();
        pedited->sharpening.method          =  method->get_active_row_number() != 2;
        pedited->sharpening.halocontrol     =  !halocontrol->get_inconsistent();
        pedited->sharpening.edgesonly       =  !edgesonly->get_inconsistent();
        pedited->sharpening.enabled         =  !get_inconsistent();
    }
}

void Sharpening::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    amount->setDefault (defParams->sharpening.amount);
    radius->setDefault (defParams->sharpening.radius);
    threshold->setDefault<int> (defParams->sharpening.threshold);
    eradius->setDefault (defParams->sharpening.edges_radius);
    etolerance->setDefault (defParams->sharpening.edges_tolerance);
    hcamount->setDefault (defParams->sharpening.halocontrol_amount);
    damount->setDefault (defParams->sharpening.deconvamount);
    dradius->setDefault (defParams->sharpening.deconvradius);
    diter->setDefault (defParams->sharpening.deconviter);
    ddamping->setDefault (defParams->sharpening.deconvdamping);

    if (pedited) {
        amount->setDefaultEditedState       (pedited->sharpening.amount ? Edited : UnEdited);
        radius->setDefaultEditedState       (pedited->sharpening.radius ? Edited : UnEdited);
        threshold->setDefaultEditedState    (pedited->sharpening.threshold ? Edited : UnEdited);
        eradius->setDefaultEditedState      (pedited->sharpening.edges_radius ? Edited : UnEdited);
        etolerance->setDefaultEditedState   (pedited->sharpening.edges_tolerance ? Edited : UnEdited);
        hcamount->setDefaultEditedState     (pedited->sharpening.halocontrol_amount ? Edited : UnEdited);
        damount->setDefaultEditedState      (pedited->sharpening.deconvamount ? Edited : UnEdited);
        dradius->setDefaultEditedState      (pedited->sharpening.deconvradius ? Edited : UnEdited);
        diter->setDefaultEditedState        (pedited->sharpening.deconviter ? Edited : UnEdited);
        ddamping->setDefaultEditedState     (pedited->sharpening.deconvdamping ? Edited : UnEdited);
    } else {
        amount->setDefaultEditedState       (Irrelevant);
        radius->setDefaultEditedState       (Irrelevant);
        threshold->setDefaultEditedState    (Irrelevant);
        eradius->setDefaultEditedState      (Irrelevant);
        etolerance->setDefaultEditedState   (Irrelevant);
        hcamount->setDefaultEditedState     (Irrelevant);
        damount->setDefaultEditedState      (Irrelevant);
        dradius->setDefaultEditedState      (Irrelevant);
        diter->setDefaultEditedState        (Irrelevant);
        ddamping->setDefaultEditedState     (Irrelevant);
    }
}

void Sharpening::adjusterChanged (Adjuster* a, double newval)
{

    if (listener && (multiImage || getEnabled()) ) {

        Glib::ustring costr;

        if (a == radius || a == dradius) {
            costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());
        } else if (a == eradius) {
            costr = Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue());
        } else {
            costr = Glib::ustring::format ((int)a->getValue());
        }

        if (a == amount) {
            listener->panelChanged (EvShrAmount, costr);
        } else if (a == radius) {
            listener->panelChanged (EvShrRadius, costr);
        } else if (a == eradius) {
            listener->panelChanged (EvShrEdgeRadius, costr);
        } else if (a == etolerance) {
            listener->panelChanged (EvShrEdgeTolerance, costr);
        } else if (a == hcamount) {
            listener->panelChanged (EvShrHaloAmount, costr);
        } else if (a == dradius) {
            listener->panelChanged (EvShrDRadius, costr);
        } else if (a == damount) {
            listener->panelChanged (EvShrDAmount, costr);
        } else if (a == ddamping) {
            listener->panelChanged (EvShrDDamping, costr);
        } else if (a == diter) {
            listener->panelChanged (EvShrDIterations, costr);
        }
    }
}

void Sharpening::adjusterChanged (ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight)
{
    if (listener && (multiImage || getEnabled()) ) {
        if(a == threshold) {
            listener->panelChanged (EvShrThresh, threshold->getHistoryString());
        }
    }
}

void Sharpening::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvShrEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvShrEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvShrEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Sharpening::edgesonly_toggled ()
{

    if (multiImage) {
        if (edgesonly->get_inconsistent()) {
            edgesonly->set_inconsistent (false);
            eonlyConn.block (true);
            edgesonly->set_active (false);
            eonlyConn.block (false);
        } else if (lastEdgesOnly) {
            edgesonly->set_inconsistent (true);
        }

        lastEdgesOnly = edgesonly->get_active ();
    }

    if (!batchMode) {
        removeIfThere (edgebin, edgebox, false);

        if (edgesonly->get_active ()) {
            edgebin->pack_start (*edgebox);
        }
    }

    if (listener && (multiImage || getEnabled()) ) {
        if (edgesonly->get_inconsistent()) {
            listener->panelChanged (EvShrEdgeOnly, M("GENERAL_INITIALVALUES"));
        } else if (edgesonly->get_active ()) {
            listener->panelChanged (EvShrEdgeOnly, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvShrEdgeOnly, M("GENERAL_DISABLED"));
        }
    }
}

void Sharpening::halocontrol_toggled ()
{

    if (multiImage) {
        if (halocontrol->get_inconsistent()) {
            halocontrol->set_inconsistent (false);
            hcConn.block (true);
            halocontrol->set_active (false);
            hcConn.block (false);
        } else if (lastHaloControl) {
            halocontrol->set_inconsistent (true);
        }

        lastHaloControl = halocontrol->get_active ();
    }

    if (!batchMode) {
        removeIfThere (hcbin, hcbox, false);

        if (halocontrol->get_active ()) {
            hcbin->pack_start (*hcbox);
        }
    }

    if (listener && (multiImage || getEnabled()) ) {
        if (halocontrol->get_inconsistent()) {
            listener->panelChanged (EvShrHaloControl, M("GENERAL_INITIALVALUES"));
        } else if (halocontrol->get_active ()) {
            listener->panelChanged (EvShrHaloControl, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvShrHaloControl, M("GENERAL_DISABLED"));
        }
    }
}

void Sharpening::method_changed ()
{

    removeIfThere (this, usm, false);
    removeIfThere (this, rld, false);

    if (method->get_active_row_number() == 0) {
        pack_start (*usm);
    } else if (method->get_active_row_number() == 1) {
        pack_start (*rld);
    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvShrMethod, method->get_active_text ());
    }

}

void Sharpening::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);

    removeIfThere (hcbin, hcbox, false);
    hcbin->pack_start (*hcbox);
    removeIfThere (edgebin, edgebox, false);
    edgebin->pack_start (*edgebox);

    radius->showEditedCB ();
    amount->showEditedCB ();
    threshold->showEditedCB ();
    eradius->showEditedCB ();
    etolerance->showEditedCB ();
    hcamount->showEditedCB ();
    dradius->showEditedCB ();
    damount->showEditedCB ();
    ddamping->showEditedCB ();
    diter->showEditedCB ();
    method->append_text (M("GENERAL_UNCHANGED"));
}

void Sharpening::setAdjusterBehavior (bool amountadd)
{

    amount->setAddMode(amountadd);
    damount->setAddMode(amountadd);
}

void Sharpening::trimValues (rtengine::procparams::ProcParams* pp)
{

    amount->trimValue(pp->sharpening.amount);
    damount->trimValue(pp->sharpening.deconvamount);
}
