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
#include <cmath>
#include "eventmapper.h"
#include "prsharpening.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring PrSharpening::TOOL_NAME = "prsharpening";

PrSharpening::PrSharpening () : FoldableToolPanel(this, TOOL_NAME, M("TP_PRSHARPENING_LABEL"), false, true)
{

    auto m = ProcEventMapper::getInstance();
    EvPrShrContrast = m->newEvent(RESIZE, "HISTORY_MSG_PRSHARPEN_CONTRAST");

    std::vector<GradientMilestone> milestones;
    milestones.push_back( GradientMilestone(0.0, 0.0, 0.0, 0.0) );
    milestones.push_back( GradientMilestone(1.0, 1.0, 1.0, 1.0) );

    setEnabledTooltipMarkup(M("TP_PRSHARPENING_TOOLTIP"));

    Gtk::Box* hb = Gtk::manage (new Gtk::Box ());
    hb->show ();
    contrast = Gtk::manage(new Adjuster (M("TP_SHARPENING_CONTRAST"), 0, 200, 1, 15));
    contrast->setAdjusterListener (this);
    pack_start(*contrast);
    contrast->show();

    Gtk::Label* ml = Gtk::manage (new Gtk::Label (M("TP_SHARPENING_METHOD") + ":"));
    ml->show ();
    method = Gtk::manage (new MyComboBoxText ());
    method->append (M("TP_SHARPENING_USM"));
    method->append (M("TP_SHARPENING_RLD"));
    method->show ();
    hb->pack_start(*ml, Gtk::PACK_SHRINK, 4);
    hb->pack_start(*method);
    pack_start (*hb);

    rld = new Gtk::Box(Gtk::ORIENTATION_VERTICAL);
    dradius = Gtk::manage (new Adjuster (M("TP_SHARPENING_EDRADIUS"), 0.4, 2.5, 0.01, 0.45));
    damount = Gtk::manage (new Adjuster (M("TP_SHARPENING_RLD_AMOUNT"), 0.0, 100, 1, 100));
    ddamping = Gtk::manage (new Adjuster (M("TP_SHARPENING_RLD_DAMPING"), 0, 100, 1, 0));
    diter = Gtk::manage (new Adjuster (M("TP_SHARPENING_RLD_ITERATIONS"), 5, 100, 1, 100));
    rld->pack_start (*dradius);
    rld->pack_start (*damount);
    rld->pack_start (*ddamping);
    rld->pack_start (*diter);
    dradius->show ();
    damount->show ();
    ddamping->show ();
    diter->show ();
    rld->show ();

    usm = new Gtk::Box(Gtk::ORIENTATION_VERTICAL);
    usm->show ();

    Gtk::Separator* hsep6a = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    amount = Gtk::manage (new Adjuster (M("TP_SHARPENING_AMOUNT"), 1, 1000, 1, 200));
    radius = Gtk::manage (new Adjuster (M("TP_SHARPENING_RADIUS"), 0.3, 3, 0.01, 0.5));
    threshold = Gtk::manage (new ThresholdAdjuster (M("TP_SHARPENING_THRESHOLD"), 0., 2000., 20., 80., 2000., 1200., 0, false));
    threshold->setAdjusterListener (this);
    threshold->setBgGradient(milestones);
    pack_start(*hsep6a, Gtk::PACK_SHRINK, 2);

    pack_start (*usm);

    usm->pack_start(*radius);
    usm->pack_start(*amount);
    usm->pack_start(*threshold);
    hsep6a->show ();
    radius->show ();
    amount->show ();
    threshold->show ();

    Gtk::Separator* hsep6 = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    edgesonly = Gtk::manage (new Gtk::CheckButton (M("TP_SHARPENING_ONLYEDGES")));
    edgesonly->set_active (false);
    edgebox = new Gtk::Box(Gtk::ORIENTATION_VERTICAL);
    eradius = Gtk::manage (new Adjuster (M("TP_SHARPENING_EDRADIUS"), 0.5, 2.5, 0.1, 1.9));
    etolerance = Gtk::manage (new Adjuster (M("TP_SHARPENING_EDTOLERANCE"), 10, 10000, 100, 1800));
    usm->pack_start(*hsep6, Gtk::PACK_SHRINK, 2);
    usm->pack_start(*edgesonly);
    edgebox->pack_start(*eradius);
    edgebox->pack_start(*etolerance);
    edgebox->show ();
    edgebin = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    usm->pack_start (*edgebin);
    edgebin->show ();
    hsep6->show();
    edgesonly->show();
    eradius->show();
    etolerance->show();

    Gtk::Separator* hsep6b = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    halocontrol = Gtk::manage (new Gtk::CheckButton (M("TP_SHARPENING_HALOCONTROL")));
    halocontrol->set_active (false);
    hcbox = new Gtk::Box(Gtk::ORIENTATION_VERTICAL);
    hcamount = Gtk::manage (new Adjuster (M("TP_SHARPENING_HCAMOUNT"), 1, 100, 1, 75));
    usm->pack_start(*hsep6b, Gtk::PACK_SHRINK, 2);
    usm->pack_start(*halocontrol);
    hcbox->pack_start(*hcamount);
    hcbox->show ();
    hcbin = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
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

    eonlyConn = edgesonly->signal_toggled().connect( sigc::mem_fun(*this, &PrSharpening::edgesonly_toggled) );
    hcConn    = halocontrol->signal_toggled().connect( sigc::mem_fun(*this, &PrSharpening::halocontrol_toggled) );
    method->signal_changed().connect( sigc::mem_fun(*this, &PrSharpening::method_changed) );
}

PrSharpening::~PrSharpening ()
{

    delete usm;
    delete rld;
    delete edgebox;
    delete hcbox;
}


void PrSharpening::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        contrast->setEditedState    (pedited->prsharpening.contrast ? Edited : UnEdited);
        amount->setEditedState      (pedited->prsharpening.amount ? Edited : UnEdited);
        radius->setEditedState      (pedited->prsharpening.radius ? Edited : UnEdited);
        threshold->setEditedState   (pedited->prsharpening.threshold ? Edited : UnEdited);
        eradius->setEditedState     (pedited->prsharpening.edges_radius ? Edited : UnEdited);
        etolerance->setEditedState  (pedited->prsharpening.edges_tolerance ? Edited : UnEdited);
        hcamount->setEditedState    (pedited->prsharpening.halocontrol_amount ? Edited : UnEdited);
        damount->setEditedState     (pedited->prsharpening.deconvamount ? Edited : UnEdited);
        dradius->setEditedState     (pedited->prsharpening.deconvradius ? Edited : UnEdited);
        diter->setEditedState       (pedited->prsharpening.deconviter ? Edited : UnEdited);
        ddamping->setEditedState    (pedited->prsharpening.deconvdamping ? Edited : UnEdited);

        halocontrol->set_inconsistent   (multiImage && !pedited->prsharpening.halocontrol);
        edgesonly->set_inconsistent     (multiImage && !pedited->prsharpening.edgesonly);
        set_inconsistent                (multiImage && !pedited->prsharpening.enabled);
    }

    setEnabled (pp->prsharpening.enabled);

    eonlyConn.block (true);
    edgesonly->set_active (pp->prsharpening.edgesonly);
    eonlyConn.block (false);
    lastEdgesOnly = pp->prsharpening.edgesonly;

    hcConn.block (true);
    halocontrol->set_active (pp->prsharpening.halocontrol);
    hcConn.block (false);
    lastHaloControl = pp->prsharpening.halocontrol;

    contrast->setValue      (pp->prsharpening.contrast);
    amount->setValue        (pp->prsharpening.amount);
    radius->setValue        (pp->prsharpening.radius);
    threshold->setValue<int>(pp->prsharpening.threshold);
    eradius->setValue       (pp->prsharpening.edges_radius);
    etolerance->setValue    (pp->prsharpening.edges_tolerance);
    hcamount->setValue      (pp->prsharpening.halocontrol_amount);

    dradius->setValue       (pp->prsharpening.deconvradius);
    damount->setValue       (pp->prsharpening.deconvamount);
    diter->setValue         (pp->prsharpening.deconviter);
    ddamping->setValue      (pp->prsharpening.deconvdamping);

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

    if (pedited && !pedited->prsharpening.method) {
        method->set_active (2);
    } else if (pp->prsharpening.method == "usm") {
        method->set_active (0);
    } else if (pp->prsharpening.method == "rld") {
        method->set_active (1);
    }

    enableListener ();
}

void PrSharpening::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->prsharpening.contrast         = contrast->getValue ();
    pp->prsharpening.amount           = (int)amount->getValue();
    pp->prsharpening.enabled          = getEnabled ();
    pp->prsharpening.radius           = radius->getValue ();
    pp->prsharpening.threshold        = threshold->getValue<int> ();
    pp->prsharpening.edgesonly        = edgesonly->get_active ();
    pp->prsharpening.edges_radius     = eradius->getValue ();
    pp->prsharpening.edges_tolerance  = (int)etolerance->getValue ();
    pp->prsharpening.halocontrol      = halocontrol->get_active ();
    pp->prsharpening.halocontrol_amount = (int)hcamount->getValue ();
    pp->prsharpening.deconvradius     = dradius->getValue ();
    pp->prsharpening.deconviter       = (int)diter->getValue ();
    pp->prsharpening.deconvamount     = (int)damount->getValue ();
    pp->prsharpening.deconvdamping    = (int)ddamping->getValue ();

    if (method->get_active_row_number() == 0) {
        pp->prsharpening.method = "usm";
    } else if (method->get_active_row_number() == 1) {
        pp->prsharpening.method = "rld";
    }

    if (pedited) {
        pedited->prsharpening.contrast           = contrast->getEditedState ();
        pedited->prsharpening.amount             = amount->getEditedState ();
        pedited->prsharpening.radius             = radius->getEditedState ();
        pedited->prsharpening.threshold          = threshold->getEditedState ();
        pedited->prsharpening.edges_radius       = eradius->getEditedState ();
        pedited->prsharpening.edges_tolerance    = etolerance->getEditedState ();
        pedited->prsharpening.halocontrol_amount = hcamount->getEditedState ();
        pedited->prsharpening.deconvamount       = damount->getEditedState ();
        pedited->prsharpening.deconvradius       = dradius->getEditedState ();
        pedited->prsharpening.deconviter         = diter->getEditedState ();
        pedited->prsharpening.deconvdamping      = ddamping->getEditedState ();
        pedited->prsharpening.method             = method->get_active_row_number() != 2;
        pedited->prsharpening.halocontrol        = !halocontrol->get_inconsistent();
        pedited->prsharpening.edgesonly          = !edgesonly->get_inconsistent();
        pedited->prsharpening.enabled            = !get_inconsistent();
    }
}

void PrSharpening::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    contrast->setDefault (defParams->prsharpening.contrast);
    amount->setDefault (defParams->prsharpening.amount);
    radius->setDefault (defParams->prsharpening.radius);
    threshold->setDefault<int> (defParams->prsharpening.threshold);
    eradius->setDefault (defParams->prsharpening.edges_radius);
    etolerance->setDefault (defParams->prsharpening.edges_tolerance);
    hcamount->setDefault (defParams->prsharpening.halocontrol_amount);
    damount->setDefault (defParams->prsharpening.deconvamount);
    dradius->setDefault (defParams->prsharpening.deconvradius);
    diter->setDefault (defParams->prsharpening.deconviter);
    ddamping->setDefault (defParams->prsharpening.deconvdamping);

    if (pedited) {
        contrast->setDefaultEditedState     (pedited->prsharpening.contrast ? Edited : UnEdited);
        amount->setDefaultEditedState       (pedited->prsharpening.amount ? Edited : UnEdited);
        radius->setDefaultEditedState       (pedited->prsharpening.radius ? Edited : UnEdited);
        threshold->setDefaultEditedState    (pedited->prsharpening.threshold ? Edited : UnEdited);
        eradius->setDefaultEditedState      (pedited->prsharpening.edges_radius ? Edited : UnEdited);
        etolerance->setDefaultEditedState   (pedited->prsharpening.edges_tolerance ? Edited : UnEdited);
        hcamount->setDefaultEditedState     (pedited->prsharpening.halocontrol_amount ? Edited : UnEdited);
        damount->setDefaultEditedState      (pedited->prsharpening.deconvamount ? Edited : UnEdited);
        dradius->setDefaultEditedState      (pedited->prsharpening.deconvradius ? Edited : UnEdited);
        diter->setDefaultEditedState        (pedited->prsharpening.deconviter ? Edited : UnEdited);
        ddamping->setDefaultEditedState     (pedited->prsharpening.deconvdamping ? Edited : UnEdited);
    } else {
        contrast->setDefaultEditedState     (Irrelevant);
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

void PrSharpening::adjusterChanged (Adjuster* a, double newval)
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

        if (a == contrast) {
            listener->panelChanged (EvPrShrContrast, costr);
        } else if (a == amount) {
            listener->panelChanged (EvPrShrAmount, costr);
        } else if (a == radius) {
            listener->panelChanged (EvPrShrRadius, costr);
        } else if (a == eradius) {
            listener->panelChanged (EvPrShrEdgeRadius, costr);
        } else if (a == etolerance) {
            listener->panelChanged (EvPrShrEdgeTolerance, costr);
        } else if (a == hcamount) {
            listener->panelChanged (EvPrShrHaloAmount, costr);
        } else if (a == dradius) {
            listener->panelChanged (EvPrShrDRadius, costr);
        } else if (a == damount) {
            listener->panelChanged (EvPrShrDAmount, costr);
        } else if (a == ddamping) {
            listener->panelChanged (EvPrShrDDamping, costr);
        } else if (a == diter) {
            listener->panelChanged (EvPrShrDIterations, costr);
        }
    }
}

void PrSharpening::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvPrShrEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvPrShrEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvPrShrEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void PrSharpening::edgesonly_toggled ()
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
            listener->panelChanged (EvPrShrEdgeOnly, M("GENERAL_INITIALVALUES"));
        } else if (edgesonly->get_active ()) {
            listener->panelChanged (EvPrShrEdgeOnly, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvPrShrEdgeOnly, M("GENERAL_DISABLED"));
        }
    }
}

void PrSharpening::halocontrol_toggled ()
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
            listener->panelChanged (EvPrShrHaloControl, M("GENERAL_INITIALVALUES"));
        } else if (halocontrol->get_active ()) {
            listener->panelChanged (EvPrShrHaloControl, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvPrShrHaloControl, M("GENERAL_DISABLED"));
        }
    }
}

void PrSharpening::method_changed ()
{

    if (!batchMode) {
        removeIfThere (this, usm, false);
        removeIfThere (this, rld, false);

        if (method->get_active_row_number() == 0) {
            pack_start (*usm);
        } else if (method->get_active_row_number() == 1) {
            pack_start (*rld);
        }
    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvPrShrMethod, method->get_active_text ());
    }

}

void PrSharpening::adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop)
{
}

void PrSharpening::adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight)
{
}

void PrSharpening::adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop)
{
}

void PrSharpening::adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight)
{
    if (listener && (multiImage || getEnabled()) ) {
        if(a == threshold) {
            listener->panelChanged (EvPrShrThresh, threshold->getHistoryString());
        }
    }
}

void PrSharpening::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
}

void PrSharpening::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);

    removeIfThere (hcbin, hcbox, false);
    hcbin->pack_start (*hcbox);
    removeIfThere (edgebin, edgebox, false);
    edgebin->pack_start (*edgebox);
    pack_start (*rld);

    contrast->showEditedCB ();
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
    method->append (M("GENERAL_UNCHANGED"));
}

void PrSharpening::setAdjusterBehavior (bool contrastadd, bool radiusadd, bool amountadd, bool dampingadd, bool iteradd, bool edgetoladd, bool haloctrladd)
{

    contrast->setAddMode(contrastadd);
    radius->setAddMode(radiusadd);
    dradius->setAddMode(radiusadd);
    amount->setAddMode(amountadd);
    damount->setAddMode(amountadd);
    ddamping->setAddMode(dampingadd);
    diter->setAddMode(iteradd);
    eradius->setAddMode(radiusadd);
    etolerance->setAddMode(edgetoladd);
    hcamount->setAddMode(haloctrladd);
}

void PrSharpening::trimValues (rtengine::procparams::ProcParams* pp)
{

    contrast->trimValue(pp->prsharpening.contrast);
    radius->trimValue(pp->prsharpening.radius);
    dradius->trimValue(pp->prsharpening.deconvradius);
    amount->trimValue(pp->prsharpening.amount);
    damount->trimValue(pp->prsharpening.deconvamount);
    ddamping->trimValue(pp->prsharpening.deconvdamping);
    diter->trimValue(pp->prsharpening.deconviter);
    eradius->trimValue(pp->prsharpening.edges_radius);
    etolerance->trimValue(pp->prsharpening.edges_tolerance);
    hcamount->trimValue(pp->prsharpening.halocontrol_amount);
}
