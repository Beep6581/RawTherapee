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
#include <iomanip>

#include "labcurve.h"
#include "eventmapper.h"

#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "options.h"

#include "../rtengine/color.h"
#include "../rtengine/improcfun.h"
#include "../rtengine/procparams.h"
#include "editcallbacks.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring LCurve::TOOL_NAME = "labcurves";

LCurve::LCurve() : FoldableToolPanel(this, TOOL_NAME, M("TP_LABCURVE_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    Evgamutmunsell = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_GAMUTMUNSEL");
    CurveListener::setMulti(true);
    brightness = Gtk::manage(new Adjuster(M("TP_LABCURVE_BRIGHTNESS"), -100., 100., 1., 0.));
    contrast   = Gtk::manage(new Adjuster(M("TP_LABCURVE_CONTRAST"), -100., 100., 1., 0.));
    chromaticity   = Gtk::manage(new Adjuster(M("TP_LABCURVE_CHROMATICITY"), -100., 100., 1., 0.));
    chromaticity->set_tooltip_markup(M("TP_LABCURVE_CHROMA_TOOLTIP"));

    pack_start(*brightness);
    brightness->show();

    pack_start(*contrast);
    contrast->show();

    pack_start(*chromaticity);
    chromaticity->show();

    brightness->setAdjusterListener(this);
    contrast->setAdjusterListener(this);
    chromaticity->setAdjusterListener(this);

    brightness->setLogScale(2, 0, true);
    contrast->setLogScale(2, 0, true);
    chromaticity->setLogScale(2, 0, true);

    //%%%%%%%%%%%%%%%%%%
    Gtk::Separator* hsep2 = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    hsep2->show();
    pack_start(*hsep2, Gtk::PACK_EXPAND_WIDGET, 4);




    Gtk::Box* metHBox = Gtk::manage(new Gtk::Box());
    metHBox->set_spacing(2);
    Gtk::Label* metLabel = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_AVOID") + ":"));
    metHBox->pack_start(*metLabel, Gtk::PACK_SHRINK);

    gamutmunselmethod =  Gtk::manage(new MyComboBoxText());
    gamutmunselmethod->append(M("TP_LOCALLAB_GAMUTNON"));
    gamutmunselmethod->append(M("TP_LOCALLAB_GAMUTLABRELA"));
    gamutmunselmethod->append(M("TP_LOCALLAB_GAMUTXYZABSO"));
    gamutmunselmethod->append(M("TP_LOCALLAB_GAMUTXYZRELA"));
    gamutmunselmethod->append(M("TP_LOCALLAB_GAMUTMUNSELL"));
    gamutmunselmethod->set_active(4);
    gamutmunselmethod->set_tooltip_text(M("TP_LOCALLAB_AVOIDCOLORSHIFT_TOOLTIP"));
    metHBox->pack_start(*gamutmunselmethod);
    pack_start(*metHBox);
    gamutmunselmethodconn = gamutmunselmethod->signal_changed().connect(sigc::mem_fun(*this, &LCurve::gamutmunselChanged));


    lcredsk = Gtk::manage(new Gtk::CheckButton(M("TP_LABCURVE_LCREDSK")));
    lcredsk->set_tooltip_markup(M("TP_LABCURVE_LCREDSK_TOOLTIP"));
    pack_start(*lcredsk);

    rstprotection = Gtk::manage(new Adjuster(M("TP_LABCURVE_RSTPROTECTION"), 0., 100., 0.1, 0.));
    pack_start(*rstprotection);
    rstprotection->show();

    rstprotection->setAdjusterListener(this);
    rstprotection->set_tooltip_text(M("TP_LABCURVE_RSTPRO_TOOLTIP"));

    lcconn = lcredsk->signal_toggled().connect(sigc::mem_fun(*this, &LCurve::lcredsk_toggled));

    //%%%%%%%%%%%%%%%%%%%

    Gtk::Separator* hsep3 = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    hsep3->show();
    pack_start(*hsep3, Gtk::PACK_EXPAND_WIDGET, 4);

    curveEditorG = new CurveEditorGroup(options.lastLabCurvesDir);
    curveEditorG->setCurveListener(this);

    lshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "L*"));
    lshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_LL_TOOLTIP"));
    lshape->setEditID(EUID_Lab_LCurve, BT_SINGLEPLANE_FLOAT);

    ashape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "a*"));
    ashape->setEditID(EUID_Lab_aCurve, BT_SINGLEPLANE_FLOAT);

    ashape->setRangeLabels(
        M("TP_LABCURVE_CURVEEDITOR_A_RANGE1"), M("TP_LABCURVE_CURVEEDITOR_A_RANGE2"),
        M("TP_LABCURVE_CURVEEDITOR_A_RANGE3"), M("TP_LABCURVE_CURVEEDITOR_A_RANGE4")
    );
    //from green to magenta
    std::vector<GradientMilestone> milestones = {
        GradientMilestone(0., 0., 1., 0.),
        GradientMilestone(1., 1., 0., 1.)
    };
    ashape->setBottomBarBgGradient(milestones);
    ashape->setLeftBarBgGradient(milestones);

    bshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "b*"));
    bshape->setRangeLabels(
        M("TP_LABCURVE_CURVEEDITOR_B_RANGE1"), M("TP_LABCURVE_CURVEEDITOR_B_RANGE2"),
        M("TP_LABCURVE_CURVEEDITOR_B_RANGE3"), M("TP_LABCURVE_CURVEEDITOR_B_RANGE4")
    );
    bshape->setEditID(EUID_Lab_bCurve, BT_SINGLEPLANE_FLOAT);

    //from blue to yellow
    milestones = {
        GradientMilestone(0., 0., 0., 1.),
        GradientMilestone(1., 1., 1., 0.)
    };
    bshape->setBottomBarBgGradient(milestones);
    bshape->setLeftBarBgGradient(milestones);

    curveEditorG->newLine();  //  ------------------------------------------------ second line

    lhshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_LABCURVE_CURVEEDITOR_LH")));
    lhshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_LH_TOOLTIP"));
    lhshape->setCurveColorProvider(this, 4);
    lhshape->setEditID(EUID_Lab_LHCurve, BT_SINGLEPLANE_FLOAT);


    chshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_LABCURVE_CURVEEDITOR_CH")));
    chshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_CH_TOOLTIP"));
    chshape->setCurveColorProvider(this, 1);
    chshape->setEditID(EUID_Lab_CHCurve, BT_SINGLEPLANE_FLOAT);


    hhshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_LABCURVE_CURVEEDITOR_HH")));
    hhshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_HH_TOOLTIP"));
    hhshape->setCurveColorProvider(this, 5);
    hhshape->setEditID(EUID_Lab_HHCurve, BT_SINGLEPLANE_FLOAT);

    curveEditorG->newLine();  //  ------------------------------------------------ 3rd line

    ccshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_LABCURVE_CURVEEDITOR_CC")));
    ccshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_CC_TOOLTIP"));
    ccshape->setEditID(EUID_Lab_CCurve, BT_SINGLEPLANE_FLOAT);
    ccshape->setRangeLabels(
        M("TP_LABCURVE_CURVEEDITOR_CC_RANGE1"), M("TP_LABCURVE_CURVEEDITOR_CC_RANGE2"),
        M("TP_LABCURVE_CURVEEDITOR_CC_RANGE3"), M("TP_LABCURVE_CURVEEDITOR_CC_RANGE4")
    );

    ccshape->setBottomBarColorProvider(this, 2);
    ccshape->setLeftBarColorProvider(this, 7);
    ccshape->setRangeDefaultMilestones(0.05, 0.2, 0.58);

    lcshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_LABCURVE_CURVEEDITOR_LC")));
    lcshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_LC_TOOLTIP"));
    lcshape->setEditID(EUID_Lab_LCCurve, BT_SINGLEPLANE_FLOAT);

    // left and bottom bar uses the same caller id because the will display the same content
    lcshape->setBottomBarColorProvider(this, 2);
    lcshape->setRangeLabels(
        M("TP_LABCURVE_CURVEEDITOR_CC_RANGE1"), M("TP_LABCURVE_CURVEEDITOR_CC_RANGE2"),
        M("TP_LABCURVE_CURVEEDITOR_CC_RANGE3"), M("TP_LABCURVE_CURVEEDITOR_CC_RANGE4")
    );
    lcshape->setRangeDefaultMilestones(0.05, 0.2, 0.58);

    clshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_LABCURVE_CURVEEDITOR_CL")));
    clshape->setTooltip(M("TP_LABCURVE_CURVEEDITOR_CL_TOOLTIP"));
    clshape->setEditID(EUID_Lab_CLCurve, BT_SINGLEPLANE_FLOAT);

    clshape->setLeftBarColorProvider(this, 7);
    clshape->setRangeDefaultMilestones(0.25, 0.5, 0.75);

    milestones = {
        GradientMilestone(0., 0., 0., 0.),
        GradientMilestone(1., 1., 1., 1.)
    };
    clshape->setBottomBarBgGradient(milestones);


    // Setting the gradient milestones

    // from black to white
    milestones.emplace_back(0., 0., 0., 0.);
    milestones.emplace_back(1., 1., 1., 1.);
    lshape->setBottomBarBgGradient(milestones);
    lshape->setLeftBarBgGradient(milestones);
    milestones.emplace_back(0., 0., 0., 0.);
    milestones.emplace_back(1., 1., 1., 1.);
    lcshape->setRangeDefaultMilestones(0.05, 0.2, 0.58);

    lcshape->setBottomBarBgGradient(milestones);

    milestones.at(0).r = milestones.at(0).g = milestones.at(0).b = 0.1;
    milestones.at(1).r = milestones.at(1).g = milestones.at(1).b = 0.8;
    lcshape->setLeftBarBgGradient(milestones);

    // whole hue range
    milestones.clear();

    for (int i = 0; i < 7; i++) {
        float R, G, B;
        float x = i / 6.0;
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        milestones.emplace_back(x, R, G, B);
    }

    chshape->setBottomBarBgGradient(milestones);
    lhshape->setBottomBarBgGradient(milestones);
    hhshape->setBottomBarBgGradient(milestones);


    // This will add the reset button at the end of the curveType buttons
    curveEditorG->curveListComplete();

    pack_start(*curveEditorG, Gtk::PACK_SHRINK, 4);
    Gtk::Separator* hsepdh = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    hsepdh->show();
    pack_start(*hsepdh, Gtk::PACK_EXPAND_WIDGET, 4);
    show_all_children();

}

LCurve::~LCurve()
{
    delete curveEditorG;

}

void LCurve::read(const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener();
    gamutmunselmethodconn.block(true);


    brightness->setValue(pp->labCurve.brightness);
    contrast->setValue(pp->labCurve.contrast);
    chromaticity->setValue(pp->labCurve.chromaticity);
    adjusterChanged(chromaticity, pp->labCurve.chromaticity); // To update the GUI sensitiveness
    //%%%%%%%%%%%%%%%%%%%%%%
    rstprotection->setValue(pp->labCurve.rstprotection);


    bwtconn.block(true);
    lcconn.block(true);
    lcredsk->set_active(pp->labCurve.lcredsk);

    bwtconn.block(false);
    lcconn.block(false);

    lastLCVal = pp->labCurve.lcredsk;
    //%%%%%%%%%%%%%%%%%%%%%%

    lshape->setCurve(pp->labCurve.lcurve);
    ashape->setCurve(pp->labCurve.acurve);
    bshape->setCurve(pp->labCurve.bcurve);
    ccshape->setCurve(pp->labCurve.cccurve);
    chshape->setCurve(pp->labCurve.chcurve);
    lhshape->setCurve(pp->labCurve.lhcurve);
    hhshape->setCurve(pp->labCurve.hhcurve);
    lcshape->setCurve(pp->labCurve.lccurve);
    clshape->setCurve(pp->labCurve.clcurve);

    if (pedited && !pedited->labCurve.gamutmunselmethod) {
        gamutmunselmethod->set_active(4);     // "Unchanged"
    } else if (pp->labCurve.gamutmunselmethod == "NONE") {
        gamutmunselmethod->set_active(0);
    } else if (pp->labCurve.gamutmunselmethod == "LAB") {
        gamutmunselmethod->set_active(1);
    } else if (pp->labCurve.gamutmunselmethod == "XYZ") {
        gamutmunselmethod->set_active(2);
    } else if (pp->labCurve.gamutmunselmethod == "XYZREL") {
        gamutmunselmethod->set_active(3);
    } else if (pp->labCurve.gamutmunselmethod == "MUN") {
        gamutmunselmethod->set_active(4);
    }

    gamutmunselChanged();

    if (pedited) {
        brightness->setEditedState(pedited->labCurve.brightness ? Edited : UnEdited);
        contrast->setEditedState(pedited->labCurve.contrast ? Edited : UnEdited);
        chromaticity->setEditedState(pedited->labCurve.chromaticity ? Edited : UnEdited);

        //%%%%%%%%%%%%%%%%%%%%%%
        rstprotection->setEditedState(pedited->labCurve.rstprotection ? Edited : UnEdited);
        lcredsk->set_inconsistent(!pedited->labCurve.lcredsk);

        //%%%%%%%%%%%%%%%%%%%%%%

        lshape->setUnChanged(!pedited->labCurve.lcurve);
        ashape->setUnChanged(!pedited->labCurve.acurve);
        bshape->setUnChanged(!pedited->labCurve.bcurve);
        ccshape->setUnChanged(!pedited->labCurve.cccurve);
        chshape->setUnChanged(!pedited->labCurve.chcurve);
        lhshape->setUnChanged(!pedited->labCurve.lhcurve);
        hhshape->setUnChanged(!pedited->labCurve.hhcurve);
        lcshape->setUnChanged(!pedited->labCurve.lccurve);
        clshape->setUnChanged(!pedited->labCurve.clcurve);

        if (!pedited->labCurve.gamutmunselmethod) {
            gamutmunselmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        set_inconsistent(multiImage && !pedited->labCurve.enabled);
    }

    gamutmunselmethodconn.block(false);


    setEnabled(pp->labCurve.enabled);

    queue_draw();

    enableListener();
}


void LCurve::autoOpenCurve()
{
    // Open up the first curve if selected
    bool active = lshape->openIfNonlinear();

    if (!active) {
        active = ashape->openIfNonlinear();
    }

    if (!active) {
        active = bshape->openIfNonlinear();
    }

    if (!active) {
        active = ccshape->openIfNonlinear();
    }

    if (!active) {
        active = chshape->openIfNonlinear();
    }

    if (!active) {
        active = lhshape->openIfNonlinear();
    }

    if (!active) {
        active = hhshape->openIfNonlinear();
    }

    if (!active) {
        active = lcshape->openIfNonlinear();
    }

    if (!active) {
        clshape->openIfNonlinear();
    }

}

void LCurve::setEditProvider(EditDataProvider *provider)
{
    lshape->setEditProvider(provider);
    ccshape->setEditProvider(provider);
    lcshape->setEditProvider(provider);
    clshape->setEditProvider(provider);
    lhshape->setEditProvider(provider);
    chshape->setEditProvider(provider);
    hhshape->setEditProvider(provider);
    ashape->setEditProvider(provider);
    bshape->setEditProvider(provider);

}


void LCurve::write(ProcParams* pp, ParamsEdited* pedited)
{
    pp->labCurve.enabled = getEnabled();

    pp->labCurve.brightness    = brightness->getValue();
    pp->labCurve.contrast      = (int)contrast->getValue();
    pp->labCurve.chromaticity  = (int)chromaticity->getValue();
    //%%%%%%%%%%%%%%%%%%%%%%
    pp->labCurve.lcredsk         = lcredsk->get_active();

    pp->labCurve.rstprotection   = rstprotection->getValue();
    //%%%%%%%%%%%%%%%%%%%%%%

    pp->labCurve.lcurve  = lshape->getCurve();
    pp->labCurve.acurve  = ashape->getCurve();
    pp->labCurve.bcurve  = bshape->getCurve();
    pp->labCurve.cccurve = ccshape->getCurve();
    pp->labCurve.chcurve = chshape->getCurve();
    pp->labCurve.lhcurve = lhshape->getCurve();
    pp->labCurve.hhcurve = hhshape->getCurve();
    pp->labCurve.lccurve = lcshape->getCurve();
    pp->labCurve.clcurve = clshape->getCurve();



    if (pedited) {
        pedited->labCurve.brightness   = brightness->getEditedState();
        pedited->labCurve.contrast     = contrast->getEditedState();
        pedited->labCurve.chromaticity = chromaticity->getEditedState();

        //%%%%%%%%%%%%%%%%%%%%%%
        pedited->labCurve.lcredsk         = !lcredsk->get_inconsistent();

        pedited->labCurve.rstprotection   = rstprotection->getEditedState();
        pedited->labCurve.gamutmunselmethod = gamutmunselmethod->get_active_text() != M("GENERAL_UNCHANGED");

        pedited->labCurve.lcurve    = !lshape->isUnChanged();
        pedited->labCurve.acurve    = !ashape->isUnChanged();
        pedited->labCurve.bcurve    = !bshape->isUnChanged();
        pedited->labCurve.cccurve   = !ccshape->isUnChanged();
        pedited->labCurve.chcurve   = !chshape->isUnChanged();
        pedited->labCurve.lhcurve   = !lhshape->isUnChanged();
        pedited->labCurve.hhcurve   = !hhshape->isUnChanged();
        pedited->labCurve.lccurve   = !lcshape->isUnChanged();
        pedited->labCurve.clcurve   = !clshape->isUnChanged();

        pedited->labCurve.enabled = !get_inconsistent();

    }

    if (gamutmunselmethod->get_active_row_number() == 0) {
        pp->labCurve.gamutmunselmethod = "NONE";
    } else if (gamutmunselmethod->get_active_row_number() == 1) {
        pp->labCurve.gamutmunselmethod = "LAB";
    } else if (gamutmunselmethod->get_active_row_number() == 2) {
        pp->labCurve.gamutmunselmethod = "XYZ";
    } else if (gamutmunselmethod->get_active_row_number() == 3) {
        pp->labCurve.gamutmunselmethod = "XYZREL";
    } else if (gamutmunselmethod->get_active_row_number() == 4) {
        pp->labCurve.gamutmunselmethod = "MUN";
    }



}

void LCurve::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited)
{

    brightness->setDefault(defParams->labCurve.brightness);
    contrast->setDefault(defParams->labCurve.contrast);
    chromaticity->setDefault(defParams->labCurve.chromaticity);
    rstprotection->setDefault(defParams->labCurve.rstprotection);

    if (pedited) {
        brightness->setDefaultEditedState(pedited->labCurve.brightness ? Edited : UnEdited);
        contrast->setDefaultEditedState(pedited->labCurve.contrast ? Edited : UnEdited);
        chromaticity->setDefaultEditedState(pedited->labCurve.chromaticity ? Edited : UnEdited);
        rstprotection->setDefaultEditedState(pedited->labCurve.rstprotection ? Edited : UnEdited);
    } else {
        brightness->setDefaultEditedState(Irrelevant);
        contrast->setDefaultEditedState(Irrelevant);
        chromaticity->setDefaultEditedState(Irrelevant);
        rstprotection->setDefaultEditedState(Irrelevant);
    }
}

//%%%%%%%%%%%%%%%%%%%%%%

void LCurve::gamutmunselChanged()
{

    if (listener && (multiImage || getEnabled())) {
        listener->panelChanged(Evgamutmunsell, gamutmunselmethod->get_active_text());
    }

}


void LCurve::lcredsk_toggled()
{

    if (batchMode) {
        if (lcredsk->get_inconsistent()) {
            lcredsk->set_inconsistent(false);
            lcconn.block(true);
            lcredsk->set_active(false);
            lcconn.block(false);
        } else if (lastLCVal) {
            lcredsk->set_inconsistent(true);
        }

        lastLCVal = lcredsk->get_active();
    } else {
        lcshape->refresh();
    }

    if (listener && getEnabled()) {
        if (lcredsk->get_active()) {
            listener->panelChanged(EvLLCredsk, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvLLCredsk, M("GENERAL_DISABLED"));
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%

/*
 * Curve listener
 *
 * If more than one curve has been added, the curve listener is automatically
 * set to 'multi=true', and send a pointer of the modified curve in a parameter
 */
void LCurve::curveChanged(CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == lshape) {
            listener->panelChanged(EvLLCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == ashape) {
            listener->panelChanged(EvLaCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == bshape) {
            listener->panelChanged(EvLbCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == ccshape) {
            listener->panelChanged(EvLCCCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == chshape) {
            listener->panelChanged(EvLCHCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == lhshape) {
            listener->panelChanged(EvLLHCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == hhshape) {
            listener->panelChanged(EvLHHCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == lcshape) {
            listener->panelChanged(EvLLCCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == clshape) {
            listener->panelChanged(EvLCLCurve, M("HISTORY_CUSTOMCURVE"));
        }


    }
}

void LCurve::adjusterChanged(Adjuster* a, double newval)
{
    Glib::ustring costr;

    if (a == brightness) {
        costr = Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(2), a->getValue());
    } else if (a == rstprotection) {
        costr = Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(1), a->getValue());
    } else {
        costr = Glib::ustring::format((int)a->getValue());
    }

    if (a == brightness) {
        if (listener && getEnabled()) {
            listener->panelChanged(EvLBrightness, costr);
        }
    } else if (a == contrast) {
        if (listener && getEnabled()) {
            listener->panelChanged(EvLContrast, costr);
        }
    } else if (a == rstprotection) {
        if (listener && getEnabled()) {
            listener->panelChanged(EvLRSTProtection, costr);
        }
    } else if (a == chromaticity) {
        if (multiImage) {
            //if chromaticity==-100 (lowest value), we enter the B&W mode and avoid color shift and rstprotection has no effect
            rstprotection->set_sensitive(true);
            lcredsk->set_sensitive(true);
        } else {
            //if chromaticity==-100 (lowest value), we enter the B&W mode and avoid color shift and rstprotection has no effect
            rstprotection->set_sensitive(int(newval) > -100);   //no reason for grey rstprotection
            lcredsk->set_sensitive(int(newval) > -100);
        }

        if (listener && getEnabled()) {
            listener->panelChanged(EvLSaturation, costr);
        }
    }
}

void LCurve::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R = 0.f, G = 0.f, B = 0.f;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve
        Color::hsv2rgb01(float(valX), float(valY), 0.5f, R, G, B);
    } else if (callerId == 2) {  // cc - bottom bar
        float value = (1.f - 0.7f) * float(valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float(valY * 0.8), float(valX), value, R, G, B);
    } else if (callerId == 6) {  // cc - left bar
        float value = (1.f - 0.7f) * float(valX) + 0.7f;
        float hue = (1.14056f - 0.92f) * float(valY) + 0.92f;

        if (hue > 1.0f) {
            hue -= 1.0f;
        }

        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(hue, float(valX), value, R, G, B);
    } else if (callerId == 3) {  // lc - bottom bar
        float value = (1.f - 0.7f) * float(valX) + 0.7f;

        if (lcredsk->get_active()) {
            // skin range
            // -0.1 rad < Hue < 1.6 rad
            // Y axis / from 0.92 up to 0.14056
            float hue = (1.14056f - 0.92f) * float(valY) + 0.92f;

            if (hue > 1.0f) {
                hue -= 1.0f;
            }

            // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
            Color::hsv2rgb01(hue, float(valX), value, R, G, B);
        } else {
            // whole hue range
            // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
            Color::hsv2rgb01(float(valY * 0.8), float(valX), value, R, G, B);
        }
    } else if (callerId == 4) {  // LH - bottom bar
        Color::hsv2rgb01(float(valX), 0.5f, float(valY), R, G, B);
    } else if (callerId == 5) {  // HH - bottom bar
        float h = float((valY - 0.5) * 0.3 + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01(h, 0.5f, 0.5f, R, G, B);
    } else if (callerId == 7) {  // cc and cl - left bar
        float value = (1.f - 0.7f) * float(valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float(valY * 0.8), 1.f - float(valX), value, R, G, B);
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
}

void LCurve::setBatchMode(bool batchMode)
{

    ToolPanel::setBatchMode(batchMode);
    brightness->showEditedCB();
    contrast->showEditedCB();
    chromaticity->showEditedCB();
    rstprotection->showEditedCB();
    curveEditorG->setBatchMode(batchMode);
    lcshape->setBottomBarColorProvider(nullptr, -1);
    lcshape->setLeftBarColorProvider(nullptr, -1);
    gamutmunselmethod->append(M("GENERAL_UNCHANGED"));

}


void LCurve::updateCurveBackgroundHistogram(
    const LUTu& histToneCurve,
    const LUTu& histLCurve,
    const LUTu& histCCurve,
    const LUTu& histLCAM,
    const LUTu& histCCAM,
    const LUTu& histRed,
    const LUTu& histGreen,
    const LUTu& histBlue,
    const LUTu& histLuma,
    const LUTu& histLRETI
)
{
    lshape->updateBackgroundHistogram(histLCurve);
    ccshape->updateBackgroundHistogram(histCCurve);
    lcshape->updateBackgroundHistogram(histCCurve);
    clshape->updateBackgroundHistogram(histLCurve);
}

void LCurve::setAdjusterBehavior(bool bradd, bool contradd, bool satadd)
{

    brightness->setAddMode(bradd);
    contrast->setAddMode(contradd);
    chromaticity->setAddMode(satadd);
}

void LCurve::trimValues(rtengine::procparams::ProcParams* pp)
{

    brightness->trimValue(pp->labCurve.brightness);
    contrast->trimValue(pp->labCurve.contrast);
    chromaticity->trimValue(pp->labCurve.chromaticity);
}

void LCurve::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvLEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvLEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvLEnabled, M("GENERAL_DISABLED"));
        }
    }
}
