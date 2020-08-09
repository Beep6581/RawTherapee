/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>frame
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
 *  2019-2020 Pierre Cabrera <pierre.cab@gmail.com>
 */
#include "locallabtools.h"

#include "options.h"
#include "../rtengine/procparams.h"
#include "locallab.h"
#include "rtimage.h"

#define MINNEIGH 0.1
#define MAXNEIGH 1500
#define CENTERNEIGH 200

using namespace rtengine;
using namespace procparams;

extern Options options;

static double retiSlider2neigh(double sval)
{
    // Slider range: 0 - 5000
    double neigh;

    if (sval <= 200) {
        // Linear below center-temp
        neigh = MINNEIGH + (sval / 200.0) * (CENTERNEIGH - MINNEIGH);
    } else {
        const double slope = (double)(CENTERNEIGH - MINNEIGH) / (MAXNEIGH - CENTERNEIGH);
        const double x = (sval - 200) / 200; // x range: 0 - 1
        const double y = x * slope + (1.0 - slope) * pow(x, 4.0);
        neigh = CENTERNEIGH + y * (MAXNEIGH - CENTERNEIGH);
    }

    if (neigh < MINNEIGH) {
        neigh = MINNEIGH;
    }

    if (neigh > MAXNEIGH) {
        neigh = MAXNEIGH;
    }

    return neigh;
}

static double retiNeigh2Slider(double neigh)
{
    double sval;

    if (neigh <= CENTERNEIGH) {
        sval = ((neigh - MINNEIGH) / (CENTERNEIGH - MINNEIGH)) * 200.0;
    } else {
        const double slope = (double)(CENTERNEIGH - MINNEIGH) / (MAXNEIGH - CENTERNEIGH);
        const double y = (neigh - CENTERNEIGH) / (MAXNEIGH - CENTERNEIGH);
        double x = pow(y, 0.25); // Rough guess of x, will be a little lower
        double k = 0.1;
        bool add = true;

        // The y=f(x) function is a mess to invert, therefore we have this trial-refinement loop instead.
        // From tests, worst case is about 20 iterations, i.e. no problem
        for (;;) {
            double y1 = x * slope + (1.0 - slope) * pow(x, 4.0);

            if (200 * fabs(y1 - y) < 0.1) {
                break;
            }

            if (y1 < y) {
                if (!add) {
                    k /= 2;
                }

                x += k;
                add = true;
            } else {
                if (add) {
                    k /= 2;
                }

                x -= k;
                add = false;
            }
        }

        sval = 200.0 + x * 200.0;
    }

    if (sval < 0.) {
        sval = 0.;
    }

    if (sval > 1500.) {
        sval = 1500.;
    }

    return sval;
}

/* ==== LocallabTone ==== */
LocallabTone::LocallabTone():
    LocallabTool(this, M("TP_LOCALLAB_TONE_TOOLNAME"), M("TP_LOCALLAB_TM"), true),

    // Tone mapping specific widgets
    amount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_AMOUNT"), 50., 100.0, 0.5, 95.))),
    stren(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STREN"), -0.5, 2.0, 0.01, 0.5))),
    equiltm(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_EQUIL")))),
    gamma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAM"), 0.4, 4.0, 0.11, 1.0))),
    satur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SATUR"), -100., 100., 0.1, 0.))), // By default satur = 0 ==> use Mantiuk value
    estop(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ESTOP"), 0.1, 4., 0.01, 1.4))),
    scaltm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALTM"), 0.1, 10.0, 0.01, 1.0))),
    rewei(Gtk::manage(new Adjuster(M("TP_LOCALLAB_REWEI"), 0, 3, 1, 0))),
    softradiustm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    sensitm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
    expmasktm(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWT")))),
    showmasktmMethod(Gtk::manage(new MyComboBoxText())),
    enatmMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    enatmMaskaft(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_AFTER_MASK")))),
    masktmCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmasktmshape(static_cast<FlatCurveEditor*>(masktmCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmasktmshape(static_cast<FlatCurveEditor*>(masktmCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmasktmshape(static_cast<FlatCurveEditor *>(masktmCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    lapmasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    radmasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.05, 5.0, 0.01, 1.))),
    slomasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2tmCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmasktmshape(static_cast<DiagonalCurveEditor*>(mask2tmCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
{
    const LocallabParams::LocallabSpot defSpot;

    // Parameter Tone Mapping specific widgets
    amount->setAdjusterListener(this);

    stren->setAdjusterListener(this);

    equiltmConn = equiltm->signal_toggled().connect(sigc::mem_fun(*this, &LocallabTone::equiltmChanged));

    gamma->setAdjusterListener(this);

    satur->setAdjusterListener(this);

    estop->setAdjusterListener(this);

    scaltm->setAdjusterListener(this);

    rewei->setAdjusterListener(this);

    softradiustm->setLogScale(10, 0);
    softradiustm->setAdjusterListener(this);

    sensitm->setAdjusterListener(this);

    setExpandAlignProperties(expmasktm, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmasktmMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasktmMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmasktmMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmasktmMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmasktmMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmasktmMethod->set_active(0);
    showmasktmMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmasktmMethodConn = showmasktmMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabTone::showmasktmMethodChanged));

    enatmMaskConn = enatmMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabTone::enatmMaskChanged));

    enatmMaskaftConn = enatmMaskaft->signal_toggled().connect(sigc::mem_fun(*this, &LocallabTone::enatmMaskaftChanged));

    masktmCurveEditorG->setCurveListener(this);

    CCmasktmshape->setIdentityValue(0.);
    CCmasktmshape->setResetCurve(FlatCurveType(defSpot.CCmasktmcurve.at(0)), defSpot.CCmasktmcurve);
    CCmasktmshape->setBottomBarColorProvider(this, 1);

    LLmasktmshape->setIdentityValue(0.);
    LLmasktmshape->setResetCurve(FlatCurveType(defSpot.LLmasktmcurve.at(0)), defSpot.LLmasktmcurve);
    LLmasktmshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmasktmshape->setIdentityValue(0.);
    HHmasktmshape->setResetCurve(FlatCurveType(defSpot.HHmasktmcurve.at(0)), defSpot.HHmasktmcurve);
    HHmasktmshape->setCurveColorProvider(this, 2);
    HHmasktmshape->setBottomBarColorProvider(this, 2);

    masktmCurveEditorG->curveListComplete();

    blendmasktm->setAdjusterListener(this);

    lapmasktm->setAdjusterListener(this);

    radmasktm->setAdjusterListener(this);

    chromasktm->setAdjusterListener(this);

    gammasktm->setAdjusterListener(this);

    slomasktm->setAdjusterListener(this);

    mask2tmCurveEditorG->setCurveListener(this);
    Lmasktmshape->setResetCurve(DiagonalCurveType(defSpot.Lmasktmcurve.at(0)), defSpot.Lmasktmcurve);
    Lmasktmshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmasktmshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2tmCurveEditorG->curveListComplete();

    // Add Tone Mapping specific widgets to GUI
    // pack_start(*amount); // To use if we change transit_shapedetect parameters
    pack_start(*stren);
    pack_start(*equiltm);
    pack_start(*gamma);
    pack_start(*satur);
    pack_start(*estop);
    pack_start(*scaltm);
    pack_start(*rewei);
    // pack_start(*softradiustm); // Always bad with TM ??
    pack_start(*sensitm);
    ToolParamBlock* const masktmBox = Gtk::manage(new ToolParamBlock());
    masktmBox->pack_start(*showmasktmMethod, Gtk::PACK_SHRINK, 4);
    masktmBox->pack_start(*enatmMask, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*enatmMaskaft, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*masktmCurveEditorG, Gtk::PACK_SHRINK, 4);
    masktmBox->pack_start(*blendmasktm, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*lapmasktm, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*radmasktm, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*chromasktm, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*gammasktm, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*slomasktm, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*mask2tmCurveEditorG, Gtk::PACK_SHRINK, 4);
    expmasktm->add(*masktmBox, false);
    pack_start(*expmasktm, false, false);
}

LocallabTone::~LocallabTone()
{
    delete masktmCurveEditorG;
    delete mask2tmCurveEditorG;
}

bool LocallabTone::isMaskViewActive()
{
    return (showmasktmMethod->get_active_row_number() != 0);
}

void LocallabTone::resetMaskView()
{
    showmasktmMethodConn.block(true);
    showmasktmMethod->set_active(0);
    showmasktmMethodConn.block(false);
}

void LocallabTone::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &maskMask)
{
    tmMask = showmasktmMethod->get_active_row_number();
}

void LocallabTone::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_TONEMAP_TOOLTIP"));
        equiltm->set_tooltip_text(M("TP_LOCALLAB_EQUILTM_TOOLTIP"));
        gamma->set_tooltip_text(M("TP_LOCALLAB_TONEMAPGAM_TOOLTIP"));
        estop->set_tooltip_text(M("TP_LOCALLAB_TONEMAPESTOP_TOOLTIP"));
        scaltm->set_tooltip_text(M("TP_LOCALLAB_TONEMASCALE_TOOLTIP"));
        rewei->set_tooltip_text(M("TP_LOCALLAB_TONEMAPREWEI_TOOLTIP"));
        sensitm->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        expmasktm->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmasktmshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmasktmshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmasktmshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        blendmasktm->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        lapmasktm->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        radmasktm->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        mask2tmCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        Lmasktmshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    } else {
        exp->set_tooltip_text("");
        equiltm->set_tooltip_text("");
        gamma->set_tooltip_text("");
        estop->set_tooltip_text("");
        scaltm->set_tooltip_text("");
        rewei->set_tooltip_text("");
        sensitm->set_tooltip_text("");
        expmasktm->set_tooltip_text("");
        CCmasktmshape->setTooltip("");
        LLmasktmshape->setTooltip("");
        HHmasktmshape->setTooltip("");
        blendmasktm->set_tooltip_text("");
        lapmasktm->set_tooltip_text("");
        radmasktm->set_tooltip_text("");
        mask2tmCurveEditorG->set_tooltip_text("");
        Lmasktmshape->setTooltip("");
    }
}

void LocallabTone::setDefaultExpanderVisibility()
{
    expmasktm->set_expanded(false);
}

void LocallabTone::disableListener()
{
    LocallabTool::disableListener();

    equiltmConn.block(true);
    showmasktmMethodConn.block(true);
    enatmMaskConn.block(true);
    enatmMaskaftConn.block(true);
}

void LocallabTone::enableListener()
{
    LocallabTool::enableListener();

    equiltmConn.block(false);
    showmasktmMethodConn.block(false);
    enatmMaskConn.block(false);
    enatmMaskaftConn.block(false);
}

void LocallabTone::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visitonemap);
        exp->setEnabled(spot.exptonemap);
        complexity->set_active(spot.complextonemap);

        amount->setValue(spot.amount);
        stren->setValue(spot.stren);
        equiltm->set_active(spot.equiltm);
        gamma->setValue(spot.gamma);
        satur->setValue(spot.satur);
        estop->setValue(spot.estop);
        scaltm->setValue(spot.scaltm);
        rewei->setValue((double)spot.rewei);
        softradiustm->setValue(spot.softradiustm);
        sensitm->setValue((double)spot.sensitm);
        enatmMask->set_active(spot.enatmMask);
        enatmMaskaft->set_active(spot.enatmMaskaft);
        CCmasktmshape->setCurve(spot.CCmasktmcurve);
        LLmasktmshape->setCurve(spot.LLmasktmcurve);
        HHmasktmshape->setCurve(spot.HHmasktmcurve);
        blendmasktm->setValue((double)spot.blendmasktm);
        lapmasktm->setValue(spot.lapmasktm);
        radmasktm->setValue(spot.radmasktm);
        chromasktm->setValue(spot.chromasktm);
        gammasktm->setValue(spot.gammasktm);
        slomasktm->setValue(spot.slomasktm);
        Lmasktmshape->setCurve(spot.Lmasktmcurve);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabTone::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.exptonemap = exp->getEnabled();
        spot.visitonemap = exp->get_visible();
        spot.complextonemap = complexity->get_active_row_number();

        spot.amount = amount->getValue();
        spot.stren = stren->getValue();
        spot.equiltm = equiltm->get_active();
        spot.gamma = gamma->getValue();
        spot.satur = satur->getValue();
        spot.estop = estop->getValue();
        spot.scaltm = scaltm->getValue();
        spot.rewei = rewei->getIntValue();
        spot.softradiustm = softradiustm->getValue();
        spot.sensitm = sensitm->getIntValue();
        spot.enatmMask = enatmMask->get_active();
        spot.enatmMaskaft = enatmMaskaft->get_active();
        spot.LLmasktmcurve = LLmasktmshape->getCurve();
        spot.CCmasktmcurve = CCmasktmshape->getCurve();
        spot.HHmasktmcurve = HHmasktmshape->getCurve();
        spot.blendmasktm = blendmasktm->getIntValue();
        spot.lapmasktm = lapmasktm->getValue();
        spot.radmasktm = radmasktm->getValue();
        spot.chromasktm = chromasktm->getValue();
        spot.gammasktm = gammasktm->getValue();
        spot.slomasktm = slomasktm->getValue();
        spot.Lmasktmcurve = Lmasktmshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabTone::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        amount->setDefault(defSpot.amount);
        stren->setDefault(defSpot.stren);
        gamma->setDefault(defSpot.gamma);
        satur->setDefault(defSpot.satur);
        estop->setDefault(defSpot.estop);
        scaltm->setDefault(defSpot.scaltm);
        rewei->setDefault((double)defSpot.rewei);
        softradiustm->setDefault(defSpot.softradiustm);
        sensitm->setDefault((double)defSpot.sensitm);
        blendmasktm->setDefault((double)defSpot.blendmasktm);
        lapmasktm->setDefault(defSpot.lapmasktm);
        radmasktm->setDefault(defSpot.radmasktm);
        chromasktm->setDefault(defSpot.chromasktm);
        gammasktm->setDefault(defSpot.gammasktm);
        slomasktm->setDefault(defSpot.slomasktm);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabTone::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled() && listener) {
        const auto spName = " (" + escapeHtmlChars(spotName) + ")";

        if (a == amount) {
            listener->panelChanged(Evlocallabamount, amount->getTextValue() + spName);
        } else if (a == stren) {
            listener->panelChanged(Evlocallabstren, stren->getTextValue() + spName);
        } else if (a == gamma) {
            listener->panelChanged(Evlocallabgamma, gamma->getTextValue() + spName);
        } else if (a == satur) {
            listener->panelChanged(Evlocallabsatur, satur->getTextValue() + spName);
        } else if (a == estop) {
            listener->panelChanged(Evlocallabestop, estop->getTextValue() + spName);
        } else if (a == scaltm) {
            listener->panelChanged(Evlocallabscaltm, scaltm->getTextValue() + spName);
        } else if (a == rewei) {
            listener->panelChanged(Evlocallabrewei, rewei->getTextValue() + spName);
        } else if (a == softradiustm) {
            listener->panelChanged(Evlocallabsoftradiustm, softradiustm->getTextValue() + spName);
        } else if (a == sensitm) {
            listener->panelChanged(Evlocallabsensitm, sensitm->getTextValue() + spName);
        } else if (a == blendmasktm) {
            listener->panelChanged(Evlocallabblendmasktm, blendmasktm->getTextValue() + spName);
        } else if (a == lapmasktm) {
            listener->panelChanged(Evlocallablapmasktm, lapmasktm->getTextValue() + spName);
        } else if (a == radmasktm) {
            listener->panelChanged(Evlocallabradmasktm, radmasktm->getTextValue() + spName);
        } else if (a == chromasktm) {
            listener->panelChanged(Evlocallabchromasktm, chromasktm->getTextValue() + spName);
        } else if (a == gammasktm) {
            listener->panelChanged(Evlocallabgammasktm, gammasktm->getTextValue() + spName);
        } else if (a == slomasktm) {
            listener->panelChanged(Evlocallabslomasktm, slomasktm->getTextValue() + spName);
        }
    }
}

void LocallabTone::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled() && listener) {
        const auto spName = M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")";

        if (ce == CCmasktmshape) {
            listener->panelChanged(EvlocallabCCmasktmshape, spName);
        } else if (ce == LLmasktmshape) {
            listener->panelChanged(EvlocallabLLmasktmshape, spName);
        } else if (ce == HHmasktmshape) {
            listener->panelChanged(EvlocallabHHmasktmshape, spName);
        } else if (ce == Lmasktmshape) {
            listener->panelChanged(EvlocallabLmasktmshape, spName);
        }
    }
}

void LocallabTone::enabledChanged()
{
    if (isLocActivated && listener) {
        listener->panelChanged(EvLocenatonemap, (exp->getEnabled() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"))
                               + " (" + escapeHtmlChars(spotName) + ")");
    }
}

void LocallabTone::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    gamma->setValue(defSpot.gamma);
    satur->setValue(defSpot.satur);
    rewei->setValue((double)defSpot.rewei);
    lapmasktm->setValue(defSpot.lapmasktm);
    gammasktm->setValue(defSpot.gammasktm);
    slomasktm->setValue(defSpot.slomasktm);

    // Enable all listeners
    enableListener();
}

void LocallabTone::convertParamToSimple()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden specific GUI widgets in Simple mode to default spot values
    showmasktmMethod->set_active(0);
    enatmMask->set_active(defSpot.enatmMask);
    enatmMaskaft->set_active(defSpot.enatmMaskaft);
    CCmasktmshape->setCurve(defSpot.CCmasktmcurve);
    LLmasktmshape->setCurve(defSpot.LLmasktmcurve);
    HHmasktmshape->setCurve(defSpot.HHmasktmcurve);
    blendmasktm->setValue((double)defSpot.blendmasktm);
    radmasktm->setValue(defSpot.radmasktm);
    chromasktm->setValue(defSpot.chromasktm);
    Lmasktmshape->setCurve(defSpot.Lmasktmcurve);

    // Enable all listeners
    enableListener();
}

void LocallabTone::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            gamma->hide();
            satur->hide();
            rewei->hide();
            expmasktm->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            gamma->hide();
            satur->hide();
            rewei->hide();
            lapmasktm->hide();
            gammasktm->hide();
            slomasktm->hide();
            // Specific Simple mode widgets are shown in Normal mode
            expmasktm->show();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            gamma->show();
            satur->show();
            rewei->show();
            expmasktm->show();
            lapmasktm->show();
            gammasktm->show();
            slomasktm->show();
    }
}

void LocallabTone::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmasktmshape->updateLocallabBackground(normChromar);
        LLmasktmshape->updateLocallabBackground(normLumar);
        HHmasktmshape->updateLocallabBackground(normHuer);
        Lmasktmshape->updateLocallabBackground(normLumar);

        return false;
    }
    );
}

void LocallabTone::equiltmChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (equiltm->get_active()) {
                listener->panelChanged(Evlocallabequiltm,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabequiltm,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabTone::showmasktmMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabTone::enatmMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enatmMask->get_active()) {
                listener->panelChanged(EvLocallabEnatmMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnatmMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabTone::enatmMaskaftChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enatmMaskaft->get_active()) {
                listener->panelChanged(EvLocallabEnatmMaskaft,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnatmMaskaft,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

/* ==== LocallabRetinex ==== */
LocallabRetinex::LocallabRetinex():
    LocallabTool(this, M("TP_LOCALLAB_RET_TOOLNAME"), M("TP_LOCALLAB_RETI"), true),

    // Retinex specific widgets
    dehaz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DEHAZ"), -100, 100, 1, 0))),
    depth(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DEPTH"), 0, 100, 1, 25))),
    lumonly(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_LUMONLY")))),
    retiFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RETIFRA")))),
    str(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STR"), 0., 100., 0.2, 0.))),
    loglin(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_LOGLIN")))),
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIH"), 0, 100, 1, 60))),
    retitoolFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RETITOOLFRA")))),
    retinexMethod(Gtk::manage(new MyComboBoxText())),
    fftwreti(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FFTW")))),
    equilret(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_EQUIL")))),
    neigh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NEIGH"), MINNEIGH, MAXNEIGH, 0.5, 50., nullptr, nullptr, &retiSlider2neigh, &retiNeigh2Slider))),
    vart(Gtk::manage(new Adjuster(M("TP_LOCALLAB_VART"), 0.1, 500., 0.1, 150.))),
    scalereti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALERETI"), 1.0, 10.0, 1., 2.))),
    limd(Gtk::manage(new Adjuster(M("TP_LOCALLAB_THRESRETI"), 1.2, 100.0, 0.1, 8.))),
    offs(Gtk::manage(new Adjuster(M("TP_LOCALLAB_OFFS"), -16386., 32768., 1., 0.))),
    expretitools(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPRETITOOLS")))),
    chrrt(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHRRT"), 0.0, 100.0, 0.1, 0.0))),
    darkness(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DARKRETI"), 0.01, 6.0, 0.01, 2.0))),
    lightnessreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LIGHTRETI"), 0.01, 4.0, 0.01, 1.))),
    cliptm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLIPTM"), 0.02, 2.0, 0.01, 1.))),
    softradiusret(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRETI"), 0.0, 100.0, 0.5, 40.))),
    LocalcurveEditortransT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONMAP"))),
    cTtransshape(static_cast<FlatCurveEditor*>(LocalcurveEditortransT->addCurve(CT_Flat, "", nullptr, false, false))),
    mMLabels(Gtk::manage(new Gtk::Label("---"))),
    transLabels(Gtk::manage(new Gtk::Label("---"))),
    transLabels2(Gtk::manage(new Gtk::Label("---"))),
    LocalcurveEditorgainT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAIN"))),
    cTgainshape(static_cast<FlatCurveEditor*>(LocalcurveEditorgainT->addCurve(CT_Flat, "", nullptr, false, false))),
    expmaskreti(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWR")))),
    showmaskretiMethod(Gtk::manage(new MyComboBoxText())),
    enaretiMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    enaretiMasktmap(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TM_MASK")))),
    maskretiCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskretishape(static_cast<FlatCurveEditor*>(maskretiCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskretishape(static_cast<FlatCurveEditor*>(maskretiCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskretishape(static_cast<FlatCurveEditor *>(maskretiCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 10.))),
    lapmaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.05, 5.0, 0.01, 1.))),
    slomaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2retiCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmaskretishape(static_cast<DiagonalCurveEditor*>(mask2retiCurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    inversret(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS"))))
{
    const LocallabParams::LocallabSpot defSpot;

    // Parameter Retinex specific widgets
    dehaz->set_tooltip_text(M("TP_LOCALLAB_DEHAZ_TOOLTIP"));
    dehaz->setAdjusterListener(this);

    depth->setAdjusterListener(this);

    lumonlyConn = lumonly->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::lumonlyChanged));

    retiFrame->set_label_align(0.025, 0.5);

    str->setAdjusterListener(this);

    loglinConn = loglin->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::loglinChanged));

    sensih->setAdjusterListener(this);

    retitoolFrame->set_label_align(0.025, 0.5);

    retinexMethod->append(M("TP_RETINEX_LOW"));
    retinexMethod->append(M("TP_RETINEX_UNIFORM"));
    retinexMethod->append(M("TP_RETINEX_HIGH"));
    retinexMethod->set_active(0);
    retinexMethod->set_tooltip_markup(M("TP_LOCRETI_METHOD_TOOLTIP"));
    retinexMethodConn = retinexMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabRetinex::retinexMethodChanged));

    fftwretiConn = fftwreti->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::fftwretiChanged));

    equilretConn = equilret->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::equilretChanged));

    neigh->setAdjusterListener(this);

    vart->setAdjusterListener(this);

    scalereti->setAdjusterListener(this);

    limd->setAdjusterListener(this);

    offs->setAdjusterListener(this);

    setExpandAlignProperties(expretitools, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    chrrt->setAdjusterListener(this);

    darkness->setAdjusterListener(this);

    lightnessreti->setAdjusterListener(this);

    cliptm->setAdjusterListener(this);

    softradiusret->setLogScale(10, 0);
    softradiusret->setAdjusterListener(this);

    LocalcurveEditortransT->setCurveListener(this);

    cTtransshape->setIdentityValue(0.);
    cTtransshape->setResetCurve(FlatCurveType(defSpot.localTtranscurve.at(0)), defSpot.localTtranscurve);

    LocalcurveEditortransT->curveListComplete();

    setExpandAlignProperties(mMLabels, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);

    setExpandAlignProperties(transLabels, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);

    setExpandAlignProperties(transLabels2, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);

    LocalcurveEditorgainT->setCurveListener(this);

    cTgainshape->setIdentityValue(0.);
    cTgainshape->setResetCurve(FlatCurveType(defSpot.localTgaincurve.at(0)), defSpot.localTgaincurve);

    LocalcurveEditorgainT->curveListComplete();

    setExpandAlignProperties(expmaskreti, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskretiMethod->set_active(0);
    showmaskretiMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskretiMethodConn = showmaskretiMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabRetinex::showmaskretiMethodChanged));

    enaretiMaskConn = enaretiMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::enaretiMaskChanged));

    enaretiMasktmapConn = enaretiMasktmap->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::enaretiMasktmapChanged));

    maskretiCurveEditorG->setCurveListener(this);

    CCmaskretishape->setIdentityValue(0.);
    CCmaskretishape->setResetCurve(FlatCurveType(defSpot.CCmaskreticurve.at(0)), defSpot.CCmaskreticurve);
    CCmaskretishape->setBottomBarColorProvider(this, 1);

    LLmaskretishape->setIdentityValue(0.);
    LLmaskretishape->setResetCurve(FlatCurveType(defSpot.LLmaskreticurve.at(0)), defSpot.LLmaskreticurve);
    LLmaskretishape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskretishape->setIdentityValue(0.);
    HHmaskretishape->setResetCurve(FlatCurveType(defSpot.HHmaskreticurve.at(0)), defSpot.HHmaskreticurve);
    HHmaskretishape->setCurveColorProvider(this, 2);
    HHmaskretishape->setBottomBarColorProvider(this, 2);

    maskretiCurveEditorG->curveListComplete();

    blendmaskreti->setAdjusterListener(this);

    radmaskreti->setAdjusterListener(this);

    lapmaskreti->setAdjusterListener(this);

    chromaskreti->setAdjusterListener(this);

    gammaskreti->setAdjusterListener(this);

    slomaskreti->setAdjusterListener(this);

    mask2retiCurveEditorG->setCurveListener(this);

    Lmaskretishape->setResetCurve(DiagonalCurveType(defSpot.Lmaskreticurve.at(0)), defSpot.Lmaskreticurve);
    Lmaskretishape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskretishape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2retiCurveEditorG->curveListComplete();

    inversretConn = inversret->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::inversretChanged));

    // Add Retinex specific widgets to GUI
    ToolParamBlock* const auxBox = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const dehaFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_DEHAFRA")));
    dehaFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const dehaBox = Gtk::manage(new ToolParamBlock());
    dehaBox->pack_start(*dehaz);
    dehaBox->pack_start(*depth);
    dehaBox->pack_start(*lumonly);
    dehaFrame->add(*dehaBox);
    auxBox->add(*dehaFrame);
    ToolParamBlock* const deharetiBox = Gtk::manage(new ToolParamBlock());
    deharetiBox->pack_start(*str);
    deharetiBox->pack_start(*loglin);
    retiFrame->add(*deharetiBox);
    auxBox->add(*retiFrame);
    ToolParamBlock* const scopeBox = Gtk::manage(new ToolParamBlock());
    scopeBox->pack_start(*sensih);
    auxBox->add(*scopeBox);
    pack_start(*auxBox);
    ToolParamBlock* const retiBox = Gtk::manage(new ToolParamBlock());
    retiBox->pack_start(*retinexMethod);
    retiBox->pack_start(*fftwreti);
    retiBox->pack_start(*equilret);
    retiBox->pack_start(*neigh);
    retiBox->pack_start(*vart);
    retiBox->pack_start(*scalereti);
    retiBox->pack_start(*limd);
    retiBox->pack_start(*offs);
    ToolParamBlock* const toolretiBox = Gtk::manage(new ToolParamBlock());
    toolretiBox->pack_start(*chrrt);
    toolretiBox->pack_start(*darkness);
    toolretiBox->pack_start(*lightnessreti);
    toolretiBox->pack_start(*cliptm);
    toolretiBox->pack_start(*softradiusret);
    toolretiBox->pack_start(*LocalcurveEditortransT, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolretiBox->pack_start(*mMLabels);
    toolretiBox->pack_start(*transLabels);
    toolretiBox->pack_start(*transLabels2);
    toolretiBox->pack_start(*LocalcurveEditorgainT, Gtk::PACK_SHRINK, 4);
    expretitools->add(*toolretiBox, false);
    retiBox->pack_start(*expretitools, false, false);
    ToolParamBlock* const maskretiBox = Gtk::manage(new ToolParamBlock());
    maskretiBox->pack_start(*showmaskretiMethod, Gtk::PACK_SHRINK, 4);
    maskretiBox->pack_start(*enaretiMask, Gtk::PACK_SHRINK, 0);
    maskretiBox->pack_start(*enaretiMasktmap, Gtk::PACK_SHRINK, 0);
    maskretiBox->pack_start(*maskretiCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskretiBox->pack_start(*blendmaskreti, Gtk::PACK_SHRINK, 0);
    maskretiBox->pack_start(*radmaskreti, Gtk::PACK_SHRINK, 0);
    maskretiBox->pack_start(*lapmaskreti, Gtk::PACK_SHRINK, 0);
    maskretiBox->pack_start(*chromaskreti, Gtk::PACK_SHRINK, 0);
    maskretiBox->pack_start(*gammaskreti, Gtk::PACK_SHRINK, 0);
    maskretiBox->pack_start(*slomaskreti, Gtk::PACK_SHRINK, 0);
    maskretiBox->pack_start(*mask2retiCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    expmaskreti->add(*maskretiBox, false);
    retiBox->pack_start(*expmaskreti, false, false);
    // retiBox->pack_start(*inversret);
    retitoolFrame->add(*retiBox);
    pack_start(*retitoolFrame);
}

LocallabRetinex::~LocallabRetinex()
{
    delete LocalcurveEditortransT;
    delete LocalcurveEditorgainT;
    delete maskretiCurveEditorG;
    delete mask2retiCurveEditorG;
}

void LocallabRetinex::updateMinMax(const double cdma, const double cdmin, const double mini, const double maxi, const double Tmean, const double Tsigma, const double Tmin, const double Tmax)
{
    idle_register.add(
    [this, cdma, cdmin, mini, maxi, Tmean, Tsigma, Tmin, Tmax]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        mMLabels->set_text(
            Glib::ustring::compose(M("TP_LOCALLAB_MLABEL"),
                                   Glib::ustring::format(std::fixed, std::setprecision(0), cdmin),
                                   Glib::ustring::format(std::fixed, std::setprecision(0), cdma))
        );
        transLabels->set_text(
            Glib::ustring::compose(M("TP_LOCALLAB_TLABEL"),
                                   Glib::ustring::format(std::fixed, std::setprecision(1), mini),
                                   Glib::ustring::format(std::fixed, std::setprecision(1), maxi),
                                   Glib::ustring::format(std::fixed, std::setprecision(1), Tmean),
                                   Glib::ustring::format(std::fixed, std::setprecision(1), Tsigma))
        );
        transLabels2->set_text(
            Glib::ustring::compose(M("TP_RETINEX_TLABEL2"),
                                   Glib::ustring::format(std::fixed, std::setprecision(1), Tmin),
                                   Glib::ustring::format(std::fixed, std::setprecision(1), Tmax))
        );

        return false;
    }
    );
}

bool LocallabRetinex::isMaskViewActive()
{
    return (showmaskretiMethod->get_active_row_number() != 0);
}

void LocallabRetinex::resetMaskView()
{
    showmaskretiMethodConn.block(true);
    showmaskretiMethod->set_active(0);
    showmaskretiMethodConn.block(false);
}

void LocallabRetinex::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &maskMask)
{
    retiMask = showmaskretiMethod->get_active_row_number();
}

void LocallabRetinex::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        loglin->set_tooltip_text(M("TP_LOCALLAB_RETI_LOGLIN_TOOLTIP"));
        sensih->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
        fftwreti->set_tooltip_text(M("TP_LOCALLAB_RETI_FFTW_TOOLTIP"));
        equilret->set_tooltip_text(M("TP_LOCALLAB_EQUILTM_TOOLTIP"));
        neigh->set_tooltip_text(M("TP_LOCALLAB_RETI_NEIGH_VART_TOOLTIP"));
        vart->set_tooltip_text(M("TP_LOCALLAB_RETI_NEIGH_VART_TOOLTIP"));
        scalereti->set_tooltip_text(M("TP_LOCALLAB_RETI_SCALE_TOOLTIP"));
        limd->set_tooltip_text(M("TP_LOCALLAB_RETI_LIMDOFFS_TOOLTIP"));
        offs->set_tooltip_text(M("TP_LOCALLAB_RETI_LIMDOFFS_TOOLTIP"));
        darkness->set_tooltip_text(M("TP_LOCALLAB_RETI_LIGHTDARK_TOOLTIP"));
        lightnessreti->set_tooltip_text(M("TP_LOCALLAB_RETI_LIGHTDARK_TOOLTIP"));
        cliptm->set_tooltip_text(M("TP_LOCALLAB_RETI_LIMDOFFS_TOOLTIP"));
        softradiusret->set_tooltip_text(M("TP_LOCALLAB_GUIDFILTER_TOOLTIP"));
        cTtransshape->setTooltip(M("TP_LOCALLAB_TRANSMISSION_TOOLTIP"));
        mMLabels->set_tooltip_markup(M("TP_LOCALLAB_MLABEL_TOOLTIP"));
        transLabels->set_tooltip_markup(M("TP_LOCALLAB_TLABEL_TOOLTIP"));
        cTgainshape->setTooltip(M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));
        expmaskreti->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        enaretiMasktmap->set_tooltip_markup(M("TP_LOCALLAB_ENARETIMASKTMAP_TOOLTIP"));
        CCmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        blendmaskreti->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        radmaskreti->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        lapmaskreti->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        mask2retiCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        Lmaskretishape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    } else {
        loglin->set_tooltip_text("");
        sensih->set_tooltip_text("");
        fftwreti->set_tooltip_text("");
        equilret->set_tooltip_text("");
        neigh->set_tooltip_text("");
        vart->set_tooltip_text("");
        scalereti->set_tooltip_text("");
        limd->set_tooltip_text("");
        offs->set_tooltip_text("");
        darkness->set_tooltip_text("");
        lightnessreti->set_tooltip_text("");
        cliptm->set_tooltip_text("");
        softradiusret->set_tooltip_text("");
        cTtransshape->setTooltip("");
        mMLabels->set_tooltip_text("");
        transLabels->set_tooltip_text("");
        cTgainshape->setTooltip("");
        expmaskreti->set_tooltip_text("");
        enaretiMasktmap->set_tooltip_text("");
        CCmaskretishape->setTooltip("");
        LLmaskretishape->setTooltip("");
        HHmaskretishape->setTooltip("");
        blendmaskreti->set_tooltip_text("");
        radmaskreti->set_tooltip_text("");
        lapmaskreti->set_tooltip_text("");
        mask2retiCurveEditorG->set_tooltip_text("");
        Lmaskretishape->setTooltip("");
    }
}

void LocallabRetinex::setDefaultExpanderVisibility()
{
    expretitools->set_expanded(false);
    expmaskreti->set_expanded(false);
}

void LocallabRetinex::disableListener()
{
    LocallabTool::disableListener();

    lumonlyConn.block(true);
    loglinConn.block(true);
    retinexMethodConn.block(true);
    fftwretiConn.block(true);
    equilretConn.block(true);
    showmaskretiMethodConn.block(true);
    enaretiMaskConn.block(true);
    enaretiMasktmapConn.block(true);
    inversretConn.block(true);
}

void LocallabRetinex::enableListener()
{
    LocallabTool::enableListener();

    lumonlyConn.block(false);
    loglinConn.block(false);
    retinexMethodConn.block(false);
    fftwretiConn.block(false);
    equilretConn.block(false);
    showmaskretiMethodConn.block(false);
    enaretiMaskConn.block(false);
    enaretiMasktmapConn.block(false);
    inversretConn.block(false);
}

void LocallabRetinex::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visireti);
        exp->setEnabled(spot.expreti);
        complexity->set_active(spot.complexreti);

        dehaz->setValue((double)spot.dehaz);
        depth->setValue((double)spot.depth);
        lumonly->set_active(spot.lumonly);
        str->setValue(spot.str);
        loglin->set_active(spot.loglin);
        sensih->setValue((double)spot.sensih);

        if (spot.retinexMethod == "low") {
            retinexMethod->set_active(0);
        } else if (spot.retinexMethod == "uni") {
            retinexMethod->set_active(1);
        } else {
            retinexMethod->set_active(2);
        }

        fftwreti->set_active(spot.fftwreti);
        equilret->set_active(spot.equilret);
        neigh->setValue(spot.neigh);
        vart->setValue(spot.vart);
        scalereti->setValue(spot.scalereti);
        limd->setValue(spot.limd);
        offs->setValue(spot.offs);
        chrrt->setValue(spot.chrrt);
        darkness->setValue(spot.darkness);
        lightnessreti->setValue(spot.lightnessreti);
        cliptm->setValue(spot.cliptm);
        softradiusret->setValue(spot.softradiusret);
        cTtransshape->setCurve(spot.localTtranscurve);
        cTgainshape->setCurve(spot.localTgaincurve);
        enaretiMask->set_active(spot.enaretiMask);
        enaretiMasktmap->set_active(spot.enaretiMasktmap);
        CCmaskretishape->setCurve(spot.CCmaskreticurve);
        LLmaskretishape->setCurve(spot.LLmaskreticurve);
        HHmaskretishape->setCurve(spot.HHmaskreticurve);
        blendmaskreti->setValue((double)spot.blendmaskreti);
        radmaskreti->setValue(spot.radmaskreti);
        lapmaskreti->setValue(spot.lapmaskreti);
        chromaskreti->setValue(spot.chromaskreti);
        gammaskreti->setValue(spot.gammaskreti);
        slomaskreti->setValue(spot.slomaskreti);
        Lmaskretishape->setCurve(spot.Lmaskreticurve);
        inversret->set_active(spot.inversret);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Update Retinex GUI according to scalereti adjuster value
    updateRetinexGUI1();

    // Update Retinex GUI according to inversret button state
    updateRetinexGUI2();

    // Update Retinex GUI according to str adjuster value
    updateRetinexGUI3();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabRetinex::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expreti = exp->getEnabled();
        spot.visireti = exp->get_visible();
        spot.complexreti = complexity->get_active_row_number();

        spot.dehaz = dehaz->getIntValue();
        spot.depth = depth->getIntValue();
        spot.lumonly = lumonly->get_active();
        spot.str = str->getValue();
        spot.loglin = loglin->get_active();
        spot.sensih = sensih->getIntValue();

        if (retinexMethod->get_active_row_number() == 0) {
            spot.retinexMethod = "low";
        } else if (retinexMethod->get_active_row_number() == 1) {
            spot.retinexMethod = "uni";
        } else if (retinexMethod->get_active_row_number() == 2) {
            spot.retinexMethod = "high";
        }

        spot.fftwreti = fftwreti->get_active();
        spot.equilret = equilret->get_active();
        spot.neigh = neigh->getValue();
        spot.vart = vart->getValue();
        spot.scalereti = scalereti->getValue();
        spot.limd = limd->getValue();
        spot.offs = offs->getValue();
        spot.chrrt = chrrt->getValue();
        spot.darkness = darkness->getValue();
        spot.lightnessreti = lightnessreti->getValue();
        spot.cliptm = cliptm->getValue();
        spot.softradiusret = softradiusret->getValue();
        spot.localTtranscurve = cTtransshape->getCurve();
        spot.localTgaincurve = cTgainshape->getCurve();
        spot.enaretiMask = enaretiMask->get_active();
        spot.enaretiMasktmap = enaretiMasktmap->get_active();
        spot.CCmaskreticurve = CCmaskretishape->getCurve();
        spot.LLmaskreticurve = LLmaskretishape->getCurve();
        spot.HHmaskreticurve = HHmaskretishape->getCurve();
        spot.blendmaskreti = blendmaskreti->getIntValue();
        spot.radmaskreti = radmaskreti->getValue();
        spot.lapmaskreti = lapmaskreti->getValue();
        spot.chromaskreti = chromaskreti->getValue();
        spot.gammaskreti = gammaskreti->getValue();
        spot.slomaskreti = slomaskreti->getValue();
        spot.Lmaskreticurve = Lmaskretishape->getCurve();
        spot.inversret = inversret->get_active();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabRetinex::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        dehaz->setDefault((double)defSpot.dehaz);
        depth->setDefault((double)defSpot.depth);
        str->setDefault(defSpot.str);
        sensih->setDefault((double)defSpot.sensih);
        neigh->setDefault(defSpot.neigh);
        vart->setDefault(defSpot.vart);
        scalereti->setDefault(defSpot.scalereti);
        limd->setDefault(defSpot.limd);
        offs->setDefault(defSpot.offs);
        chrrt->setDefault(defSpot.chrrt);
        darkness->setDefault(defSpot.darkness);
        lightnessreti->setDefault(defSpot.lightnessreti);
        cliptm->setDefault(defSpot.cliptm);
        softradiusret->setDefault(defSpot.softradiusret);
        blendmaskreti->setDefault((double)defSpot.blendmaskreti);
        radmaskreti->setDefault(defSpot.radmaskreti);
        lapmaskreti->setDefault(defSpot.lapmaskreti);
        chromaskreti->setDefault(defSpot.chromaskreti);
        gammaskreti->setDefault(defSpot.gammaskreti);
        slomaskreti->setDefault(defSpot.slomaskreti);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabRetinex::adjusterChanged(Adjuster* a, double newval)
{
    // Update Retinex GUI according to scalereti adjuster value
    if (a == scalereti) {
        updateRetinexGUI1();
    }

    // Update Retinex GUI according to str adjuster value
    if (a == str) {
        updateRetinexGUI3();
    }

    if (isLocActivated && exp->getEnabled()) {
        if (a == dehaz) {
            if (listener) {
                listener->panelChanged(Evlocallabdehaz,
                                       dehaz->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == depth) {
            if (listener) {
                listener->panelChanged(Evlocallabdepth,
                                       depth->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == str) {
            if (listener) {
                listener->panelChanged(Evlocallabstr,
                                       str->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensih) {
            if (listener) {
                listener->panelChanged(Evlocallabsensih,
                                       sensih->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == neigh) {
            if (listener) {
                listener->panelChanged(Evlocallabneigh,
                                       neigh->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == vart) {
            if (listener) {
                listener->panelChanged(Evlocallabvart,
                                       vart->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == scalereti) {
            if (listener) {
                listener->panelChanged(Evlocallabscalereti,
                                       scalereti->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == limd) {
            if (listener) {
                listener->panelChanged(Evlocallablimd,
                                       limd->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == offs) {
            if (listener) {
                listener->panelChanged(Evlocallaboffs,
                                       offs->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chrrt) {
            if (listener) {
                listener->panelChanged(Evlocallabchrrt,
                                       chrrt->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == darkness) {
            if (listener) {
                listener->panelChanged(Evlocallabdarkness,
                                       darkness->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lightnessreti) {
            if (listener) {
                listener->panelChanged(Evlocallablightnessreti,
                                       lightnessreti->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == cliptm) {
            if (listener) {
                listener->panelChanged(Evlocallabcliptm,
                                       cliptm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == softradiusret) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusret,
                                       softradiusret->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskreti,
                                       blendmaskreti->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskreti,
                                       radmaskreti->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskreti,
                                       lapmaskreti->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskreti,
                                       chromaskreti->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskreti,
                                       gammaskreti->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slomaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskreti,
                                       slomaskreti->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == cTtransshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCTtransCurve,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == cTgainshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCTgainCurve,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == CCmaskretishape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskretishape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmaskretishape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmaskretishape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenareti,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenareti,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    str->setValue(defSpot.str);
    loglin->set_active(defSpot.loglin);

    if (defSpot.retinexMethod == "low") {
        retinexMethod->set_active(0);
    } else if (defSpot.retinexMethod == "uni") {
        retinexMethod->set_active(1);
    } else {
        retinexMethod->set_active(2);
    }

    fftwreti->set_active(defSpot.fftwreti);
    equilret->set_active(defSpot.equilret);
    neigh->setValue(defSpot.neigh);
    vart->setValue(defSpot.vart);
    scalereti->setValue(defSpot.scalereti);
    limd->setValue(defSpot.limd);
    offs->setValue(defSpot.offs);
    chrrt->setValue(defSpot.chrrt);
    darkness->setValue(defSpot.darkness);
    lightnessreti->setValue(defSpot.lightnessreti);
    cliptm->setValue(defSpot.cliptm);
    softradiusret->setValue(defSpot.softradiusret);
    cTtransshape->setCurve(defSpot.localTtranscurve);
    cTgainshape->setCurve(defSpot.localTgaincurve);
    showmaskretiMethod->set_active(0);
    enaretiMask->set_active(defSpot.enaretiMask);
    enaretiMasktmap->set_active(defSpot.enaretiMasktmap);
    CCmaskretishape->setCurve(defSpot.CCmaskreticurve);
    LLmaskretishape->setCurve(defSpot.LLmaskreticurve);
    HHmaskretishape->setCurve(defSpot.HHmaskreticurve);
    blendmaskreti->setValue((double)defSpot.blendmaskreti);
    radmaskreti->setValue(defSpot.radmaskreti);
    lapmaskreti->setValue(defSpot.lapmaskreti);
    chromaskreti->setValue(defSpot.chromaskreti);
    gammaskreti->setValue(defSpot.gammaskreti);
    slomaskreti->setValue(defSpot.slomaskreti);
    Lmaskretishape->setCurve(defSpot.Lmaskreticurve);
    inversret->set_active(defSpot.inversret);

    // Enable all listeners
    enableListener();

    // Update GUI based on converted widget parameters:
    // - Update Retinex GUI according to scalereti adjuster value
    updateRetinexGUI1();
    // - Update Retinex GUI according to inversret button state
    updateRetinexGUI2();
    // - Update Retinex GUI according to str adjuster value
    updateRetinexGUI3();
}

void LocallabRetinex::convertParamToSimple()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden specific GUI widgets in Simple mode to default spot values
    lumonly->set_active(defSpot.lumonly);

    // Enable all listeners
    enableListener();
}

void LocallabRetinex::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            lumonly->hide();
            retiFrame->hide();
            retitoolFrame->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            retiFrame->hide();
            retitoolFrame->hide();
            // Specific Simple mode widgets are shown in Normal mode
            lumonly->show();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            lumonly->show();
            retiFrame->show();
            retitoolFrame->show();
    }
}

void LocallabRetinex::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskretishape->updateLocallabBackground(normChromar);
        LLmaskretishape->updateLocallabBackground(normLumar);
        HHmaskretishape->updateLocallabBackground(normHuer);
        Lmaskretishape->updateLocallabBackground(normLumar);

        return false;
    }
    );
}

void LocallabRetinex::lumonlyChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (lumonly->get_active()) {
                listener->panelChanged(Evlocallablumonly,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallablumonly,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::loglinChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (loglin->get_active()) {
                listener->panelChanged(Evlocallabloglin,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabloglin,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::retinexMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabretinexMethod,
                                   retinexMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabRetinex::fftwretiChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (fftwreti->get_active()) {
                listener->panelChanged(Evlocallabfftwreti,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabfftwreti,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::equilretChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (inversret->get_active()) {
                listener->panelChanged(Evlocallabequilret,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabequilret,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::showmaskretiMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabRetinex::enaretiMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaretiMask->get_active()) {
                listener->panelChanged(EvLocallabEnaretiMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaretiMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::enaretiMasktmapChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaretiMasktmap->get_active()) {
                listener->panelChanged(EvLocallabEnaretiMasktmap,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaretiMasktmap,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::inversretChanged()
{
    const bool maskPreviewActivated = isMaskViewActive();

    // Update Retinex GUI according to inversret button state
    updateRetinexGUI2();

    if (maskPreviewActivated) {
        // This event is called to transmit reset mask state
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
    }

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (inversret->get_active()) {
                listener->panelChanged(Evlocallabinversret,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabinversret,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::updateRetinexGUI1()
{
    // Update Retinex GUI according to scalereti adjuster value
    if (scalereti->getValue() == 1) {
        retinexMethod->hide();
        softradiusret->hide();
        LocalcurveEditortransT->hide();
        LocalcurveEditorgainT->hide();
    } else {
        retinexMethod->show();
        softradiusret->show();
        LocalcurveEditortransT->show();
        LocalcurveEditorgainT->show();
    }
}

void LocallabRetinex::updateRetinexGUI2()
{
    // Update Retinex GUI according to inversret button state
    if (inversret->get_active()) {
        expmaskreti->hide();
        showmaskretiMethodConn.block(true);
        showmaskretiMethod->set_active(0);
        showmaskretiMethodConn.block(false);
    } else {
        expmaskreti->show();
    }
}

void LocallabRetinex::updateRetinexGUI3()
{
    if (str->getValue() >= 0.1f) {
        retitoolFrame->show();
    } else {
        retitoolFrame->hide();
    }
}

/* ==== LocallabSharp ==== */
LocallabSharp::LocallabSharp():
    LocallabTool(this, M("TP_LOCALLAB_SHARP_TOOLNAME"), M("TP_LOCALLAB_SHARP"), true),

    // Sharpening specific widgets
    sharcontrast(Gtk::manage(new Adjuster(M("TP_SHARPENING_CONTRAST"), 0, 200, 1, 20))),
    sharblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARBLUR"), 0.2, 2.0, 0.05, 0.2))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 100))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 0))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 0.4, 2.5, 0.01, 0.75))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 40))),
    inverssha(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    sharFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHARFRAME")))),
    showmasksharMethod(Gtk::manage(new MyComboBoxText()))
{
    // Parameter Sharpening specific widgets
    sharcontrast->setAdjusterListener(this);

    sharradius->setAdjusterListener(this);

    sharamount->setAdjusterListener(this);

    shardamping->setAdjusterListener(this);

    shariter->setAdjusterListener(this);

    sharblur->setAdjusterListener(this);

    sensisha->setAdjusterListener(this);

    inversshaConn = inverssha->signal_toggled().connect(sigc::mem_fun(*this, &LocallabSharp::inversshaChanged));

    showmasksharMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasksharMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmasksharMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmasksharMethod->set_active(0);
    showmasksharMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmasksharMethodConn = showmasksharMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabSharp::showmasksharMethodChanged));

    // Add Sharpening specific widgets to GUI
    pack_start(*sharcontrast);
    pack_start(*sharblur);
    pack_start(*sharradius);
    pack_start(*sharamount);
    pack_start(*shardamping);
    pack_start(*shariter);
    pack_start(*sensisha);
    pack_start(*inverssha);
    sharFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const sharfBox = Gtk::manage(new ToolParamBlock());
    sharfBox->pack_start(*showmasksharMethod);
    sharFrame->add(*sharfBox);
    pack_start(*sharFrame);
}

bool LocallabSharp::isMaskViewActive()
{
    return (showmasksharMethod->get_active_row_number() != 0);
}

void LocallabSharp::resetMaskView()
{
    showmasksharMethodConn.block(true);
    showmasksharMethod->set_active(0);
    showmasksharMethodConn.block(false);
}

void LocallabSharp::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &maskMask)
{
    sharMask = showmasksharMethod->get_active_row_number();
}

void LocallabSharp::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPSHARP_TOOLTIP"));
        sensisha->set_tooltip_text(M("TP_LOCALLAB_SENSIS_TOOLTIP"));
    } else {
        exp->set_tooltip_text("");
        sensisha->set_tooltip_text("");
    }
}

void LocallabSharp::disableListener()
{
    LocallabTool::disableListener();

    inversshaConn.block(true);
    showmasksharMethodConn.block(true);
}

void LocallabSharp::enableListener()
{
    LocallabTool::enableListener();

    inversshaConn.block(false);
    showmasksharMethodConn.block(false);
}

void LocallabSharp::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visisharp);
        exp->setEnabled(spot.expsharp);
        complexity->set_active(spot.complexsharp);

        sharcontrast->setValue((double)spot.sharcontrast);
        sharradius->setValue(spot.sharradius);
        sharamount->setValue((double)spot.sharamount);
        shardamping->setValue((double)spot.shardamping);
        shariter->setValue((double)spot.shariter);
        sharblur->setValue(spot.sharblur);
        sensisha->setValue((double)spot.sensisha);
        inverssha->set_active(spot.inverssha);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSharp::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expsharp = exp->getEnabled();
        spot.visisharp = exp->get_visible();
        spot.complexsharp = complexity->get_active_row_number();

        spot.sharcontrast = sharcontrast->getIntValue();
        spot.sharradius = sharradius->getValue();
        spot.sharamount = sharamount->getIntValue();
        spot.shardamping = shardamping->getIntValue();
        spot.shariter = shariter->getIntValue();
        spot.sharblur = sharblur->getValue();
        spot.sensisha = sensisha->getIntValue();
        spot.inverssha = inverssha->get_active();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSharp::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        sharcontrast->setDefault((double)defSpot.sharcontrast);
        sharradius->setDefault(defSpot.sharradius);
        sharamount->setDefault((double)defSpot.sharamount);
        shardamping->setDefault((double)defSpot.shardamping);
        shariter->setDefault((double)defSpot.shariter);
        sharblur->setDefault(defSpot.sharblur);
        sensisha->setDefault((double)defSpot.sensisha);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSharp::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == sharcontrast) {
            if (listener) {
                listener->panelChanged(Evlocallabsharcontrast,
                                       sharcontrast->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sharradius) {
            if (listener) {
                listener->panelChanged(Evlocallabsharradius,
                                       sharradius->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sharamount) {
            if (listener) {
                listener->panelChanged(Evlocallabsharamount,
                                       sharamount->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == shardamping) {
            if (listener) {
                listener->panelChanged(Evlocallabshardamping,
                                       shardamping->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == shariter) {
            if (listener) {
                listener->panelChanged(Evlocallabshariter,
                                       shariter->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sharblur) {
            if (listener) {
                listener->panelChanged(Evlocallabsharblur,
                                       sharblur->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensisha) {
            if (listener) {
                listener->panelChanged(Evlocallabsensis,
                                       sensisha->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabSharp::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenasharp,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenasharp,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabSharp::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    sharcontrast->setValue((double)defSpot.sharcontrast);
    sharblur->setValue(defSpot.sharblur);
    sharamount->setValue(defSpot.sharamount);
    shardamping->setValue((double)defSpot.shardamping);
    shariter->setValue((double)defSpot.shariter);

    // Enable all listeners
    enableListener();
}

void LocallabSharp::convertParamToSimple()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden specific GUI widgets in Simple mode to default spot values
    showmasksharMethod->set_active(0);

    // Enable all listeners
    enableListener();

    // Update GUI based on converted widget parameters:
}

void LocallabSharp::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            sharcontrast->hide();
            sharblur->hide();
            sharamount->hide();
            shardamping->hide();
            shariter->hide();
            sharFrame->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            sharcontrast->hide();
            sharblur->hide();
            sharamount->hide();
            shardamping->hide();
            shariter->hide();
            // Specific Simple mode widgets are shown in Normal mode
            sharFrame->show();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            sharcontrast->show();
            sharblur->show();
            sharamount->show();
            shardamping->show();
            shariter->show();
            sharFrame->show();
    }
}

void LocallabSharp::inversshaChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (inverssha->get_active()) {
                listener->panelChanged(Evlocallabinverssha,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabinverssha,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabSharp::showmasksharMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

/* ==== LocallabContrast ==== */
LocallabContrast::LocallabContrast():
    LocallabTool(this, M("TP_LOCALLAB_LC_TOOLNAME"), M("TP_LOCALLAB_LOC_CONTRAST"), true),

    // Local contrast specific widgets
    localcontMethod(Gtk::manage(new MyComboBoxText())),
    lcradius(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_RADIUS"), 10, 100, 1, 80))),
    lcamount(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_AMOUNT"), 0, 1.0, 0.01, 0))),
    lcdarkness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_DARKNESS"), 0, 3.0, 0.01, 1.0))),
    lclightness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_LIGHTNESS"), 0, 3.0, 0.01, 1.0))),
    contFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CONTWFRA")))),
    sigmalc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    LocalcurveEditorwav(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAV"))),
    wavshape(static_cast<FlatCurveEditor*>(LocalcurveEditorwav->addCurve(CT_Flat, "", nullptr, false, false))),
    csThreshold(Gtk::manage(new ThresholdAdjuster(M("TP_LOCALLAB_CSTHRESHOLD"), 0, 9, 0, 0, 6, 6, 0, false))),
    levelwav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LEVELWAV"), 1, 9, 1, 4))),
    expresidpyr(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::HBox())))),
    residcont(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCONT"), -100, 100, 1, 0))),
    residchro(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCHRO"), -100., 100., 1., 0.))),
    residsha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDSHA"), -100., 100., 1., 0.))),
    residshathr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDSHATHR"), 0., 100., 1., 30.))),
    residhi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDHI"), -100., 100., 1., 0.))),
    residhithr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDHITHR"), 0., 100., 1., 70.))),
    sensilc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 60))),
    clariFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CLARIFRA")))),
    clarilres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARILRES"), -20., 100., 0.5, 0.))),
    claricres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARICRES"), -20., 100., 0.5, 0.))),
    clarisoft(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.5, 1.))),
    origlc(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ORIGLC")))),
    expcontrastpyr(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::HBox())))),
    wavgradl(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_GRALWFRA")))),
    sigmalc2(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    strwav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -4.0, 4.0, 0.05, 0.))),
    angwav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    wavedg(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_EDGFRA")))),
    strengthw(Gtk::manage(new Adjuster(M("TP_WAVELET_EDVAL"), 0., 100.0, 0.5, 0.))),
    sigmaed(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    LocalcurveEditorwavedg(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVEDG"))),
    wavshapeedg(static_cast<FlatCurveEditor*>(LocalcurveEditorwavedg->addCurve(CT_Flat, "", nullptr, false, false))),
    gradw(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGEDETECT"), 0., 100.0, 0.5, 90.))),
    waveshow(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_EDGSHOW")))),
    edgsBoxshow(Gtk::manage(new ToolParamBlock())),
    radiusw(Gtk::manage(new Adjuster(M("TP_WAVELET_EDRAD"), 5., 100.0, 0.5, 15.))),
    detailw(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGTHRESH"), -50., 100.0, 1., 10.))),
    localedgMethod(Gtk::manage(new MyComboBoxText())),
    tloww(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGEDETECTTHR"), 0., 100.0, 1., 20.))),
    thigw(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGEDETECTTHR2"), -10., 100.0, 1., 0.))),
    edgw(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGESENSI"), 0., 100.0, 1., 60.))),
    basew(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGEAMPLI"), 0., 100.0, 1., 10.))),
    localneiMethod(Gtk::manage(new MyComboBoxText())),
    wavblur(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_BLURLEVELFRA")))),
    levelblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LEVELBLUR"), 0., 100., 0.5, 0.))),
    sigmabl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    chromablu(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMABLU"), 0.0, 5., 0.1, 0.))),
    LocalcurveEditorwavlev(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVLEV"))),
    wavshapelev(static_cast<FlatCurveEditor*>(LocalcurveEditorwavlev->addCurve(CT_Flat, "", nullptr, false, false))),
    residblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDBLUR"), 0., 100., 0.5, 0.))),
    blurlc(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_BLURLC")))),
    expcontrastpyr2(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::HBox())))),
    wavcont(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CONTFRA")))),
    sigma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    offset(Gtk::manage(new Adjuster(M("TP_LOCALLAB_OFFSETWAV"), 0.33, 1.66, 0.01, 1., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    chromalev(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMALEV"), 0.1, 5., 0.1, 1.))),
    LocalcurveEditorwavcon(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVCON"))),
    wavshapecon(static_cast<FlatCurveEditor*>(LocalcurveEditorwavcon->addCurve(CT_Flat, "", nullptr, false, false))),
    wavcompre(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_COMPREFRA")))),
    LocalcurveEditorwavcompre(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVCOMPRE"))),
    wavshapecompre(static_cast<FlatCurveEditor*>(LocalcurveEditorwavcompre->addCurve(CT_Flat, "", nullptr, false, false))),
    sigmadr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    threswav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_THRESWAV"), 0.9, 2., 0.01, 1.4))),
    residcomp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCOMP"), -1., 1., 0.01, 0.))),
    wavcomp(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_COMPFRA")))),
    sigmadc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 3., 0.01, 1.))),
    deltad(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DELTAD"), -3., 3., 0.1, 0.))),//, Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    LocalcurveEditorwavcomp(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVCOMP"))),
    wavshapecomp(static_cast<FlatCurveEditor*>(LocalcurveEditorwavcomp->addCurve(CT_Flat, "", nullptr, false, false))),
    fatres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATRES"), 0., 100., 1., 0.))),
    fftwlc(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FFTW")))),
    expmasklc(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWLC")))),
    showmasklcMethod(Gtk::manage(new MyComboBoxText())),
    enalcMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    masklcCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmasklcshape(static_cast<FlatCurveEditor*>(masklcCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmasklcshape(static_cast<FlatCurveEditor*>(masklcCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmasklcshape(static_cast<FlatCurveEditor *>(masklcCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmasklc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmasklc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromasklc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    mask2lcCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmasklcshape(static_cast<DiagonalCurveEditor*>(mask2lcCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
{
    const LocallabParams::LocallabSpot defSpot;

    // Parameter Local contrast specific widgets
    localcontMethod->append(M("TP_LOCALLAB_LOCCONT"));
    localcontMethod->append(M("TP_LOCALLAB_WAVE"));
    localcontMethod->set_active(0);
    localcontMethodConn = localcontMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabContrast::localcontMethodChanged));

    lcradius->setAdjusterListener(this);

    lcamount->setAdjusterListener(this);

    lcdarkness->setAdjusterListener(this);

    lclightness->setAdjusterListener(this);

    contFrame->set_label_align(0.025, 0.5);

    sigmalc->setAdjusterListener(this);

    LocalcurveEditorwav->setCurveListener(this);

    wavshape->setIdentityValue(0.);
    wavshape->setResetCurve(FlatCurveType(defSpot.locwavcurve.at(0)), defSpot.locwavcurve);
    wavshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    LocalcurveEditorwav->curveListComplete();

    csThreshold->setAdjusterListener(this);

    levelwav->setAdjusterListener(this);

    Gtk::HBox* const LresTitleHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const LresLabel = Gtk::manage(new Gtk::Label());
    LresLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_RESIDPYR")) + Glib::ustring("</b>"));
    LresLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    LresTitleHBox->pack_start(*LresLabel, Gtk::PACK_EXPAND_WIDGET, 0);
    expresidpyr->setLabel(LresTitleHBox);
    setExpandAlignProperties(expresidpyr, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    residcont->setAdjusterListener(this);

    residchro->setAdjusterListener(this);

    residsha->setAdjusterListener(this);

    residshathr->setAdjusterListener(this);

    residhi->setAdjusterListener(this);

    residhithr->setAdjusterListener(this);

    sensilc->setAdjusterListener(this);

    clariFrame->set_label_align(0.025, 0.5);

    clarilres->setAdjusterListener(this);

    claricres->setAdjusterListener(this);

    clarisoft->setLogScale(10, 0);
    clarisoft->setAdjusterListener(this);

    origlcConn = origlc->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::origlcChanged));

    Gtk::HBox* const LCTitleHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const LCLabel = Gtk::manage(new Gtk::Label());
    LCLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYR")) + Glib::ustring("</b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYRLAB")));
    LCLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    LCTitleHBox->pack_start(*LCLabel, Gtk::PACK_EXPAND_WIDGET, 0);
    expcontrastpyr->setLabel(LCTitleHBox);
    setExpandAlignProperties(expcontrastpyr, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    wavgradlConn = wavgradl->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavgradlChanged));

    sigmalc2->setAdjusterListener(this);

    strwav->setAdjusterListener(this);

    angwav->setAdjusterListener(this);

    wavedgConn = wavedg->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavedgChanged));

    strengthw->setAdjusterListener(this);

    sigmaed->setAdjusterListener(this);

    LocalcurveEditorwavedg->setCurveListener(this);

    wavshapeedg->setIdentityValue(0.);
    wavshapeedg->setResetCurve(FlatCurveType(defSpot.locedgwavcurve.at(0)), defSpot.locedgwavcurve);
    wavshapeedg->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    LocalcurveEditorwavedg->curveListComplete();

    gradw->setAdjusterListener(this);

    waveshowConn = waveshow->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::waveshowChanged));

    radiusw->setAdjusterListener(this);

    detailw->setAdjusterListener(this);

    localedgMethod->append(M("TP_WAVELET_RE1"));
    localedgMethod->append(M("TP_WAVELET_RE2"));
    localedgMethod->append(M("TP_WAVELET_RE3"));
    localedgMethod->set_active(0);
    localedgMethodConn = localedgMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabContrast::localedgMethodChanged));

    tloww->setAdjusterListener(this);

    thigw->setAdjusterListener(this);

    edgw->setAdjusterListener(this);

    basew->setAdjusterListener(this);

    localneiMethod->append(M("TP_WAVELET_NPNONE"));
    localneiMethod->append(M("TP_WAVELET_NPLOW"));
    localneiMethod->append(M("TP_WAVELET_NPHIGH"));
    localneiMethod->set_active(0);
    localneiMethodConn = localneiMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabContrast::localneiMethodChanged));

    wavblurConn = wavblur->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavblurChanged));

    levelblur->setAdjusterListener(this);

    sigmabl->setAdjusterListener(this);

    chromablu->setAdjusterListener(this);

    LocalcurveEditorwavlev->setCurveListener(this);

    wavshapelev->setIdentityValue(0.);
    wavshapelev->setResetCurve(FlatCurveType(defSpot.loclevwavcurve.at(0)), defSpot.loclevwavcurve);

    LocalcurveEditorwavlev->curveListComplete();

    residblur->setAdjusterListener(this);

    blurlcConn = blurlc->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::blurlcChanged));

    Gtk::HBox* const LCTitleHBox2 = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const LCLabel2 = Gtk::manage(new Gtk::Label());
    LCLabel2->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYR2")) + Glib::ustring("</b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYR2LAB")));
    LCLabel2->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    LCTitleHBox2->pack_start(*LCLabel2, Gtk::PACK_EXPAND_WIDGET, 0);
    expcontrastpyr2->setLabel(LCTitleHBox2);
    setExpandAlignProperties(expcontrastpyr2, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    wavcontConn = wavcont->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavcontChanged));

    sigma->setAdjusterListener(this);

    offset->setAdjusterListener(this);

    chromalev->setAdjusterListener(this);

    LocalcurveEditorwavcon->setCurveListener(this);

    wavshapecon->setIdentityValue(0.);
    wavshapecon->setResetCurve(FlatCurveType(defSpot.locconwavcurve.at(0)), defSpot.locconwavcurve);

    LocalcurveEditorwavcon->curveListComplete();

    wavcompreConn = wavcompre->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavcompreChanged));

    LocalcurveEditorwavcompre->setCurveListener(this);

    wavshapecompre->setIdentityValue(0.);
    wavshapecompre->setResetCurve(FlatCurveType(defSpot.loccomprewavcurve.at(0)), defSpot.loccomprewavcurve);
    wavshapecompre->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));

    LocalcurveEditorwavcompre->curveListComplete();

    sigmadr->setAdjusterListener(this);

    threswav->setAdjusterListener(this);

    residcomp->setAdjusterListener(this);

    wavcompConn = wavcomp->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavcompChanged));

    sigmadc->setAdjusterListener(this);

    deltad->setAdjusterListener(this);

    LocalcurveEditorwavcomp->setCurveListener(this);

    wavshapecomp->setIdentityValue(0.);
    wavshapecomp->setResetCurve(FlatCurveType(defSpot.loccompwavcurve.at(0)), defSpot.loccompwavcurve);
    wavshapecomp->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    LocalcurveEditorwavcomp->curveListComplete();

    fatres->setAdjusterListener(this);

    fftwlcConn = fftwlc->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::fftwlcChanged));

    setExpandAlignProperties(expmasklc, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmasklcMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasklcMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmasklcMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmasklcMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmasklcMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmasklcMethod->set_active(0);
    showmasklcMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmasklcMethodConn = showmasklcMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabContrast::showmasklcMethodChanged));

    enalcMaskConn = enalcMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::enalcMaskChanged));

    masklcCurveEditorG->setCurveListener(this);

    CCmasklcshape->setIdentityValue(0.);
    CCmasklcshape->setResetCurve(FlatCurveType(defSpot.CCmasklccurve.at(0)), defSpot.CCmasklccurve);
    CCmasklcshape->setBottomBarColorProvider(this, 1);

    LLmasklcshape->setIdentityValue(0.);
    LLmasklcshape->setResetCurve(FlatCurveType(defSpot.LLmasklccurve.at(0)), defSpot.LLmasklccurve);
    LLmasklcshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmasklcshape->setIdentityValue(0.);
    HHmasklcshape->setResetCurve(FlatCurveType(defSpot.HHmasklccurve.at(0)), defSpot.HHmasklccurve);
    HHmasklcshape->setCurveColorProvider(this, 2);
    HHmasklcshape->setBottomBarColorProvider(this, 2);

    masklcCurveEditorG->curveListComplete();

    blendmasklc->setAdjusterListener(this);

    radmasklc->setAdjusterListener(this);

    chromasklc->setAdjusterListener(this);

    mask2lcCurveEditorG->setCurveListener(this);

    Lmasklcshape->setResetCurve(DiagonalCurveType(defSpot.Lmasklccurve.at(0)), defSpot.Lmasklccurve);
    Lmasklcshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmasklcshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2lcCurveEditorG->curveListComplete();

    // Add Local contrast specific widgets to GUI
    pack_start(*localcontMethod);
    pack_start(*lcradius);
    pack_start(*lcamount);
    pack_start(*lcdarkness);
    pack_start(*lclightness);
    ToolParamBlock* const coBox = Gtk::manage(new ToolParamBlock());
    coBox->pack_start(*sigmalc);
    coBox->pack_start(*LocalcurveEditorwav, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    coBox->pack_start(*csThreshold);
    contFrame->add(*coBox);
    pack_start(*contFrame);
    // pack_start(*levelwav);
    ToolParamBlock* const resiBox = Gtk::manage(new ToolParamBlock());
    resiBox->pack_start(*residcont);
    resiBox->pack_start(*residchro);
    Gtk::Frame* const shresFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHRESFRA")));
    shresFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const shresBox = Gtk::manage(new ToolParamBlock());
    shresBox->pack_start(*residsha);
    shresBox->pack_start(*residshathr);
    shresBox->pack_start(*residhi);
    shresBox->pack_start(*residhithr);
    shresFrame->add(*shresBox);
    resiBox->pack_start(*shresFrame);
    expresidpyr->add(*resiBox, false);
    pack_start(*expresidpyr);
    pack_start(*sensilc);
    Gtk::HSeparator* const separatorcontr = Gtk::manage(new  Gtk::HSeparator());
    pack_start(*separatorcontr);
    ToolParamBlock* const clariBox = Gtk::manage(new ToolParamBlock());
    clariBox->pack_start(*clarilres);
    clariBox->pack_start(*claricres);
    clariBox->pack_start(*clarisoft);
    clariBox->pack_start(*origlc);
    clariFrame->add(*clariBox);
    pack_start(*clariFrame);
    ToolParamBlock* const blurcontBox = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const gradwavFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADWAVFRA")));
    gradwavFrame->set_label_align(0.025, 0.5);
    gradwavFrame->set_label_widget(*wavgradl);
    ToolParamBlock* const gradwavBox = Gtk::manage(new ToolParamBlock());
    gradwavBox->pack_start(*sigmalc2);
    gradwavBox->pack_start(*strwav);
    gradwavBox->pack_start(*angwav);
    gradwavFrame->add(*gradwavBox);
    blurcontBox->pack_start(*gradwavFrame);
    Gtk::Frame* const edgFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_EDGSHARPFRA")));
    edgFrame->set_label_align(0.025, 0.5);
    edgFrame->set_label_widget(*wavedg);
    ToolParamBlock* const edgsBox = Gtk::manage(new ToolParamBlock());
    edgsBox->pack_start(*strengthw);
    edgsBox->pack_start(*sigmaed);
    edgsBox->pack_start(*LocalcurveEditorwavedg, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    edgsBox->pack_start(*gradw);
    edgsBox->pack_start(*waveshow);
    edgsBoxshow->pack_start(*radiusw);
    edgsBoxshow->pack_start(*detailw);
    Gtk::HBox* const edbox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const labmedgr = Gtk::manage(new Gtk::Label(M("TP_WAVELET_MEDGREINF") + ":"));
    edbox->pack_start(*labmedgr, Gtk::PACK_SHRINK, 1);
    edbox->pack_start(*localedgMethod);
    edgsBoxshow->pack_start(*edbox);
    Gtk::HSeparator* const separatoredg2 = Gtk::manage(new  Gtk::HSeparator());
    edgsBoxshow->pack_start(*separatoredg2);
    edgsBoxshow->pack_start(*tloww);
    edgsBoxshow->pack_start(*thigw);
    Gtk::HSeparator* const separatoredg = Gtk::manage(new  Gtk::HSeparator());
    edgsBoxshow->pack_start(*separatoredg);
    edgsBoxshow->pack_start(*edgw);
    edgsBoxshow->pack_start(*basew);
    Gtk::HBox* const ctboxNP = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const labmNP = Gtk::manage(new Gtk::Label(M("TP_WAVELET_NPTYPE") + ":"));
    ctboxNP->pack_start(*labmNP, Gtk::PACK_SHRINK, 1);
    ctboxNP->pack_start(*localneiMethod);
    edgsBoxshow->pack_start(*ctboxNP);
    edgsBox->pack_start(*edgsBoxshow);
    edgFrame->add(*edgsBox);
    blurcontBox->pack_start(*edgFrame);
    Gtk::Frame* const blurlevelFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_BLURLEVELFRA")));
    blurlevelFrame->set_label_align(0.025, 0.5);
    blurlevelFrame->set_label_widget(*wavblur);
    Gtk::VBox* const blurlevcontBox = Gtk::manage(new Gtk::VBox());
    blurlevcontBox->set_spacing(2);
    blurlevcontBox->pack_start(*levelblur);
    blurlevcontBox->pack_start(*sigmabl);
    blurlevcontBox->pack_start(*chromablu);
    blurlevcontBox->pack_start(*LocalcurveEditorwavlev, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    Gtk::HSeparator* const separatorblu = Gtk::manage(new  Gtk::HSeparator());
    blurlevcontBox->pack_start(*separatorblu);
    blurlevcontBox->pack_start(*residblur);
    blurlevcontBox->pack_start(*blurlc);
    blurlevelFrame->add(*blurlevcontBox);
    blurcontBox->pack_start(*blurlevelFrame);
    expcontrastpyr->add(*blurcontBox, false);
    pack_start(*expcontrastpyr);
    ToolParamBlock* const blurcontBox2 = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const contFrame2 = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CONTFRA")));
    contFrame2->set_label_align(0.025, 0.5);
    contFrame2->set_label_widget(*wavcont);
    Gtk::VBox* const contlevBox = Gtk::manage(new Gtk::VBox());
    contlevBox->set_spacing(2);
    contlevBox->pack_start(*sigma);
    contlevBox->pack_start(*offset);
    contlevBox->pack_start(*chromalev);
    contlevBox->pack_start(*LocalcurveEditorwavcon, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    contFrame2->add(*contlevBox);
    blurcontBox2->pack_start(*contFrame2);
    Gtk::Frame* const compreFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_COMPREFRA")));
    compreFrame->set_label_align(0.025, 0.5);
    compreFrame->set_label_widget(*wavcompre);
    Gtk::VBox* const compreBox = Gtk::manage(new Gtk::VBox());
    compreBox->set_spacing(2);
    compreBox->pack_start(*LocalcurveEditorwavcompre, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    compreBox->pack_start(*sigmadr);
    compreBox->pack_start(*threswav);
    compreBox->pack_start(*residcomp);
    compreFrame->add(*compreBox);
    blurcontBox2->pack_start(*compreFrame);
    Gtk::Frame* const compFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_COMPFRA")));
    compFrame->set_label_align(0.025, 0.5);
    compFrame->set_label_widget(*wavcomp);
    Gtk::VBox* const compBox = Gtk::manage(new Gtk::VBox());
    compBox->set_spacing(2);
    compBox->pack_start(*sigmadc);
    compBox->pack_start(*deltad);
    compBox->pack_start(*LocalcurveEditorwavcomp, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    // Gtk::HSeparator* const separatorcomp = Gtk::manage(new  Gtk::HSeparator());
    // compBox->pack_start(*separatorcomp);
    // compBox->pack_start(*fatres);
    compFrame->add(*compBox);
    blurcontBox2->pack_start(*compFrame);
    expcontrastpyr2->add(*blurcontBox2, false);
    pack_start(*expcontrastpyr2);
    pack_start(*fftwlc);
    ToolParamBlock* const masklcBox = Gtk::manage(new ToolParamBlock());
    masklcBox->pack_start(*showmasklcMethod, Gtk::PACK_SHRINK, 4);
    masklcBox->pack_start(*enalcMask, Gtk::PACK_SHRINK, 0);
    masklcBox->pack_start(*masklcCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    masklcBox->pack_start(*blendmasklc, Gtk::PACK_SHRINK, 0);
    masklcBox->pack_start(*radmasklc, Gtk::PACK_SHRINK, 0);
    masklcBox->pack_start(*chromasklc, Gtk::PACK_SHRINK, 0);
    masklcBox->pack_start(*mask2lcCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    expmasklc->add(*masklcBox, false);
    pack_start(*expmasklc, false, false);
}

LocallabContrast::~LocallabContrast()
{
    delete LocalcurveEditorwav;
    delete LocalcurveEditorwavedg;
    delete LocalcurveEditorwavlev;
    delete LocalcurveEditorwavcon;
    delete LocalcurveEditorwavcompre;
    delete LocalcurveEditorwavcomp;
    delete masklcCurveEditorG;
    delete mask2lcCurveEditorG;
}

bool LocallabContrast::isMaskViewActive()
{
    return (showmasklcMethod->get_active_row_number() != 0);
}

void LocallabContrast::resetMaskView()
{
    showmasklcMethodConn.block(true);
    showmasklcMethod->set_active(0);
    showmasklcMethodConn.block(false);
}

void LocallabContrast::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &maskMask)
{
    lcMask = showmasklcMethod->get_active_row_number();
}

void LocallabContrast::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        contFrame->set_tooltip_text(M("TP_LOCALLAB_EXPCONTRAST_TOOLTIP"));
        LocalcurveEditorwav->set_tooltip_markup(M("TP_LOCALLAB_LEVELLOCCONTRAST_TOOLTIP"));
        levelwav->set_tooltip_markup(M("TP_LOCALLAB_LEVELWAV_TOOLTIP"));
        clariFrame->set_tooltip_markup(M("TP_LOCALLAB_CLARI_TOOLTIP"));
        clarisoft->set_tooltip_markup(M("TP_LOCALLAB_CLARISOFT_TOOLTIP"));
        expcontrastpyr->set_tooltip_text(M("TP_LOCALLAB_EXPCONTRASTPYR_TOOLTIP"));
        wavgradl->set_tooltip_text(M("TP_LOCALLAB_WAVGRAD_TOOLTIP"));
        wavedg->set_tooltip_text(M("TP_LOCALLAB_WAVEEDG_TOOLTIP"));
        wavblur->set_tooltip_text(M("TP_LOCALLAB_WAVBLUR_TOOLTIP"));
        chromablu->set_tooltip_text(M("TP_LOCALLAB_CHROMABLU_TOOLTIP"));
        expcontrastpyr2->set_tooltip_text(M("TP_LOCALLAB_EXPCONTRASTPYR_TOOLTIP"));
        wavcont->set_tooltip_text(M("TP_LOCALLAB_WAVCONTF_TOOLTIP"));
        chromalev->set_tooltip_text(M("TP_LOCALLAB_CHROMABLU_TOOLTIP"));
        wavcompre->set_tooltip_text(M("TP_LOCALLAB_WAVCOMPRE_TOOLTIP"));
        wavcomp->set_tooltip_text(M("TP_LOCALLAB_WAVCOMP_TOOLTIP"));
        fftwlc->set_tooltip_text(M("TP_LOCALLAB_LC_FFTW_TOOLTIP"));
        expmasklc->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmasklcshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmasklcshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmasklcshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        blendmasklc->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        mask2lcCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        Lmasklcshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    } else {
        contFrame->set_tooltip_text("");
        LocalcurveEditorwav->set_tooltip_text("");
        levelwav->set_tooltip_text("");
        clariFrame->set_tooltip_text("");
        clarisoft->set_tooltip_text("");
        expcontrastpyr->set_tooltip_text("");
        wavgradl->set_tooltip_text("");
        wavedg->set_tooltip_text("");
        wavblur->set_tooltip_text("");
        chromablu->set_tooltip_text("");
        expcontrastpyr2->set_tooltip_text("");
        wavcont->set_tooltip_text("");
        chromalev->set_tooltip_text("");
        wavcompre->set_tooltip_text("");
        wavcomp->set_tooltip_text("");
        fftwlc->set_tooltip_text("");
        expmasklc->set_tooltip_text("");
        CCmasklcshape->setTooltip("");
        LLmasklcshape->setTooltip("");
        HHmasklcshape->setTooltip("");
        blendmasklc->set_tooltip_text("");
        mask2lcCurveEditorG->set_tooltip_text("");
        Lmasklcshape->setTooltip("");
    }
}

void LocallabContrast::setDefaultExpanderVisibility()
{
    expresidpyr->set_expanded(false);
    expcontrastpyr->set_expanded(false);
    expcontrastpyr2->set_expanded(false);
    expmasklc->set_expanded(false);
}

void LocallabContrast::disableListener()
{
    LocallabTool::disableListener();

    localcontMethodConn.block(true);
    origlcConn.block(true);
    wavgradlConn.block(true);
    wavedgConn.block(true);
    localedgMethodConn.block(true);
    waveshowConn.block(true);
    localneiMethodConn.block(true);
    wavblurConn.block(true);
    blurlcConn.block(true);
    wavcontConn.block(true);
    wavcompreConn.block(true);
    wavcompConn.block(true);
    fftwlcConn.block(true);
    showmasklcMethodConn.block(true);
    enalcMaskConn.block(true);
}

void LocallabContrast::enableListener()
{
    LocallabTool::enableListener();

    localcontMethodConn.block(false);
    origlcConn.block(false);
    wavgradlConn.block(false);
    wavedgConn.block(false);
    localedgMethodConn.block(false);
    waveshowConn.block(false);
    localneiMethodConn.block(false);
    wavblurConn.block(false);
    blurlcConn.block(false);
    wavcontConn.block(false);
    wavcompreConn.block(false);
    wavcompConn.block(false);
    fftwlcConn.block(false);
    showmasklcMethodConn.block(false);
    enalcMaskConn.block(false);
}

void LocallabContrast::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visicontrast);
        exp->setEnabled(spot.expcontrast);
        complexity->set_active(spot.complexcontrast);

        if (spot.localcontMethod == "loc") {
            localcontMethod->set_active(0);
        } else if (spot.localcontMethod == "wav") {
            localcontMethod->set_active(1);
        }

        fftwlc->set_active(spot.fftwlc);
        // Update Local contrast GUI according to fftwlc button state
        // Note: Contrary to the others, shall be called before setting lcradius value
        updateContrastGUI3();
        lcradius->setValue((double)spot.lcradius);
        lcamount->setValue(spot.lcamount);
        lcdarkness->setValue(spot.lcdarkness);
        lclightness->setValue(spot.lclightness);
        sigmalc->setValue(spot.sigmalc);
        wavshape->setCurve(spot.locwavcurve);
        csThreshold->setValue<int>(spot.csthreshold);
        levelwav->setValue((double)spot.levelwav);
        residcont->setValue(spot.residcont);
        residchro->setValue(spot.residchro);
        residsha->setValue(spot.residsha);
        residshathr->setValue(spot.residshathr);
        residhi->setValue(spot.residhi);
        residhithr->setValue(spot.residhithr);
        sensilc->setValue((double)spot.sensilc);
        clarilres->setValue(spot.clarilres);
        claricres->setValue(spot.claricres);
        clarisoft->setValue(spot.clarisoft);
        origlc->set_active(spot.origlc);
        wavgradl->set_active(spot.wavgradl);
        sigmalc2->setValue(spot.sigmalc2);
        strwav->setValue(spot.strwav);
        angwav->setValue(spot.angwav);
        wavedg->set_active(spot.wavedg);
        strengthw->setValue(spot.strengthw);
        sigmaed->setValue(spot.sigmaed);
        wavshapeedg->setCurve(spot.locedgwavcurve);
        gradw->setValue(spot.gradw);
        waveshow->set_active(spot.waveshow);
        radiusw->setValue(spot.radiusw);
        detailw->setValue(spot.detailw);

        if (spot.localedgMethod == "fir") {
            localedgMethod->set_active(0);
        } else if (spot.localedgMethod == "sec") {
            localedgMethod->set_active(1);
        } else if (spot.localedgMethod == "thr") {
            localedgMethod->set_active(2);
        }

        tloww->setValue(spot.tloww);
        thigw->setValue(spot.thigw);
        edgw->setValue(spot.edgw);
        basew->setValue(spot.basew);

        if (spot.localneiMethod == "none") {
            localneiMethod->set_active(0);
        } else if (spot.localneiMethod == "low") {
            localneiMethod->set_active(1);
        } else if (spot.localneiMethod == "high") {
            localneiMethod->set_active(2);
        }

        wavblur->set_active(spot.wavblur);
        levelblur->setValue(spot.levelblur);
        sigmabl->setValue(spot.sigmabl);
        chromablu->setValue(spot.chromablu);
        wavshapelev->setCurve(spot.loclevwavcurve);
        residblur->setValue(spot.residblur);
        blurlc->set_active(spot.blurlc);
        wavcont->set_active(spot.wavcont);
        sigma->setValue(spot.sigma);
        offset->setValue(spot.offset);
        chromalev->setValue(spot.chromalev);
        wavshapecon->setCurve(spot.locconwavcurve);
        wavcompre->set_active(spot.wavcompre);
        wavshapecompre->setCurve(spot.loccomprewavcurve);
        sigmadr->setValue(spot.sigmadr);
        threswav->setValue(spot.threswav);
        residcomp->setValue(spot.residcomp);
        wavcomp->set_active(spot.wavcomp);
        sigmadc->setValue(spot.sigmadc);
        deltad->setValue(spot.deltad);
        wavshapecomp->setCurve(spot.loccompwavcurve);
        fatres->setValue(spot.fatres);
        enalcMask->set_active(spot.enalcMask);
        CCmasklcshape->setCurve(spot.CCmasklccurve);
        LLmasklcshape->setCurve(spot.LLmasklccurve);
        HHmasklcshape->setCurve(spot.HHmasklccurve);
        blendmasklc->setValue((double)spot.blendmasklc);
        radmasklc->setValue(spot.radmasklc);
        chromasklc->setValue(spot.chromasklc);
        Lmasklcshape->setCurve(spot.Lmasklccurve);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Update Local contrast GUI according to localcontMethod combobox value
    updateContrastGUI1();

    // Update Local contrast GUI according to waveshow button state
    updateContrastGUI2();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabContrast::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expcontrast = exp->getEnabled();
        spot.visicontrast = exp->get_visible();
        spot.complexcontrast = complexity->get_active_row_number();

        if (localcontMethod->get_active_row_number() == 0) {
            spot.localcontMethod = "loc";
        } else if (localcontMethod->get_active_row_number() == 1) {
            spot.localcontMethod = "wav";
        }

        spot.lcradius = lcradius->getIntValue();
        spot.lcamount = lcamount->getValue();
        spot.lcdarkness = lcdarkness->getValue();
        spot.lclightness = lclightness->getValue();
        spot.sigmalc = sigmalc->getValue();
        spot.locwavcurve = wavshape->getCurve();
        spot.csthreshold = csThreshold->getValue<int>();
        spot.levelwav = levelwav->getIntValue();
        spot.residcont = residcont->getValue();
        spot.residchro = residchro->getValue();
        spot.residsha = residsha->getValue();
        spot.residshathr = residshathr->getValue();
        spot.residhi = residhi->getValue();
        spot.residhithr = residhithr->getValue();
        spot.sensilc = sensilc->getIntValue();
        spot.clarilres = clarilres->getValue();
        spot.claricres = claricres->getValue();
        spot.clarisoft = clarisoft->getValue();
        spot.origlc = origlc->get_active();
        spot.wavgradl = wavgradl->get_active();
        spot.sigmalc2 = sigmalc2->getValue();
        spot.strwav = strwav->getValue();
        spot.angwav = angwav->getValue();
        spot.wavedg = wavedg->get_active();
        spot.strengthw = strengthw->getValue();
        spot.sigmaed = sigmaed->getValue();
        spot.locedgwavcurve = wavshapeedg->getCurve();
        spot.gradw = gradw->getValue();
        spot.waveshow = waveshow->get_active();
        spot.radiusw = radiusw->getValue();
        spot.detailw = detailw->getValue();

        if (localedgMethod->get_active_row_number() == 0) {
            spot.localedgMethod = "fir";
        } else if (localedgMethod->get_active_row_number() == 1) {
            spot.localedgMethod = "sec";
        } else if (localedgMethod->get_active_row_number() == 2) {
            spot.localedgMethod = "thr";
        }

        spot.tloww = tloww->getValue();
        spot.thigw = thigw->getValue();
        spot.edgw = edgw->getValue();
        spot.basew = basew->getValue();

        if (localneiMethod->get_active_row_number() == 0) {
            spot.localneiMethod = "none";
        } else if (localneiMethod->get_active_row_number() == 1) {
            spot.localneiMethod = "low";
        } else if (localneiMethod->get_active_row_number() == 2) {
            spot.localneiMethod = "high";
        }

        spot.wavblur = wavblur->get_active();
        spot.levelblur = levelblur->getValue();
        spot.sigmabl = sigmabl->getValue();
        spot.chromablu = chromablu->getValue();
        spot.loclevwavcurve = wavshapelev->getCurve();
        spot.residblur = residblur->getValue();
        spot.blurlc = blurlc->get_active();
        spot.wavcont = wavcont->get_active();
        spot.sigma = sigma->getValue();
        spot.offset = offset->getValue();
        spot.chromalev = chromalev->getValue();
        spot.locconwavcurve = wavshapecon->getCurve();
        spot.wavcompre = wavcompre->get_active();
        spot.loccomprewavcurve = wavshapecompre->getCurve();
        spot.sigmadr = sigmadr->getValue();
        spot.threswav = threswav->getValue();
        spot.residcomp = residcomp->getValue();
        spot.wavcomp = wavcomp->get_active();
        spot.sigmadc = sigmadc->getValue();
        spot.deltad = deltad->getValue();
        spot.loccompwavcurve = wavshapecomp->getCurve();
        spot.fatres = fatres->getValue();
        spot.fftwlc = fftwlc->get_active();
        spot.enalcMask = enalcMask->get_active();
        spot.CCmasklccurve = CCmasklcshape->getCurve();
        spot.LLmasklccurve = LLmasklcshape->getCurve();
        spot.HHmasklccurve = HHmasklcshape->getCurve();
        spot.blendmasklc = blendmasklc->getIntValue();
        spot.radmasklc = radmasklc->getValue();
        spot.chromasklc = chromasklc->getValue();
        spot.Lmasklccurve = Lmasklcshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabContrast::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster and threshold adjuster widgets
        lcradius->setDefault((double)defSpot.lcradius);
        lcamount->setDefault(defSpot.lcamount);
        lcdarkness->setDefault(defSpot.lcdarkness);
        lclightness->setDefault(defSpot.lclightness);
        sigmalc->setDefault(defSpot.sigmalc);
        levelwav->setDefault((double)defSpot.levelwav);
        csThreshold->setDefault<int>(defSpot.csthreshold);
        residcont->setDefault(defSpot.residcont);
        residchro->setDefault(defSpot.residchro);
        residsha->setDefault(defSpot.residsha);
        residshathr->setDefault(defSpot.residshathr);
        residhi->setDefault(defSpot.residhi);
        residhithr->setDefault(defSpot.residhithr);
        sensilc->setDefault((double)defSpot.sensilc);
        clarilres->setDefault(defSpot.clarilres);
        claricres->setDefault(defSpot.claricres);
        clarisoft->setDefault(defSpot.clarisoft);
        sigmalc2->setDefault(defSpot.sigmalc2);
        strwav->setDefault(defSpot.strwav);
        angwav->setDefault(defSpot.angwav);
        strengthw->setDefault(defSpot.strengthw);
        sigmaed->setDefault(defSpot.sigmaed);
        gradw->setDefault(defSpot.gradw);
        radiusw->setDefault(defSpot.radiusw);
        detailw->setDefault(defSpot.detailw);
        tloww->setDefault(defSpot.tloww);
        thigw->setDefault(defSpot.thigw);
        edgw->setDefault(defSpot.edgw);
        basew->setDefault(defSpot.basew);
        levelblur->setDefault(defSpot.levelblur);
        sigmabl->setDefault(defSpot.sigmabl);
        chromablu->setDefault(defSpot.chromablu);
        residblur->setDefault(defSpot.residblur);
        sigma->setDefault(defSpot.sigma);
        offset->setDefault(defSpot.offset);
        chromalev->setDefault(defSpot.chromalev);
        sigmadr->setDefault(defSpot.sigmadr);
        threswav->setDefault(defSpot.threswav);
        residcomp->setDefault(defSpot.residcomp);
        sigmadc->setDefault(defSpot.sigmadc);
        deltad->setDefault(defSpot.deltad);
        fatres->setDefault(defSpot.fatres);
        blendmasklc->setDefault((double)defSpot.blendmasklc);
        radmasklc->setDefault(defSpot.radmasklc);
        chromasklc->setDefault(defSpot.chromasklc);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabContrast::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == lcradius) {
            if (listener) {
                listener->panelChanged(Evlocallablcradius,
                                       lcradius->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lcamount) {
            if (listener) {
                listener->panelChanged(Evlocallablcamount,
                                       lcamount->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lcdarkness) {
            if (listener) {
                listener->panelChanged(Evlocallablcdarkness,
                                       lcdarkness->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lclightness) {
            if (listener) {
                listener->panelChanged(Evlocallablclightness,
                                       lclightness->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sigmalc) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmalc,
                                       sigmalc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == levelwav) {
            if (listener) {
                listener->panelChanged(Evlocallablevelwav,
                                       levelwav->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == residcont) {
            if (listener) {
                listener->panelChanged(Evlocallabresidcont,
                                       residcont->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == residchro) {
            if (listener) {
                listener->panelChanged(Evlocallabresidchro,
                                       residchro->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == residsha) {
            if (listener) {
                listener->panelChanged(Evlocallabresidsha,
                                       residsha->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == residshathr) {
            if (listener) {
                listener->panelChanged(Evlocallabresidshathr,
                                       residshathr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == residhi) {
            if (listener) {
                listener->panelChanged(Evlocallabresidhi,
                                       residhi->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == residhithr) {
            if (listener) {
                listener->panelChanged(Evlocallabresidhithr,
                                       residhithr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensilc) {
            if (listener) {
                listener->panelChanged(Evlocallabsensilc,
                                       sensilc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == clarilres) {
            if (listener) {
                listener->panelChanged(Evlocallabclarilres,
                                       clarilres->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == claricres) {
            if (listener) {
                listener->panelChanged(Evlocallabclaricres,
                                       claricres->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == clarisoft) {
            if (listener) {
                listener->panelChanged(Evlocallabclarisoft,
                                       clarisoft->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sigmalc2) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmalc2,
                                       sigmalc2->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strwav) {
            if (listener) {
                listener->panelChanged(Evlocallabstrwav,
                                       strwav->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == angwav) {
            if (listener) {
                listener->panelChanged(Evlocallabangwav,
                                       angwav->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strengthw) {
            if (listener) {
                listener->panelChanged(Evlocallabstrengthw,
                                       strengthw->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sigmaed) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmaed,
                                       sigmaed->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gradw) {
            if (listener) {
                listener->panelChanged(Evlocallabgradw,
                                       gradw->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radiusw) {
            if (listener) {
                listener->panelChanged(Evlocallabradiusw,
                                       radiusw->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == detailw) {
            if (listener) {
                listener->panelChanged(Evlocallabdetailw,
                                       detailw->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == tloww) {
            if (listener) {
                listener->panelChanged(Evlocallabtloww,
                                       tloww->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == thigw) {
            if (listener) {
                listener->panelChanged(Evlocallabthigw,
                                       thigw->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == edgw) {
            if (listener) {
                listener->panelChanged(Evlocallabedgw,
                                       edgw->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == basew) {
            if (listener) {
                listener->panelChanged(Evlocallabbasew,
                                       basew->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == levelblur) {
            if (listener) {
                listener->panelChanged(Evlocallablevelblur,
                                       levelblur->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sigmabl) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmabl,
                                       sigmabl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromablu) {
            if (listener) {
                listener->panelChanged(Evlocallabchromablu,
                                       chromablu->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == residblur) {
            if (listener) {
                listener->panelChanged(Evlocallabresidblur,
                                       residblur->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sigma) {
            if (listener) {
                listener->panelChanged(Evlocallabsigma,
                                       sigma->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == offset) {
            if (listener) {
                listener->panelChanged(Evlocallaboffset,
                                       offset->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromalev) {
            if (listener) {
                listener->panelChanged(Evlocallabchromalev,
                                       chromalev->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sigmadr) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmadr,
                                       sigmadr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }


        if (a == threswav) {
            if (listener) {
                listener->panelChanged(Evlocallabthreswav,
                                       threswav->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == residcomp) {
            if (listener) {
                listener->panelChanged(Evlocallabresidcomp,
                                       residcomp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sigmadc) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmadc,
                                       sigmadc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == deltad) {
            if (listener) {
                listener->panelChanged(Evlocallabdeltad,
                                       deltad->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == fatres) {
            if (listener) {
                listener->panelChanged(Evlocallabfatres,
                                       fatres->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmasklc) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmasklc,
                                       blendmasklc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmasklc) {
            if (listener) {
                listener->panelChanged(Evlocallabradmasklc,
                                       radmasklc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromasklc) {
            if (listener) {
                listener->panelChanged(Evlocallabchromasklc,
                                       chromasklc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabcsThreshold,
                                   csThreshold->getHistoryString() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabContrast::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == wavshape) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurve,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == wavshapeedg) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurveedg,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == wavshapelev) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvelev,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == wavshapecon) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvecon,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == wavshapecompre) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvecompre,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == wavshapecomp) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvecomp,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == CCmasklcshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmasklcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmasklcshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmasklcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmasklcshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmasklcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmasklcshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmasklcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenacontrast,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenacontrast,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    origlc->set_active(defSpot.origlc);
    wavgradl->set_active(defSpot.wavgradl);
    sigmalc2->setValue(defSpot.sigmalc2);
    strwav->setValue(defSpot.strwav);
    angwav->setValue(defSpot.angwav);
    wavedg->set_active(defSpot.wavedg);
    strengthw->setValue(defSpot.strengthw);
    sigmaed->setValue(defSpot.sigmaed);
    wavshapeedg->setCurve(defSpot.locedgwavcurve);
    gradw->setValue(defSpot.gradw);
    waveshow->set_active(defSpot.waveshow);
    radiusw->setValue(defSpot.radiusw);
    detailw->setValue(defSpot.detailw);

    if (defSpot.localedgMethod == "fir") {
        localedgMethod->set_active(0);
    } else if (defSpot.localedgMethod == "sec") {
        localedgMethod->set_active(1);
    } else if (defSpot.localedgMethod == "thr") {
        localedgMethod->set_active(2);
    }

    tloww->setValue(defSpot.tloww);
    thigw->setValue(defSpot.thigw);
    edgw->setValue(defSpot.edgw);
    basew->setValue(defSpot.basew);

    if (defSpot.localneiMethod == "none") {
        localneiMethod->set_active(0);
    } else if (defSpot.localneiMethod == "low") {
        localneiMethod->set_active(1);
    } else if (defSpot.localneiMethod == "high") {
        localneiMethod->set_active(2);
    }

    wavblur->set_active(defSpot.wavblur);
    levelblur->setValue(defSpot.levelblur);
    sigmabl->setValue(defSpot.sigmabl);
    chromablu->setValue(defSpot.chromablu);
    wavshapelev->setCurve(defSpot.loclevwavcurve);
    residblur->setValue(defSpot.residblur);
    blurlc->set_active(defSpot.blurlc);
    wavcont->set_active(defSpot.wavcont);
    sigma->setValue(defSpot.sigma);
    offset->setValue(defSpot.offset);
    chromalev->setValue(defSpot.chromalev);
    wavshapecon->setCurve(defSpot.locconwavcurve);
    wavcompre->set_active(defSpot.wavcompre);
    wavshapecompre->setCurve(defSpot.loccomprewavcurve);
    sigmadr->setValue(defSpot.sigmadr);
    threswav->setValue(defSpot.threswav);
    residcomp->setValue(defSpot.residcomp);
    wavcomp->set_active(defSpot.wavcomp);
    sigmadc->setValue(defSpot.sigmadc);
    deltad->setValue(defSpot.deltad);
    wavshapecomp->setCurve(defSpot.loccompwavcurve);
    fatres->setValue(defSpot.fatres);
    fftwlc->set_active(defSpot.fftwlc);

    // Enable all listeners
    enableListener();

    // Update GUI based on converted widget parameters:
    // - Update Local contrast GUI according to fftwlc button state
    updateContrastGUI3();
}

void LocallabContrast::convertParamToSimple()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden specific GUI widgets in Simple mode to default spot values
    if (defSpot.localcontMethod == "loc") {
        localcontMethod->set_active(0);
    } else if (defSpot.localcontMethod == "wav") {
        localcontMethod->set_active(1);
    }

    showmasklcMethod->set_active(0);
    enalcMask->set_active(defSpot.enalcMask);
    CCmasklcshape->setCurve(defSpot.CCmasklccurve);
    LLmasklcshape->setCurve(defSpot.LLmasklccurve);
    HHmasklcshape->setCurve(defSpot.HHmasklccurve);
    blendmasklc->setValue((double)defSpot.blendmasklc);
    radmasklc->setValue(defSpot.radmasklc);
    chromasklc->setValue(defSpot.chromasklc);
    Lmasklcshape->setCurve(defSpot.Lmasklccurve);

    // Enable all listeners
    enableListener();

    // Update GUI based on converted widget parameters:
    // - Update Local contrast GUI according to localcontMethod combobox value
    updateContrastGUI1();
}

void LocallabContrast::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            localcontMethod->hide();
            origlc->hide();
            expcontrastpyr->hide();
            expcontrastpyr2->hide();
            fftwlc->hide();
            expmasklc->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            origlc->hide();
            expcontrastpyr->hide();
            expcontrastpyr2->hide();
            fftwlc->hide();
            // Specific Simple mode widgets are shown in Normal mode
            localcontMethod->show();
            expmasklc->show();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            localcontMethod->show();
            origlc->show();

            if (localcontMethod->get_active_row_number() != 0) { // Keep widgets hidden when localcontMethod is equal to 0
                expcontrastpyr->show();
                expcontrastpyr2->show();
            }

            if (localcontMethod->get_active_row_number() != 1) { // Keep widget hidden when localcontMethod is equal to 1
                fftwlc->show();
            }

            expmasklc->show();
    }
}

void LocallabContrast::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmasklcshape->updateLocallabBackground(normChromar);
        LLmasklcshape->updateLocallabBackground(normLumar);
        HHmasklcshape->updateLocallabBackground(normHuer);
        Lmasklcshape->updateLocallabBackground(normLumar);

        return false;
    }
    );
}

void LocallabContrast::localcontMethodChanged()
{
    // Update Local contrast GUI according to localcontMethod combobox value
    updateContrastGUI1();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallablocalcontMethod,
                                   localcontMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabContrast::origlcChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (origlc->get_active()) {
                listener->panelChanged(Evlocallaboriglc,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallaboriglc,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::wavgradlChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (wavgradl->get_active()) {
                listener->panelChanged(Evlocallabwavgradl,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabwavgradl,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::wavedgChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (wavedg->get_active()) {
                listener->panelChanged(Evlocallabwavedg,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabwavedg,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::localedgMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallablocaledgMethod,
                                   localedgMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabContrast::waveshowChanged()
{
    // Update Local contrast GUI according to waveshow button state
    updateContrastGUI2();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (waveshow->get_active()) {
                listener->panelChanged(Evlocallabwaveshow,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabwaveshow,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::localneiMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallablocalneiMethod,
                                   localneiMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabContrast::wavblurChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (wavblur->get_active()) {
                listener->panelChanged(Evlocallabwavblur,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabwavblur,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::blurlcChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (blurlc->get_active()) {
                listener->panelChanged(Evlocallabblurlc,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabblurlc,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::wavcontChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (wavcont->get_active()) {
                listener->panelChanged(Evlocallabwavcont,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabwavcont,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::wavcompreChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (wavcompre->get_active()) {
                listener->panelChanged(Evlocallabwavcompre,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabwavcompre,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::wavcompChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (wavcomp->get_active()) {
                listener->panelChanged(Evlocallabwavcomp,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabwavcomp,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::fftwlcChanged()
{
    // Update Local contrast GUI according to fftwlc button state
    updateContrastGUI3();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (fftwlc->get_active()) {
                listener->panelChanged(Evlocallabfftwlc,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabfftwlc,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::showmasklcMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabContrast::enalcMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enalcMask->get_active()) {
                listener->panelChanged(EvLocallabEnalcMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnalcMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabContrast::updateContrastGUI1()
{
    const int mode = complexity->get_active_row_number();

    // Update Local contrast GUI according to localcontMethod combobox value
    if (localcontMethod->get_active_row_number() == 0) {
        lcradius->show();
        lcamount->show();
        lcdarkness->show();
        lclightness->show();
        contFrame->hide();
        levelwav->hide();
        expresidpyr->hide();
        clariFrame->hide();
        expcontrastpyr->hide();
        expcontrastpyr2->hide();

        if (mode == Expert) { // Keep widget hidden in Normal and Simple mode
            fftwlc->show();
        }
    } else if (localcontMethod->get_active_row_number() == 1) {
        lcradius->hide();
        lcamount->hide();
        lcdarkness->hide();
        lclightness->hide();
        contFrame->show();
        levelwav->show();
        expresidpyr->show();
        clariFrame->show();

        if (mode == Expert) { // Keep widget hidden in Normal and Simple mode
            expcontrastpyr->show();
            expcontrastpyr2->show();
        }

        fftwlc->hide();
    }
}
void LocallabContrast::updateContrastGUI2()
{
    // Update Local contrast GUI according to waveshow button state
    if (waveshow->get_active()) {
        edgsBoxshow->show();
    } else {
        edgsBoxshow->hide();
    }
}

void LocallabContrast::updateContrastGUI3()
{
    // Update Local contrast GUI according to fftwlc button state
    const double temp = lcradius->getValue();

    if (fftwlc->get_active()) {
        lcradius->setLimits(20, 1000, 1, 80);
    } else {
        lcradius->setLimits(20, 100, 1, 80);
    }

    lcradius->setValue(temp);
}

/* ==== LocallabCBDL ==== */
LocallabCBDL::LocallabCBDL():
    LocallabTool(this, M("TP_LOCALLAB_CBDL_TOOLNAME"), M("TP_LOCALLAB_CBDL"), true),

    // CBDL specific widgets
    levFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LEVFRA")))),
    multiplier([]() -> std::array<Adjuster *, 6>
    {
    std::array<Adjuster*, 6> res = {};

    for (unsigned int i = 0; i < res.size(); ++i) {
        Glib::ustring ss = Glib::ustring::format(i);

        if (i == 0) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMAFINEST"));
        } else if (i == 5) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        }

        res[i] = Gtk::manage(new Adjuster(std::move(ss), 0.0, 4.0, 0.01, 1.0));
    }

    return res;
    }
    ()),
    chromacbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMACBDL"), 0., 1.5, 0.01, 0.))),
    threshold(Gtk::manage(new Adjuster(M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 1., 0.01, 0.2))),
    clarityml(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARITYML"), 0.1, 100., 0.1, 0.1))),
    contresid(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRESID"), -100, 100, 1, 0))),
    softradiuscb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.5, 0.))),
    sensicb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSICB"), 0, 100, 1, 60))),
    expmaskcb(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWCB")))),
    showmaskcbMethod(Gtk::manage(new MyComboBoxText())),
    enacbMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    maskcbCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskcbshape(static_cast<FlatCurveEditor*>(maskcbCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskcbshape(static_cast<FlatCurveEditor*>(maskcbCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskcbshape(static_cast<FlatCurveEditor *>(maskcbCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    lapmaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2cbCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmaskcbshape(static_cast<DiagonalCurveEditor*>(mask2cbCurveEditorG->addCurve(CT_Diagonal, "L(L)"))),

    lumacontrastMinusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")))),
    lumaneutralButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")))),
    lumacontrastPlusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS"))))
{
    const LocallabParams::LocallabSpot defSpot;

    // Parameter CBDL specific widgets
    for (const auto adj : multiplier) {
        adj->setAdjusterListener(this);
    }

    chromacbdl->setAdjusterListener(this);

    threshold->setAdjusterListener(this);

    clarityml->setAdjusterListener(this);

    contresid->setAdjusterListener(this);

    softradiuscb->setLogScale(10, 0);
    softradiuscb->setAdjusterListener(this);

    sensicb->setAdjusterListener(this);

    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskcbMethod->set_active(0);
    showmaskcbMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskcbMethodConn = showmaskcbMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabCBDL::showmaskcbMethodChanged));

    enacbMaskConn = enacbMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabCBDL::enacbMaskChanged));

    maskcbCurveEditorG->setCurveListener(this);

    CCmaskcbshape->setIdentityValue(0.);
    CCmaskcbshape->setResetCurve(FlatCurveType(defSpot.CCmaskcbcurve.at(0)), defSpot.CCmaskcbcurve);
    CCmaskcbshape->setBottomBarColorProvider(this, 1);

    LLmaskcbshape->setIdentityValue(0.);
    LLmaskcbshape->setResetCurve(FlatCurveType(defSpot.LLmaskcbcurve.at(0)), defSpot.LLmaskcbcurve);
    LLmaskcbshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskcbshape->setIdentityValue(0.);
    HHmaskcbshape->setResetCurve(FlatCurveType(defSpot.HHmaskcbcurve.at(0)), defSpot.HHmaskcbcurve);
    HHmaskcbshape->setCurveColorProvider(this, 2);
    HHmaskcbshape->setBottomBarColorProvider(this, 2);

    maskcbCurveEditorG->curveListComplete();

    blendmaskcb->setAdjusterListener(this);

    radmaskcb->setAdjusterListener(this);

    lapmaskcb->setAdjusterListener(this);

    chromaskcb->setAdjusterListener(this);

    gammaskcb->setAdjusterListener(this);

    slomaskcb->setAdjusterListener(this);

    mask2cbCurveEditorG->setCurveListener(this);

    Lmaskcbshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskcbcurve.at(0)), defSpot.Lmaskcbcurve);
    Lmaskcbshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskcbshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2cbCurveEditorG->curveListComplete();

    lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumacontrastMinusPressed));

    lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumaneutralPressed));

    lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumacontrastPlusPressed));

    // Add CBDL specific widgets to GUI
    ToolParamBlock* const levBox = Gtk::manage(new ToolParamBlock());
    Gtk::HBox* buttonBox = Gtk::manage(new Gtk::HBox(true, 10));
    buttonBox->pack_start(*lumacontrastMinusButton);
    buttonBox->pack_start(*lumaneutralButton);
    buttonBox->pack_start(*lumacontrastPlusButton);
    levBox->pack_start(*buttonBox);

    for (const auto adj : multiplier) {
        levBox->pack_start(*adj);
    }

    Gtk::HSeparator* const separator = Gtk::manage(new  Gtk::HSeparator());
    levBox->pack_start(*separator, Gtk::PACK_SHRINK, 2);
    levBox->pack_start(*chromacbdl);
    levBox->pack_start(*threshold);
    levFrame->add(*levBox);
    pack_start(*levFrame);
    Gtk::Frame* const residFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RESID")));
    residFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const residBox = Gtk::manage(new ToolParamBlock());
    residBox->pack_start(*clarityml);
    residBox->pack_start(*contresid);
    residFrame->add(*residBox);
    pack_start(*residFrame);
    pack_start(*softradiuscb);
    pack_start(*sensicb);
    ToolParamBlock* const maskcbBox = Gtk::manage(new ToolParamBlock());
    maskcbBox->pack_start(*showmaskcbMethod, Gtk::PACK_SHRINK, 4);
    maskcbBox->pack_start(*enacbMask, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*maskcbCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskcbBox->pack_start(*blendmaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*radmaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*lapmaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*chromaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*gammaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*slomaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*mask2cbCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    expmaskcb->add(*maskcbBox, false);
    pack_start(*expmaskcb, false, false);
}

LocallabCBDL::~LocallabCBDL()
{
    delete maskcbCurveEditorG;
    delete mask2cbCurveEditorG;
}

bool LocallabCBDL::isMaskViewActive()
{
    return (showmaskcbMethod->get_active_row_number() != 0);
}


void LocallabCBDL::resetMaskView()
{
    showmaskcbMethodConn.block(true);
    showmaskcbMethod->set_active(0);
    showmaskcbMethodConn.block(false);
}

void LocallabCBDL::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &maskMask)
{
    cbMask = showmaskcbMethod->get_active_row_number();
}

void LocallabCBDL::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        levFrame->set_tooltip_text(M("TP_LOCALLAB_EXPCBDL_TOOLTIP"));

        for (const auto adj : multiplier) {
            adj->set_tooltip_text(M("TP_LOCALLAB_CBDL_ADJ_TOOLTIP"));
        }

        chromacbdl->set_tooltip_text(M("TP_LOCALLAB_CHROMACB_TOOLTIP"));
        threshold->set_tooltip_text(M("TP_LOCALLAB_CBDL_THRES_TOOLTIP"));
        clarityml->set_tooltip_text(M("TP_LOCALLAB_CBDLCLARI_TOOLTIP"));
        sensicb->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
        expmaskcb->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        blendmaskcb->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        radmaskcb->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        lapmaskcb->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        mask2cbCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        Lmaskcbshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    } else {
        levFrame->set_tooltip_text("");

        for (const auto adj : multiplier) {
            adj->set_tooltip_text("");
        }

        chromacbdl->set_tooltip_text("");
        threshold->set_tooltip_text("");
        clarityml->set_tooltip_text("");
        sensicb->set_tooltip_text("");
        expmaskcb->set_tooltip_text("");
        CCmaskcbshape->setTooltip("");
        LLmaskcbshape->setTooltip("");
        HHmaskcbshape->setTooltip("");
        blendmaskcb->set_tooltip_text("");
        radmaskcb->set_tooltip_text("");
        lapmaskcb->set_tooltip_text("");
        mask2cbCurveEditorG->set_tooltip_text("");
        Lmaskcbshape->setTooltip("");
    }
}

void LocallabCBDL::setDefaultExpanderVisibility()
{
    expmaskcb->set_expanded(false);
}

void LocallabCBDL::disableListener()
{
    LocallabTool::disableListener();

    showmaskcbMethodConn.block(true);
    enacbMaskConn.block(true);

    lumacontrastMinusPressedConn.block(true);
    lumaneutralPressedConn.block(true);
    lumacontrastPlusPressedConn.block(true);
}

void LocallabCBDL::enableListener()
{
    LocallabTool::enableListener();

    showmaskcbMethodConn.block(false);
    enacbMaskConn.block(false);

    lumacontrastMinusPressedConn.block(false);
    lumaneutralPressedConn.block(false);
    lumacontrastPlusPressedConn.block(false);
}

void LocallabCBDL::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visicbdl);
        exp->setEnabled(spot.expcbdl);
        complexity->set_active(spot.complexcbdl);

        for (int i = 0; i < 6; i++) {
            multiplier[i]->setValue(spot.mult[i]);
        }

        chromacbdl->setValue(spot.chromacbdl);
        threshold->setValue(spot.threshold);
        clarityml->setValue(spot.clarityml);
        contresid->setValue((double)spot.contresid);
        softradiuscb->setValue(spot.softradiuscb);
        sensicb->setValue((double)spot.sensicb);
        enacbMask->set_active(spot.enacbMask);
        CCmaskcbshape->setCurve(spot.CCmaskcbcurve);
        LLmaskcbshape->setCurve(spot.LLmaskcbcurve);
        HHmaskcbshape->setCurve(spot.HHmaskcbcurve);
        blendmaskcb->setValue((double)spot.blendmaskcb);
        radmaskcb->setValue(spot.radmaskcb);
        lapmaskcb->setValue(spot.lapmaskcb);
        chromaskcb->setValue(spot.chromaskcb);
        gammaskcb->setValue(spot.gammaskcb);
        slomaskcb->setValue(spot.slomaskcb);
        Lmaskcbshape->setCurve(spot.Lmaskcbcurve);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabCBDL::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expcbdl = exp->getEnabled();
        spot.visicbdl = exp->get_visible();
        spot.complexcbdl = complexity->get_active_row_number();

        for (int i = 0; i < 6; i++) {
            spot.mult[i] = multiplier[i]->getValue();
        }

        spot.chromacbdl = chromacbdl->getValue();
        spot.threshold = threshold->getValue();
        spot.clarityml = clarityml->getValue();
        spot.contresid = contresid->getIntValue();
        spot.softradiuscb = softradiuscb->getValue();
        spot.sensicb = sensicb->getIntValue();
        spot.enacbMask = enacbMask->get_active();
        spot.LLmaskcbcurve = LLmaskcbshape->getCurve();
        spot.CCmaskcbcurve = CCmaskcbshape->getCurve();
        spot.HHmaskcbcurve = HHmaskcbshape->getCurve();
        spot.blendmaskcb = blendmaskcb->getIntValue();
        spot.radmaskcb = radmaskcb->getValue();
        spot.lapmaskcb = lapmaskcb->getValue();
        spot.chromaskcb = chromaskcb->getValue();
        spot.gammaskcb = gammaskcb->getValue();
        spot.slomaskcb = slomaskcb->getValue();
        spot.Lmaskcbcurve = Lmaskcbshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabCBDL::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        for (int i = 0; i < 6; i++) {
            multiplier[i]->setDefault(defSpot.mult[i]);
        }

        chromacbdl->setDefault(defSpot.chromacbdl);
        threshold->setDefault(defSpot.threshold);
        clarityml->setDefault(defSpot.clarityml);
        contresid->setDefault((double)defSpot.contresid);
        softradiuscb->setDefault(defSpot.softradiuscb);
        sensicb->setDefault((double)defSpot.sensicb);
        blendmaskcb->setDefault((double)defSpot.blendmaskcb);
        radmaskcb->setDefault(defSpot.radmaskcb);
        lapmaskcb->setDefault(defSpot.lapmaskcb);
        chromaskcb->setDefault(defSpot.chromaskcb);
        gammaskcb->setDefault(defSpot.gammaskcb);
        slomaskcb->setDefault(defSpot.slomaskcb);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabCBDL::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == multiplier[0] || a == multiplier[1] || a == multiplier[2] || a == multiplier[3] || a == multiplier[4] || a == multiplier[5]) {
            if (listener) {
                listener->panelChanged(EvlocallabEqualizer,
                                       Glib::ustring::compose("%1, %2, %3, %4, %5, %6",
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[0]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[1]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[2]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[3]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[4]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[5]->getValue()))
                                       + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromacbdl) {
            if (listener) {
                listener->panelChanged(Evlocallabchromacbdl,
                                       chromacbdl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == threshold) {
            if (listener) {
                listener->panelChanged(EvlocallabThresho,
                                       threshold->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == clarityml) {
            if (listener) {
                listener->panelChanged(EvLocallabclarityml,
                                       clarityml->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == contresid) {
            if (listener) {
                listener->panelChanged(EvLocallabcontresid,
                                       contresid->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == softradiuscb) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiuscb,
                                       softradiuscb->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensicb) {
            if (listener) {
                listener->panelChanged(Evlocallabsensicb,
                                       sensicb->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskcb,
                                       blendmaskcb->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskcb,
                                       radmaskcb->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskcb,
                                       lapmaskcb->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskcb,
                                       chromaskcb->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskcb,
                                       gammaskcb->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slomaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskcb,
                                       slomaskcb->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabCBDL::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == CCmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabCBDL::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenacbdl,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenacbdl,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabCBDL::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    lapmaskcb->setValue(defSpot.lapmaskcb);

    // Enable all listeners
    enableListener();
}

void LocallabCBDL::convertParamToSimple()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden specific GUI widgets in Simple mode to default spot values
    softradiuscb->setValue(defSpot.softradiuscb);
    showmaskcbMethod->set_active(0);
    enacbMask->set_active(defSpot.enacbMask);
    CCmaskcbshape->setCurve(defSpot.CCmaskcbcurve);
    LLmaskcbshape->setCurve(defSpot.LLmaskcbcurve);
    HHmaskcbshape->setCurve(defSpot.HHmaskcbcurve);
    blendmaskcb->setValue((double)defSpot.blendmaskcb);
    radmaskcb->setValue(defSpot.radmaskcb);
    chromaskcb->setValue(defSpot.chromaskcb);
    gammaskcb->setValue(defSpot.gammaskcb);
    slomaskcb->setValue(defSpot.slomaskcb);
    Lmaskcbshape->setCurve(defSpot.Lmaskcbcurve);

    // Enable all listers
    enableListener();
}

void LocallabCBDL::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            softradiuscb->hide();
            expmaskcb->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            lapmaskcb->hide();
            // Specific Simple mode widgets are shown in Normal mode
            softradiuscb->show();
            expmaskcb->show();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            softradiuscb->show();
            expmaskcb->show();
            lapmaskcb->show();
    }
}

void LocallabCBDL::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskcbshape->updateLocallabBackground(normChromar);
        LLmaskcbshape->updateLocallabBackground(normLumar);
        HHmaskcbshape->updateLocallabBackground(normHuer);
        Lmaskcbshape->updateLocallabBackground(normLumar);

        return false;
    }
    );
}

void LocallabCBDL::showmaskcbMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabCBDL::enacbMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enacbMask->get_active()) {
                listener->panelChanged(EvLocallabEnacbMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnacbMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabCBDL::lumacontrastMinusPressed()
{
    for (int i = 0; i < 6; i++) {
        float inc = - (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + 0.01f * inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void LocallabCBDL::lumaneutralPressed()
{
    for (int i = 0; i < 6; i++) {
        multiplier[i]->setValue(1.0);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void LocallabCBDL::lumacontrastPlusPressed()
{
    for (int i = 0; i < 6; i++) {
        float inc = (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + 0.01f * inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

/* ==== LocallabLog ==== */
LocallabLog::LocallabLog():
    LocallabTool(this, M("TP_LOCALLAB_LOG_TOOLNAME"), M("TP_LOCALLAB_LOG"), false, false),

    // Log encoding specific widgets
    autocompute(Gtk::manage(new Gtk::ToggleButton(M("TP_LOCALLAB_LOGAUTO")))),
    logPFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOGPFRA")))),
    blackEv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLACK_EV"), -16.0, 0.0, 0.1, -5.0))),
    whiteEv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WHITE_EV"), 0.0, 32.0, 0.1, 10.0))),
    fullimage(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FULLIMAGE")))),
    Autogray(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_AUTOGRAY")))),
    sourceGray(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOURCE_GRAY"), 1.0, 100.0, 0.1, 10.0))),
    targetGray(Gtk::manage(new Adjuster(M("TP_LOCALLAB_TARGET_GRAY"), 5.0, 80.0, 0.1, 18.0))),
    detail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DETAIL"), 0., 1., 0.01, 0.6))),
    baselog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BASELOG"), 1.3, 8., 0.05, 2., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    sensilog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSILOG"), 0, 100, 1, 60))),
    strlog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -2.0, 2.0, 0.05, 0.))),
    anglog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.)))
{
    // Parameter Log encoding specific widgets
    autoconn = autocompute->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::autocomputeToggled));

    blackEv->setLogScale(2, -8);
    blackEv->setAdjusterListener(this);

    whiteEv->setLogScale(16, 0);
    whiteEv->setAdjusterListener(this);

    fullimageConn = fullimage->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::fullimageChanged));

    AutograyConn = Autogray->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::AutograyChanged));

    sourceGray->setAdjusterListener(this);

    targetGray->setAdjusterListener(this);

    detail->setAdjusterListener(this);

    baselog->setAdjusterListener(this);

    sensilog->setAdjusterListener(this);

    strlog->setAdjusterListener(this);

    anglog->setAdjusterListener(this);

    // Add Log encoding specific widgets to GUI
    logPFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const logPBox = Gtk::manage(new ToolParamBlock());
    logPBox->pack_start(*autocompute);
    logPBox->pack_start(*blackEv);
    logPBox->pack_start(*whiteEv);
    logPBox->pack_start(*fullimage);
    logPFrame->add(*logPBox);
    pack_start(*logPFrame);
    Gtk::Frame* const logFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOGFRA")));
    logFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const logFBox = Gtk::manage(new ToolParamBlock());
    logFBox->pack_start(*Autogray);
    logFBox->pack_start(*sourceGray);
    logFrame->add(*logFBox);
    pack_start(*logFrame);
    pack_start(*targetGray);
    pack_start(*detail);
    pack_start(*baselog);
    pack_start(*sensilog);
    Gtk::Frame* const gradlogFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADLOGFRA")));
    gradlogFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const gradlogBox = Gtk::manage(new ToolParamBlock());
    gradlogBox->pack_start(*strlog);
    gradlogBox->pack_start(*anglog);
    gradlogFrame->add(*gradlogBox);
    pack_start(*gradlogFrame);
}

void LocallabLog::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_LOGENCOD_TOOLTIP"));
        logPFrame->set_tooltip_text(M("TP_LOCALLAB_LOGFRAME_TOOLTIP"));
        autocompute->set_tooltip_text(M("TP_LOCALLAB_LOGAUTO_TOOLTIP"));
        blackEv->set_tooltip_text(M("TP_LOCALLAB_LOGBLACKWHEV_TOOLTIP"));
        whiteEv->set_tooltip_text(M("TP_LOCALLAB_LOGBLACKWHEV_TOOLTIP"));
        sourceGray->set_tooltip_text(M("TP_LOCALLAB_LOGSRCGREY_TOOLTIP"));
        targetGray->set_tooltip_text(M("TP_LOCALLAB_LOGTARGGREY_TOOLTIP"));
        baselog->set_tooltip_text(M("TP_LOCALLAB_LOGBASE_TOOLTIP"));
        strlog->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
        anglog->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
    } else {
        exp->set_tooltip_text("");
        logPFrame->set_tooltip_text("");
        autocompute->set_tooltip_text("");
        blackEv->set_tooltip_text("");
        whiteEv->set_tooltip_text("");
        sourceGray->set_tooltip_text("");
        targetGray->set_tooltip_text("");
        baselog->set_tooltip_text("");
        strlog->set_tooltip_text("");
        anglog->set_tooltip_text("");
    }
}

void LocallabLog::disableListener()
{
    LocallabTool::disableListener();

    autoconn.block(true);
    fullimageConn.block(true);
    AutograyConn.block(true);
}

void LocallabLog::enableListener()
{
    LocallabTool::enableListener();

    autoconn.block(false);
    fullimageConn.block(false);
    AutograyConn.block(false);
}

void LocallabLog::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visilog);
        exp->setEnabled(spot.explog);

        autocompute->set_active(spot.autocompute);
        blackEv->setValue(spot.blackEv);
        whiteEv->setValue(spot.whiteEv);
        fullimage->set_active(spot.fullimage);
        Autogray->set_active(spot.Autogray);
        sourceGray->setValue(spot.sourceGray);
        targetGray->setValue(spot.targetGray);
        detail->setValue(spot.detail);
        baselog->setValue(spot.baselog);
        sensilog->setValue((double)spot.sensilog);
        strlog->setValue(spot.strlog);
        anglog->setValue(spot.anglog);
    }

    // Enable all listeners
    enableListener();

    // Update Log Encoding GUI according to autocompute button state
    updateLogGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabLog::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.explog = exp->getEnabled();
        spot.visilog = exp->get_visible();

        spot.autocompute = autocompute->get_active();
        spot.blackEv = blackEv->getValue();
        spot.whiteEv = whiteEv->getValue();
        spot.fullimage = fullimage->get_active();
        spot.Autogray = Autogray->get_active();
        spot.sourceGray = sourceGray->getValue();
        spot.targetGray = targetGray->getValue();
        spot.detail = detail->getValue();
        spot.baselog = baselog->getValue();
        spot.sensilog = sensilog->getIntValue();
        spot.strlog = strlog->getValue();
        spot.anglog = anglog->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabLog::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        // Set default value for adjuster widgets
        blackEv->setDefault(defSpot.blackEv);
        whiteEv->setDefault(defSpot.whiteEv);
        sourceGray->setDefault(defSpot.sourceGray);
        targetGray->setDefault(defSpot.targetGray);
        detail->setDefault(defSpot.detail);
        baselog->setDefault(defSpot.baselog);
        sensilog->setDefault((double)defSpot.sensilog);
        strlog->setDefault(defSpot.strlog);
        anglog->setDefault(defSpot.anglog);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabLog::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == blackEv) {
            if (listener) {
                listener->panelChanged(EvlocallabblackEv,
                                       blackEv->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == whiteEv) {
            if (listener) {
                listener->panelChanged(EvlocallabwhiteEv,
                                       whiteEv->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sourceGray) {
            if (listener) {
                listener->panelChanged(EvlocallabsourceGray,
                                       sourceGray->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == targetGray) {
            if (listener) {
                listener->panelChanged(EvlocallabtargetGray,
                                       targetGray->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == detail) {
            if (listener) {
                listener->panelChanged(Evlocallabdetail,
                                       detail->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == baselog) {
            if (listener) {
                listener->panelChanged(Evlocallabbaselog,
                                       baselog->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensilog) {
            if (listener) {
                listener->panelChanged(Evlocallabsensilog,
                                       sensilog->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strlog) {
            if (listener) {
                listener->panelChanged(Evlocallabstrlog,
                                       strlog->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == anglog) {
            if (listener) {
                listener->panelChanged(Evlocallabanglog,
                                       anglog->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabLog::updateAutocompute(const float blackev, const float whiteev, const float sourceg, const float targetg)
{
    idle_register.add(
    [this, blackev, whiteev, sourceg, targetg]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update adjuster values according to autocomputed ones
        disableListener();

        blackEv->setValue(blackev);
        whiteEv->setValue(whiteev);
        sourceGray->setValue(sourceg);
        targetGray->setValue(targetg);

        enableListener();

        return false;
    }
    );
}

void LocallabLog::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenalog,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenalog,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabLog::autocomputeToggled()
{
    // Update Log Encoding GUI according to autocompute button state
    updateLogGUI();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (autocompute->get_active()) {
                listener->panelChanged(EvLocallabAuto,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabAuto,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabLog::fullimageChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (fullimage->get_active()) {
                listener->panelChanged(Evlocallabfullimage,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabfullimage,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabLog::AutograyChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (Autogray->get_active()) {
                listener->panelChanged(EvlocallabAutogray,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvlocallabAutogray,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabLog::updateLogGUI()
{
    if (autocompute->get_active()) {
        blackEv->set_sensitive(false);
        whiteEv->set_sensitive(false);
        sourceGray->set_sensitive(false);
    } else {
        blackEv->set_sensitive(true);
        whiteEv->set_sensitive(true);
        sourceGray->set_sensitive(true);
    }
}


/* ==== LocallabMask ==== */
LocallabMask::LocallabMask():
    LocallabTool(this, M("TP_LOCALLAB_MASKCOM_TOOLNAME"), M("TP_LOCALLAB_MASKCOM"), false),

    // Common mask specific widgets
    sensimask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
    blendmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKMASK"), -100., 100., 0.1, -10.))),
    blendmaskab(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKMASKAB"), -100., 100., 0.1, -10.))),
    softradiusmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.5, 1.))),
    showmask_Method(Gtk::manage(new MyComboBoxText())),
    enamask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    mask_CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASKCOL"))),
    CCmask_shape(static_cast<FlatCurveEditor*>(mask_CurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmask_shape(static_cast<FlatCurveEditor*>(mask_CurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmask_shape(static_cast<FlatCurveEditor *>(mask_CurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    struFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABSTRUM")))),
    strumaskmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUMASKCOL"), 0., 200., 0.1, 0.))),
    toolmask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TOOLCOL")))),
    blurFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABBLURM")))),
    fftmask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FFTCOL_MASK")))),
    contmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTCOL"), 0., 200., 0.5, 0.))),
    blurmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURCOL"), 0.2, 100., 0.5, 0.2))),
    toolmaskFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_TOOLMASK")))),
    radmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    lapmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slopmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    shadmask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHAMASKCOL"), 0, 100, 1, 0))),
    mask_HCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASKH"))),
    HHhmask_shape(static_cast<FlatCurveEditor *>(mask_HCurveEditorG->addCurve(CT_Flat, "H(H)", nullptr, false, true))),
    mask2CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmask_shape(static_cast<DiagonalCurveEditor*>(mask2CurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    mask2CurveEditorGwav(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVMASK"))),
    LLmask_shapewav(static_cast<FlatCurveEditor*>(mask2CurveEditorGwav->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    csThresholdmask(Gtk::manage(new ThresholdAdjuster(M("TP_LOCALLAB_CSTHRESHOLDBLUR"), 0, 9, 0, 0, 6, 5, 0, false))),
    gradFramemask(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADFRA")))),
    str_mask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -2., 2., 0.05, 0.))),
    ang_mask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180., 180., 0.1, 0.)))
{
    const LocallabParams::LocallabSpot defSpot;

    // Parameter Mask common specific widgets
    sensimask->setAdjusterListener(this);

    blendmask->setLogScale(10, 0);
    blendmask->setAdjusterListener(this);

    blendmaskab->setLogScale(10, 0);
    blendmaskab->setAdjusterListener(this);

    softradiusmask->setAdjusterListener(this);

    showmask_Method->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmask_Method->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmask_Method->append(M("TP_LOCALLAB_SHOWMASK"));
    showmask_Method->append(M("TP_LOCALLAB_SHOWREF"));
    showmask_Method->set_active(0);
    showmask_Method->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmask_MethodConn = showmask_Method->signal_changed().connect(sigc::mem_fun(*this, &LocallabMask::showmask_MethodChanged));

    enamaskConn = enamask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabMask::enamaskChanged));

    mask_CurveEditorG->setCurveListener(this);

    CCmask_shape->setIdentityValue(0.);
    CCmask_shape->setResetCurve(FlatCurveType(defSpot.CCmask_curve.at(0)), defSpot.CCmask_curve);
    CCmask_shape->setBottomBarColorProvider(this, 1);

    LLmask_shape->setIdentityValue(0.);
    LLmask_shape->setResetCurve(FlatCurveType(defSpot.LLmask_curve.at(0)), defSpot.LLmask_curve);
    LLmask_shape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmask_shape->setIdentityValue(0.);
    HHmask_shape->setResetCurve(FlatCurveType(defSpot.HHmask_curve.at(0)), defSpot.HHmask_curve);
    HHmask_shape->setCurveColorProvider(this, 2);
    HHmask_shape->setBottomBarColorProvider(this, 2);

    mask_CurveEditorG->curveListComplete();

    struFrame->set_label_align(0.025, 0.5);

    strumaskmask->setAdjusterListener(this);

    toolmaskConn  = toolmask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabMask::toolmaskChanged));

    blurFrame->set_label_align(0.025, 0.5);

    fftmaskConn = fftmask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabMask::fftmaskChanged));

    contmask->setAdjusterListener(this);

    blurmask->setAdjusterListener(this);

    toolmaskFrame->set_label_align(0.025, 0.5);

    radmask->setAdjusterListener(this);

    lapmask->setAdjusterListener(this);

    chromask->setAdjusterListener(this);

    gammask->setAdjusterListener(this);

    slopmask->setAdjusterListener(this);

    shadmask->setAdjusterListener(this);

    mask_HCurveEditorG->setCurveListener(this);

    HHhmask_shape->setIdentityValue(0.);
    HHhmask_shape->setResetCurve(FlatCurveType(defSpot.HHhmask_curve.at(0)), defSpot.HHhmask_curve);
    HHhmask_shape->setCurveColorProvider(this, 2);
    HHhmask_shape->setBottomBarColorProvider(this, 2);

    mask_HCurveEditorG->curveListComplete();

    mask2CurveEditorG->setCurveListener(this);

    Lmask_shape->setResetCurve(DiagonalCurveType(defSpot.Lmask_curve.at(0)), defSpot.Lmask_curve);
    Lmask_shape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmask_shape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2CurveEditorG->curveListComplete();

    mask2CurveEditorGwav->setCurveListener(this);

    LLmask_shapewav->setIdentityValue(0.);
    LLmask_shapewav->setResetCurve(FlatCurveType(defSpot.LLmask_curvewav.at(0)), defSpot.LLmask_curvewav);
    LLmask_shapewav->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2CurveEditorGwav->curveListComplete();

    csThresholdmask->setAdjusterListener(this);

    gradFramemask->set_label_align(0.025, 0.5);

    str_mask->setAdjusterListener(this);

    ang_mask->setAdjusterListener(this);
    ang_mask->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));

    // Add Common mask specific widgets to GUI
    pack_start(*sensimask);
    pack_start(*blendmask);
    pack_start(*blendmaskab);
    pack_start(*softradiusmask);
    pack_start(*showmask_Method);
    pack_start(*enamask);
    pack_start(*mask_CurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    ToolParamBlock* const strumBox = Gtk::manage(new ToolParamBlock());
    strumBox->pack_start(*strumaskmask);
    strumBox->pack_start(*toolmask);
    struFrame->add(*strumBox);
    pack_start(*struFrame);
    ToolParamBlock* const blurmBox = Gtk::manage(new ToolParamBlock());
    blurmBox->pack_start(*fftmask, Gtk::PACK_SHRINK, 0);
    blurmBox->pack_start(*contmask);
    blurmBox->pack_start(*blurmask);
    blurFrame->add(*blurmBox);
    pack_start(*blurFrame);
    ToolParamBlock* const toolmaskBox = Gtk::manage(new ToolParamBlock());
    toolmaskBox->pack_start(*radmask, Gtk::PACK_SHRINK, 0);
    toolmaskBox->pack_start(*lapmask, Gtk::PACK_SHRINK, 0);
    toolmaskBox->pack_start(*chromask, Gtk::PACK_SHRINK, 0);
    toolmaskBox->pack_start(*gammask, Gtk::PACK_SHRINK, 0);
    toolmaskBox->pack_start(*slopmask, Gtk::PACK_SHRINK, 0);
    toolmaskBox->pack_start(*shadmask, Gtk::PACK_SHRINK, 0);
    toolmaskBox->pack_start(*mask_HCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolmaskBox->pack_start(*mask2CurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolmaskBox->pack_start(*mask2CurveEditorGwav, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolmaskBox->pack_start(*csThresholdmask, Gtk::PACK_SHRINK, 0);
    ToolParamBlock* const gradmaskBox = Gtk::manage(new ToolParamBlock());
    gradmaskBox->pack_start(*str_mask);
    gradmaskBox->pack_start(*ang_mask);
    gradFramemask->add(*gradmaskBox);
    toolmaskBox->pack_start(*gradFramemask, Gtk::PACK_SHRINK, 0);
    toolmaskFrame->add(*toolmaskBox);
    pack_start(*toolmaskFrame);
}

LocallabMask::~LocallabMask()
{
    delete mask_CurveEditorG;
    delete mask_HCurveEditorG;
    delete mask2CurveEditorG;
    delete mask2CurveEditorGwav;
}

bool LocallabMask::isMaskViewActive()
{
    return ((showmask_Method->get_active_row_number() != 0));
}

void LocallabMask::resetMaskView()
{
    showmask_MethodConn.block(true);
    showmask_Method->set_active(0);
    showmask_MethodConn.block(false);
}

void LocallabMask::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &maskMask)
{
    maskMask = showmask_Method->get_active_row_number();
}

void LocallabMask::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_MASKCOM_TOOLTIP"));
        sensimask->set_tooltip_text(M("TP_LOCALLAB_SENSIMASK_TOOLTIP"));
        blendmask->set_tooltip_text(M("TP_LOCALLAB_BLENDMASKMASK_TOOLTIP"));
        blendmaskab->set_tooltip_text(M("TP_LOCALLAB_BLENDMASKMASK_TOOLTIP"));
        CCmask_shape->setTooltip(M("TP_LOCALLAB_CURVEEDITORM_CC_TOOLTIP"));
        LLmask_shape->setTooltip(M("TP_LOCALLAB_CURVEEDITORM_CC_TOOLTIP"));
        HHmask_shape->setTooltip(M("TP_LOCALLAB_CURVEEDITORM_CC_TOOLTIP"));
        struFrame->set_tooltip_text(M("TP_LOCALLAB_STRUMASK_TOOLTIP"));
        radmask->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        lapmask->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        mask_HCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_HHMASK_TOOLTIP"));
        mask2CurveEditorG->set_tooltip_text(M("TP_LOCALLAB_WAVMASK_TOOLTIP"));
        Lmask_shape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        mask2CurveEditorGwav->set_tooltip_text(M("TP_LOCALLAB_WAVMASK_TOOLTIP"));
        LLmask_shapewav->setTooltip(M("TP_LOCALLAB_LMASK_LEVEL_TOOLTIP"));
    } else {
        exp->set_tooltip_text("");
        sensimask->set_tooltip_text("");
        blendmask->set_tooltip_text("");
        blendmaskab->set_tooltip_text("");
        CCmask_shape->setTooltip("");
        LLmask_shape->setTooltip("");
        HHmask_shape->setTooltip("");
        struFrame->set_tooltip_text("");
        radmask->set_tooltip_text("");
        lapmask->set_tooltip_text("");
        mask_HCurveEditorG->set_tooltip_text("");
        mask2CurveEditorG->set_tooltip_text("");
        Lmask_shape->setTooltip("");
        mask2CurveEditorGwav->set_tooltip_text("");
        LLmask_shapewav->setTooltip("");
    }
}

void LocallabMask::disableListener()
{
    LocallabTool::disableListener();

    showmask_MethodConn.block(true);
    enamaskConn.block(true);
    toolmaskConn.block(true);
    fftmaskConn.block(true);
}

void LocallabMask::enableListener()
{
    LocallabTool::enableListener();

    showmask_MethodConn.block(false);
    enamaskConn.block(false);
    toolmaskConn.block(false);
    fftmaskConn.block(false);
}

void LocallabMask::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visimask);
        exp->setEnabled(spot.expmask);
        complexity->set_active(spot.complexmask);

        sensimask->setValue((double)spot.sensimask);
        blendmask->setValue(spot.blendmask);
        blendmaskab->setValue(spot.blendmaskab);
        softradiusmask->setValue(spot.softradiusmask);
        enamask->set_active(spot.enamask);
        CCmask_shape->setCurve(spot.CCmask_curve);
        LLmask_shape->setCurve(spot.LLmask_curve);
        HHmask_shape->setCurve(spot.HHmask_curve);
        strumaskmask->setValue(spot.strumaskmask);
        toolmask->set_active(spot.toolmask);
        fftmask->set_active(spot.fftmask);
        contmask->setValue(spot.contmask);
        // Update Common mask GUI according to fftmask button state
        // Note: Contrary to the others, shall be called before setting blurmask value
        updateMaskGUI();
        blurmask->setValue(spot.blurmask);
        radmask->setValue(spot.radmask);
        lapmask->setValue(spot.lapmask);
        chromask->setValue(spot.chromask);
        gammask->setValue(spot.gammask);
        slopmask->setValue(spot.slopmask);
        shadmask->setValue(spot.shadmask);
        HHhmask_shape->setCurve(spot.HHhmask_curve);
        Lmask_shape->setCurve(spot.Lmask_curve);
        LLmask_shapewav->setCurve(spot.LLmask_curvewav);
        csThresholdmask->setValue<int>(spot.csthresholdmask);
        str_mask->setValue((double)spot.str_mask);
        ang_mask->setValue((double)spot.ang_mask);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabMask::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expmask = exp->getEnabled();
        spot.visimask = exp->get_visible();
        spot.complexmask = complexity->get_active_row_number();

        spot.sensimask = sensimask->getIntValue();
        spot.blendmask = blendmask->getValue();
        spot.blendmaskab = blendmaskab->getValue();
        spot.softradiusmask = softradiusmask->getValue();
        spot.enamask = enamask->get_active();
        spot.CCmask_curve = CCmask_shape->getCurve();
        spot.LLmask_curve = LLmask_shape->getCurve();
        spot.HHmask_curve = HHmask_shape->getCurve();
        spot.strumaskmask = strumaskmask->getValue();
        spot.toolmask = toolmask->get_active();
        spot.fftmask = fftmask->get_active();
        spot.contmask = contmask->getValue();
        spot.blurmask = blurmask->getValue();
        spot.radmask = radmask->getValue();
        spot.lapmask = lapmask->getValue();
        spot.chromask = chromask->getValue();
        spot.gammask = gammask->getValue();
        spot.slopmask = slopmask->getValue();
        spot.shadmask = shadmask->getValue();
        spot.HHhmask_curve = HHhmask_shape->getCurve();
        spot.Lmask_curve = Lmask_shape->getCurve();
        spot.LLmask_curvewav = LLmask_shapewav->getCurve();
        spot.csthresholdmask = csThresholdmask->getValue<int>();
        spot.str_mask = str_mask->getIntValue();
        spot.ang_mask = ang_mask->getIntValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabMask::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        // Set default value for adjuster widgets
        sensimask->setDefault((double)defSpot.sensimask);
        blendmask->setDefault(defSpot.blendmask);
        blendmaskab->setDefault(defSpot.blendmaskab);
        softradiusmask->setDefault(defSpot.softradiusmask);
        strumaskmask->setDefault(defSpot.strumaskmask);
        contmask->setDefault(defSpot.contmask);
        blurmask->setDefault(defSpot.blurmask);
        radmask->setDefault(defSpot.radmask);
        lapmask->setDefault(defSpot.lapmask);
        chromask->setDefault(defSpot.chromask);
        gammask->setDefault(defSpot.lapmask);
        slopmask->setDefault(defSpot.slopmask);
        shadmask->setDefault(defSpot.shadmask);
        csThresholdmask->setDefault<int>(defSpot.csthresholdmask);
        str_mask->setDefault((double)defSpot.str_mask);
        ang_mask->setDefault((double)defSpot.ang_mask);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabMask::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {

        if (a == sensimask) {
            if (listener) {
                listener->panelChanged(Evlocallabsensimask,
                                       sensimask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmask) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmask,
                                       blendmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmaskab) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskab,
                                       blendmaskab->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == softradiusmask) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusmask,
                                       softradiusmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strumaskmask) {
            if (listener) {
                listener->panelChanged(Evlocallabstrumaskmask,
                                       strumaskmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == contmask) {
            if (listener) {
                listener->panelChanged(Evlocallabcontmask,
                                       contmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blurmask) {
            if (listener) {
                listener->panelChanged(Evlocallabblurmask,
                                       blurmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmask) {
            if (listener) {
                listener->panelChanged(Evlocallabradmask,
                                       radmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmask) {
            if (listener) {
                listener->panelChanged(Evlocallablapmask,
                                       lapmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromask) {
            if (listener) {
                listener->panelChanged(Evlocallabchromask,
                                       chromask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammask) {
            if (listener) {
                listener->panelChanged(Evlocallabgammask,
                                       gammask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slopmask) {
            if (listener) {
                listener->panelChanged(Evlocallabslopmask,
                                       slopmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == shadmask) {
            if (listener) {
                listener->panelChanged(Evlocallabshadmask,
                                       shadmask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == str_mask) {
            if (listener) {
                listener->panelChanged(Evlocallabstr_mask,
                                       str_mask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == ang_mask) {
            if (listener) {
                listener->panelChanged(Evlocallabang_mask,
                                       ang_mask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

    }
}

void LocallabMask::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == csThresholdmask) {
            if (listener) {
                listener->panelChanged(EvlocallabcsThresholdmask,
                                       csThresholdmask->getHistoryString() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabMask::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == CCmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHhmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHhmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmask_shapewav) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmask_shapewav,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

    }
}

void LocallabMask::complexityModeChanged()
{
    if (complexity->get_active_row_number() == Simple) { // New selected mode is Simple one
        // Convert tool widget parameters
        convertParamToNormal(); // From Expert mode to Normal mode
        convertParamToSimple(); // From Normal mode to Simple mode
        // Update GUI based on new mode
        updateGUIToMode(Simple);

        if (listener && isLocActivated) {
            listener->panelChanged(EvlocallabcomplexityWithRefresh,
                                   M("TP_LOCALLAB_MODE_SIMPLE") + " (" + escapeHtmlChars(spotName) + ")");
        }
    } else if (complexity->get_active_row_number() == Normal) { // New selected mode is Normal one
        // Convert tool widget parameters
        convertParamToNormal();
        // Update GUI based on new mode
        updateGUIToMode(Normal);

        if (listener && isLocActivated) {
            listener->panelChanged(EvlocallabcomplexityWithRefresh,
                                   M("TP_LOCALLAB_MODE_NORMAL") + " (" + escapeHtmlChars(spotName) + ")");
        }
    } else if (complexity->get_active_row_number() == Expert) { // New selected mode is Expert one
        // Update GUI based on new mode
        updateGUIToMode(Expert);

        if (listener && isLocActivated) {
            listener->panelChanged(EvlocallabcomplexityWithRefresh,
                                   M("TP_LOCALLAB_MODE_EXPERT") + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabMask::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocena_mask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocena_mask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabMask::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    softradiusmask->setValue(defSpot.softradiusmask);
    strumaskmask->setValue(defSpot.strumaskmask);
    toolmask->set_active(defSpot.toolmask);
    fftmask->set_active(defSpot.fftmask);
    contmask->setValue(defSpot.contmask);
    blurmask->setValue(defSpot.blurmask);
    lapmask->setValue(defSpot.lapmask);
    gammask->setValue(defSpot.gammask);
    slopmask->setValue(defSpot.slopmask);
    shadmask->setValue(defSpot.shadmask);
    HHhmask_shape->setCurve(defSpot.HHhmask_curve);
    LLmask_shapewav->setCurve(defSpot.LLmask_curvewav);
    csThresholdmask->setValue<int>(defSpot.csthresholdmask);
    str_mask->setValue((double)defSpot.str_mask);
    ang_mask->setValue((double)defSpot.ang_mask);

    // Enable all listeners
    enableListener();

    // Update GUI based on converted widget parameters:
    // - Update Common mask GUI according to fftmask button state
    updateMaskGUI();
}

void LocallabMask::convertParamToSimple()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden specific GUI widgets in Simple mode to default spot values
    radmask->setValue(defSpot.radmask);
    chromask->setValue(defSpot.chromask);
    Lmask_shape->setCurve(defSpot.Lmask_curve);

    // Enable all listeners
    enableListener();
}

void LocallabMask::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            softradiusmask->hide();
            struFrame->hide();
            blurFrame->hide();
            toolmaskFrame->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            softradiusmask->hide();
            struFrame->hide();
            blurFrame->hide();
            lapmask->hide();
            gammask->hide();
            slopmask->hide();
            shadmask->hide();
            mask_HCurveEditorG->hide();
            mask2CurveEditorGwav->hide();
            csThresholdmask->hide();
            gradFramemask->hide();
            // Specific Simple mode widgets are shown in Normal mode
            toolmaskFrame->show();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            softradiusmask->hide();
            struFrame->hide();
            blurFrame->hide();
            toolmaskFrame->hide();
            lapmask->hide();
            gammask->hide();
            slopmask->hide();
            shadmask->hide();
            mask_HCurveEditorG->hide();
            mask2CurveEditorGwav->hide();
            csThresholdmask->hide();
            gradFramemask->hide();
    }
}

void LocallabMask::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmask_shape->updateLocallabBackground(normChromar);
        LLmask_shape->updateLocallabBackground(normLumar);
        HHmask_shape->updateLocallabBackground(normHuer);
        HHhmask_shape->updateLocallabBackground(normHuer);
        Lmask_shape->updateLocallabBackground(normLumar);

        return false;
    }
    );
}

void LocallabMask::showmask_MethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabMask::enamaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enamask->get_active()) {
                listener->panelChanged(EvLocallabEnaMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabMask::toolmaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (toolmask->get_active()) {
                listener->panelChanged(EvLocallabtoolmask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabtoolmask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabMask::fftmaskChanged()
{
    // Update Common mask GUI according to fftmask button state
    updateMaskGUI();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (fftmask->get_active()) {
                listener->panelChanged(EvLocallabfftmask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabfftmask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabMask::updateMaskGUI()
{
    const double temp = blurmask->getValue();

    if (fftmask->get_active()) {
        blurmask->setLimits(0.2, 1000., 0.5, 0.2);
    } else {
        blurmask->setLimits(0.2, 100., 0.5, 0.2);
    }

    blurmask->setValue(temp);
}
