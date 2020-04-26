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
 *  2019 Pierre Cabrera <pierre.cab@gmail.com>
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
    softradiustm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), -10.0, 1000.0, 0.1, 0.))),
    sensitm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
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
    radmasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
    chromasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.05, 5.0, 0.01, 1.))),
    slomasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2tmCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmasktmshape(static_cast<DiagonalCurveEditor*>(mask2tmCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
{
    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    const LocallabParams::LocallabSpot defSpot;

    // Parameter Tone Mapping specific widgets
    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_TONEMAP_TOOLTIP"));
    }

    amount->setAdjusterListener(this);

    stren->setAdjusterListener(this);

    equiltmConn = equiltm->signal_toggled().connect(sigc::mem_fun(*this, &LocallabTone::equiltmChanged));

    gamma->setAdjusterListener(this);

    satur->setAdjusterListener(this);

    estop->setAdjusterListener(this);

    if (showtooltip) {
        estop->set_tooltip_text(M("TP_LOCALLAB_TONEMAPESTOP_TOOLTIP"));
    }

    scaltm->setAdjusterListener(this);

    rewei->setAdjusterListener(this);

    if (showtooltip) {
        rewei->set_tooltip_text(M("TP_LOCALLAB_TONEMAPESTOP_TOOLTIP"));
    }

    softradiustm->setLogScale(10, -10);
    softradiustm->setAdjusterListener(this);

    sensitm->setAdjusterListener(this);

    if (showtooltip) {
        sensitm->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    }

    setExpandAlignProperties(expmasktm, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        expmasktm->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    showmasktmMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasktmMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmasktmMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmasktmMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmasktmMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmasktmMethod->set_active(0);

    if (showtooltip) {
        showmasktmMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmasktmMethodConn = showmasktmMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabTone::showmasktmMethodChanged));

    enatmMaskConn = enatmMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabTone::enatmMaskChanged));

    enatmMaskaftConn = enatmMaskaft->signal_toggled().connect(sigc::mem_fun(*this, &LocallabTone::enatmMaskaftChanged));

    masktmCurveEditorG->setCurveListener(this);

    CCmasktmshape->setIdentityValue(0.);
    CCmasktmshape->setResetCurve(FlatCurveType(defSpot.CCmasktmcurve.at(0)), defSpot.CCmasktmcurve);

    if (showtooltip) {
        CCmasktmshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmasktmshape->setBottomBarColorProvider(this, 1);

    LLmasktmshape->setIdentityValue(0.);
    LLmasktmshape->setResetCurve(FlatCurveType(defSpot.LLmasktmcurve.at(0)), defSpot.LLmasktmcurve);

    if (showtooltip) {
        LLmasktmshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmasktmshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmasktmshape->setIdentityValue(0.);
    HHmasktmshape->setResetCurve(FlatCurveType(defSpot.HHmasktmcurve.at(0)), defSpot.HHmasktmcurve);

    if (showtooltip) {
        HHmasktmshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmasktmshape->setCurveColorProvider(this, 2);
    HHmasktmshape->setBottomBarColorProvider(this, 2);

    masktmCurveEditorG->curveListComplete();

    blendmasktm->setAdjusterListener(this);

    lapmasktm->setAdjusterListener(this);

    if (showtooltip) {
        lapmasktm->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    radmasktm->setLogScale(10, -10);
    radmasktm->setAdjusterListener(this);

    if (showtooltip) {
        radmasktm->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    chromasktm->setAdjusterListener(this);

    gammasktm->setAdjusterListener(this);

    slomasktm->setAdjusterListener(this);

    mask2tmCurveEditorG->setCurveListener(this);
    Lmasktmshape->setResetCurve(DiagonalCurveType(defSpot.Lmasktmcurve.at(0)), defSpot.Lmasktmcurve);

    if (showtooltip) {
        Lmasktmshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    Lmasktmshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmasktmshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2tmCurveEditorG->curveListComplete();

    // Add Tone Mapping specific widgets to GUI
    // pack_start(*amount); // To use if we change transit_shapedetect parameters
    pack_start(*stren);
    pack_start(*equiltm);

    if (complexsoft < 2) {
        pack_start(*gamma);
        pack_start(*satur);
    }

    pack_start(*estop);
    pack_start(*scaltm);

    if (complexsoft < 2) {
        pack_start(*rewei);
    }

    // pack_start(*softradiustm); // Always bad with TM ??
    pack_start(*sensitm);
    ToolParamBlock* const masktmBox = Gtk::manage(new ToolParamBlock());
    masktmBox->pack_start(*showmasktmMethod, Gtk::PACK_SHRINK, 4);
    masktmBox->pack_start(*enatmMask, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*enatmMaskaft, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*masktmCurveEditorG, Gtk::PACK_SHRINK, 4);
    masktmBox->pack_start(*blendmasktm, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 1) {
        masktmBox->pack_start(*lapmasktm, Gtk::PACK_SHRINK, 0);
    }

    masktmBox->pack_start(*radmasktm, Gtk::PACK_SHRINK, 0);
    masktmBox->pack_start(*chromasktm, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        masktmBox->pack_start(*gammasktm, Gtk::PACK_SHRINK, 0);
        masktmBox->pack_start(*slomasktm, Gtk::PACK_SHRINK, 0);
    }

    masktmBox->pack_start(*mask2tmCurveEditorG, Gtk::PACK_SHRINK, 4);
    expmasktm->add(*masktmBox, false);
    pack_start(*expmasktm, false, false);
}

LocallabTone::~LocallabTone()
{
    delete masktmCurveEditorG;
    delete mask2tmCurveEditorG;
}

void LocallabTone::resetMaskView()
{
    showmasktmMethodConn.block(true);
    showmasktmMethod->set_active(0);
    showmasktmMethodConn.block(false);
}

void LocallabTone::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    tmMask = showmasktmMethod->get_active_row_number();
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
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visitonemap);
        exp->setEnabled(pp->locallab.spots.at(index).exptonemap);

        amount->setValue(pp->locallab.spots.at(index).amount);
        stren->setValue(pp->locallab.spots.at(index).stren);
        equiltm->set_active(pp->locallab.spots.at(index).equiltm);

        if (complexsoft < 2) {
            gamma->setValue(pp->locallab.spots.at(index).gamma);
            satur->setValue(pp->locallab.spots.at(index).satur);
        } else {
            gamma->setValue(1.);
            satur->setValue(0.);
        }

        estop->setValue(pp->locallab.spots.at(index).estop);
        scaltm->setValue(pp->locallab.spots.at(index).scaltm);

        if (complexsoft < 2) {
            rewei->setValue((double)pp->locallab.spots.at(index).rewei);
        } else {
            rewei->setValue(0.);
        }

        softradiustm->setValue(pp->locallab.spots.at(index).softradiustm);
        sensitm->setValue((double)pp->locallab.spots.at(index).sensitm);
        enatmMask->set_active(pp->locallab.spots.at(index).enatmMask);
        enatmMaskaft->set_active(pp->locallab.spots.at(index).enatmMaskaft);
        CCmasktmshape->setCurve(pp->locallab.spots.at(index).CCmasktmcurve);
        LLmasktmshape->setCurve(pp->locallab.spots.at(index).LLmasktmcurve);
        HHmasktmshape->setCurve(pp->locallab.spots.at(index).HHmasktmcurve);
        blendmasktm->setValue((double)pp->locallab.spots.at(index).blendmasktm);

        if (complexsoft == 0) {
            lapmasktm->setValue(pp->locallab.spots.at(index).lapmasktm);
        } else {
            lapmasktm->setValue(0.);
        }

        radmasktm->setValue(pp->locallab.spots.at(index).radmasktm);
        chromasktm->setValue(pp->locallab.spots.at(index).chromasktm);

        if (complexsoft < 2) {
            gammasktm->setValue(pp->locallab.spots.at(index).gammasktm);
            slomasktm->setValue(pp->locallab.spots.at(index).slomasktm);
        } else {
            gammasktm->setValue(1.);
            slomasktm->setValue(0.);
        }

        Lmasktmshape->setCurve(pp->locallab.spots.at(index).Lmasktmcurve);
    }

    // Enable all listeners
    enableListener();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabTone::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).exptonemap = exp->getEnabled();
        pp->locallab.spots.at(index).visitonemap = exp->get_visible();

        pp->locallab.spots.at(index).amount = amount->getValue();
        pp->locallab.spots.at(index).stren = stren->getValue();
        pp->locallab.spots.at(index).equiltm = equiltm->get_active();
        pp->locallab.spots.at(index).gamma = gamma->getValue();
        pp->locallab.spots.at(index).satur = satur->getValue();
        pp->locallab.spots.at(index).estop = estop->getValue();
        pp->locallab.spots.at(index).scaltm = scaltm->getValue();
        pp->locallab.spots.at(index).rewei = rewei->getIntValue();
        pp->locallab.spots.at(index).softradiustm = softradiustm->getValue();
        pp->locallab.spots.at(index).sensitm = sensitm->getIntValue();
        pp->locallab.spots.at(index).enatmMask = enatmMask->get_active();
        pp->locallab.spots.at(index).enatmMaskaft = enatmMaskaft->get_active();
        pp->locallab.spots.at(index).LLmasktmcurve = LLmasktmshape->getCurve();
        pp->locallab.spots.at(index).CCmasktmcurve = CCmasktmshape->getCurve();
        pp->locallab.spots.at(index).HHmasktmcurve = HHmasktmshape->getCurve();
        pp->locallab.spots.at(index).blendmasktm = blendmasktm->getIntValue();
        pp->locallab.spots.at(index).lapmasktm = lapmasktm->getValue();
        pp->locallab.spots.at(index).radmasktm = radmasktm->getValue();
        pp->locallab.spots.at(index).chromasktm = chromasktm->getValue();
        pp->locallab.spots.at(index).gammasktm = gammasktm->getValue();
        pp->locallab.spots.at(index).slomasktm = slomasktm->getValue();
        pp->locallab.spots.at(index).Lmasktmcurve = Lmasktmshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabTone::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

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
    if (isLocActivated && exp->getEnabled()) {
        if (a == amount) {
            if (listener) {
                listener->panelChanged(Evlocallabamount,
                                       amount->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == stren) {
            if (listener) {
                listener->panelChanged(Evlocallabstren,
                                       stren->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamma) {
            if (listener) {
                listener->panelChanged(Evlocallabgamma,
                                       gamma->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == satur) {
            if (listener) {
                listener->panelChanged(Evlocallabsatur,
                                       satur->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == estop) {
            if (listener) {
                listener->panelChanged(Evlocallabestop,
                                       estop->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == scaltm) {
            if (listener) {
                listener->panelChanged(Evlocallabscaltm,
                                       scaltm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == rewei) {
            if (listener) {
                listener->panelChanged(Evlocallabrewei,
                                       rewei->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == softradiustm) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiustm,
                                       softradiustm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensitm) {
            if (listener) {
                listener->panelChanged(Evlocallabsensitm,
                                       sensitm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmasktm) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmasktm,
                                       blendmasktm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmasktm) {
            if (listener) {
                listener->panelChanged(Evlocallablapmasktm,
                                       lapmasktm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmasktm) {
            if (listener) {
                listener->panelChanged(Evlocallabradmasktm,
                                       radmasktm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromasktm) {
            if (listener) {
                listener->panelChanged(Evlocallabchromasktm,
                                       chromasktm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammasktm) {
            if (listener) {
                listener->panelChanged(Evlocallabgammasktm,
                                       gammasktm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slomasktm) {
            if (listener) {
                listener->panelChanged(Evlocallabslomasktm,
                                       slomasktm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabTone::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == CCmasktmshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmasktmshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmasktmshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmasktmshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmasktmshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmasktmshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmasktmshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmasktmshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabTone::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenatonemap,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenatonemap,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
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
    dehaFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_DEHAFRA")))),
    dehaz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DEHAZ"), -100, 100, 1, 0))),
    depth(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DEPTH"), 0, 100, 1, 25))),
    lumonly(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_LUMONLY")))),
    retiFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RETIFRA")))),
    str(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STR"), 0., 100., 0.1, 0.2))),
    loglin(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_LOGLIN")))),
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIH"), 0, 100, 1, 60))),
    retitoolFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RETITOOLFRA")))),
    retiBox(Gtk::manage(new ToolParamBlock())),
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
    softradiusret(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRETI"), -10.0, 1000.0, 0.5, 40.))),
    LocalcurveEditortransT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONMAP"))),
    cTtransshape(static_cast<FlatCurveEditor*>(LocalcurveEditortransT->addCurve(CT_Flat, "", nullptr, false, false))),
    mMLabels(Gtk::manage(new Gtk::Label("---"))),
    transLabels(Gtk::manage(new Gtk::Label("---"))),
    transLabels2(Gtk::manage(new Gtk::Label("---"))),
    LocalcurveEditorgainT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAIN"))),
    cTgainshape(static_cast<FlatCurveEditor*>(LocalcurveEditorgainT->addCurve(CT_Flat, "", nullptr, false, false))),
    expmaskreti(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWR")))),
    maskretiBox(Gtk::manage(new ToolParamBlock())),
    showmaskretiMethod(Gtk::manage(new MyComboBoxText())),
    enaretiMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    enaretiMasktmap(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TM_MASK")))),
    maskretiCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskretishape(static_cast<FlatCurveEditor*>(maskretiCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskretishape(static_cast<FlatCurveEditor*>(maskretiCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskretishape(static_cast<FlatCurveEditor *>(maskretiCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 10.))),
    lapmaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.05, 5.0, 0.01, 1.))),
    slomaskreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2retiCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmaskretishape(static_cast<DiagonalCurveEditor*>(mask2retiCurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    inversret(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS"))))
{
    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    const LocallabParams::LocallabSpot defSpot;

    // Parameter Retinex specific widgets
    dehaFrame->set_label_align(0.025, 0.5);

    if (showtooltip) {
        dehaz->set_tooltip_text(M("TP_LOCALLAB_DEHAZ_TOOLTIP"));
    }

    dehaz->setAdjusterListener(this);

    depth->setAdjusterListener(this);

    lumonlyConn = lumonly->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::lumonlyChanged));

    retiFrame->set_label_align(0.025, 0.5);

    str->setAdjusterListener(this);

    loglinConn = loglin->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::loglinChanged));

    if (showtooltip) {
        sensih->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    }

    sensih->setAdjusterListener(this);

    retitoolFrame->set_label_align(0.025, 0.5);

    retinexMethod->append(M("TP_RETINEX_LOW"));
    retinexMethod->append(M("TP_RETINEX_UNIFORM"));
    retinexMethod->append(M("TP_RETINEX_HIGH"));
    retinexMethod->set_active(0);

    if (showtooltip) {
        retinexMethod->set_tooltip_markup(M("TP_LOCRETI_METHOD_TOOLTIP"));
    }

    retinexMethodConn = retinexMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabRetinex::retinexMethodChanged));

    if (showtooltip) {
        fftwreti->set_tooltip_text(M("TP_LOCALLAB_RETI_FFTW_TOOLTIP"));
    }

    fftwretiConn = fftwreti->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::fftwretiChanged));

    equilretConn = equilret->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::equilretChanged));

    if (showtooltip) {
        neigh->set_tooltip_text(M("TP_LOCALLAB_RETI_NEIGH_VART_TOOLTIP"));
    }

    neigh->setAdjusterListener(this);

    if (showtooltip) {
        vart->set_tooltip_text(M("TP_LOCALLAB_RETI_NEIGH_VART_TOOLTIP"));
    }

    vart->setAdjusterListener(this);

    scalereti->setAdjusterListener(this);

    limd->setAdjusterListener(this);

    offs->setAdjusterListener(this);

    setExpandAlignProperties(expretitools, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    chrrt->setAdjusterListener(this);

    darkness->setAdjusterListener(this);

    lightnessreti->setAdjusterListener(this);

    cliptm->setAdjusterListener(this);

    if (showtooltip) {
        softradiusret->set_tooltip_text(M("TP_LOCALLAB_GUIDFILTER_TOOLTIP"));
    }

    softradiusret->setLogScale(10, -10);
    softradiusret->setAdjusterListener(this);

    LocalcurveEditortransT->setCurveListener(this);

    cTtransshape->setIdentityValue(0.);
    cTtransshape->setResetCurve(FlatCurveType(defSpot.localTtranscurve.at(0)), defSpot.localTtranscurve);

    if (showtooltip) {
        cTtransshape->setTooltip(M("TP_LOCALLAB_TRANSMISSION_TOOLTIP"));
    }

    LocalcurveEditortransT->curveListComplete();

    if (showtooltip) {
        mMLabels->set_tooltip_markup(M("TP_LOCALLAB_MLABEL_TOOLTIP"));
    }

    setExpandAlignProperties(mMLabels, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);

    if (showtooltip) {
        transLabels->set_tooltip_markup(M("TP_LOCALLAB_TLABEL_TOOLTIP"));
    }

    setExpandAlignProperties(transLabels, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);

    setExpandAlignProperties(transLabels2, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);

    LocalcurveEditorgainT->setCurveListener(this);

    cTgainshape->setIdentityValue(0.);
    cTgainshape->setResetCurve(FlatCurveType(defSpot.localTgaincurve.at(0)), defSpot.localTgaincurve);

    if (showtooltip) {
        cTgainshape->setTooltip(M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));
    }

    LocalcurveEditorgainT->curveListComplete();

    if (showtooltip) {
        expmaskreti->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    setExpandAlignProperties(expmaskreti, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskretiMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskretiMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmaskretiMethod->set_active(0);

    if (showtooltip) {
        showmaskretiMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskretiMethodConn = showmaskretiMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabRetinex::showmaskretiMethodChanged));

    enaretiMaskConn = enaretiMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::enaretiMaskChanged));

    if (showtooltip) {
        enaretiMasktmap->set_tooltip_markup(M("TP_LOCALLAB_ENARETIMASKTMAP_TOOLTIP"));
    }

    enaretiMasktmapConn = enaretiMasktmap->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::enaretiMasktmapChanged));

    maskretiCurveEditorG->setCurveListener(this);

    CCmaskretishape->setIdentityValue(0.);
    CCmaskretishape->setResetCurve(FlatCurveType(defSpot.CCmaskreticurve.at(0)), defSpot.CCmaskreticurve);

    if (showtooltip) {
        CCmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmaskretishape->setBottomBarColorProvider(this, 1);

    LLmaskretishape->setIdentityValue(0.);
    LLmaskretishape->setResetCurve(FlatCurveType(defSpot.LLmaskreticurve.at(0)), defSpot.LLmaskreticurve);

    if (showtooltip) {
        LLmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmaskretishape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskretishape->setIdentityValue(0.);
    HHmaskretishape->setResetCurve(FlatCurveType(defSpot.HHmaskreticurve.at(0)), defSpot.HHmaskreticurve);

    if (showtooltip) {
        HHmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmaskretishape->setCurveColorProvider(this, 2);
    HHmaskretishape->setBottomBarColorProvider(this, 2);

    maskretiCurveEditorG->curveListComplete();

    blendmaskreti->setAdjusterListener(this);

    if (showtooltip) {
        radmaskreti->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    radmaskreti->setAdjusterListener(this);

    if (showtooltip) {
        lapmaskreti->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    lapmaskreti->setAdjusterListener(this);

    chromaskreti->setAdjusterListener(this);

    gammaskreti->setAdjusterListener(this);

    slomaskreti->setAdjusterListener(this);

    mask2retiCurveEditorG->setCurveListener(this);

    Lmaskretishape->setResetCurve(DiagonalCurveType(defSpot.Lmaskreticurve.at(0)), defSpot.Lmaskreticurve);

    if (showtooltip) {
        Lmaskretishape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    Lmaskretishape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskretishape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2retiCurveEditorG->curveListComplete();

    inversretConn = inversret->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::inversretChanged));

    // Add Retinex specific widgets to GUI
    ToolParamBlock* const auxBox = Gtk::manage(new ToolParamBlock());
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

    if (complexsoft < 1) {
        auxBox->add(*retiFrame);
    }

    ToolParamBlock* const scopeBox = Gtk::manage(new ToolParamBlock());
    scopeBox->pack_start(*sensih);
    auxBox->add(*scopeBox);
    pack_start(*auxBox);
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

void LocallabRetinex::resetMaskView()
{
    showmaskretiMethodConn.block(true);
    showmaskretiMethod->set_active(0);
    showmaskretiMethodConn.block(false);
}

void LocallabRetinex::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    retiMask = showmaskretiMethod->get_active_row_number();
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
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visireti);
        exp->setEnabled(pp->locallab.spots.at(index).expreti);

        dehaz->setValue((double)pp->locallab.spots.at(index).dehaz);
        depth->setValue((double)pp->locallab.spots.at(index).depth);
        lumonly->set_active(pp->locallab.spots.at(index).lumonly);

        if (complexsoft < 2) {
            str->setValue(pp->locallab.spots.at(index).str);
        } else {
            str->setValue(0.);
        }

        loglin->set_active(pp->locallab.spots.at(index).loglin);
        sensih->setValue((double)pp->locallab.spots.at(index).sensih);

        if (pp->locallab.spots.at(index).retinexMethod == "low") {
            retinexMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).retinexMethod == "uni") {
            retinexMethod->set_active(1);
        } else {
            retinexMethod->set_active(2);
        }

        fftwreti->set_active(pp->locallab.spots.at(index).fftwreti);
        equilret->set_active(pp->locallab.spots.at(index).equilret);
        neigh->setValue(pp->locallab.spots.at(index).neigh);
        vart->setValue(pp->locallab.spots.at(index).vart);
        scalereti->setValue(pp->locallab.spots.at(index).scalereti);
        limd->setValue(pp->locallab.spots.at(index).limd);
        offs->setValue(pp->locallab.spots.at(index).offs);
        chrrt->setValue(pp->locallab.spots.at(index).chrrt);
        darkness->setValue(pp->locallab.spots.at(index).darkness);
        lightnessreti->setValue(pp->locallab.spots.at(index).lightnessreti);
        cliptm->setValue(pp->locallab.spots.at(index).cliptm);
        softradiusret->setValue(pp->locallab.spots.at(index).softradiusret);
        cTtransshape->setCurve(pp->locallab.spots.at(index).localTtranscurve);
        cTgainshape->setCurve(pp->locallab.spots.at(index).localTgaincurve);
        enaretiMask->set_active(pp->locallab.spots.at(index).enaretiMask);
        enaretiMasktmap->set_active(pp->locallab.spots.at(index).enaretiMasktmap);
        CCmaskretishape->setCurve(pp->locallab.spots.at(index).CCmaskreticurve);
        LLmaskretishape->setCurve(pp->locallab.spots.at(index).LLmaskreticurve);
        HHmaskretishape->setCurve(pp->locallab.spots.at(index).HHmaskreticurve);
        blendmaskreti->setValue((double)pp->locallab.spots.at(index).blendmaskreti);
        radmaskreti->setValue(pp->locallab.spots.at(index).radmaskreti);
        lapmaskreti->setValue(pp->locallab.spots.at(index).lapmaskreti);
        chromaskreti->setValue(pp->locallab.spots.at(index).chromaskreti);
        gammaskreti->setValue(pp->locallab.spots.at(index).gammaskreti);
        slomaskreti->setValue(pp->locallab.spots.at(index).slomaskreti);
        Lmaskretishape->setCurve(pp->locallab.spots.at(index).Lmaskreticurve);
        inversret->set_active(pp->locallab.spots.at(index).inversret);
    }

    // Enable all listeners
    enableListener();

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
        pp->locallab.spots.at(index).expreti = exp->getEnabled();
        pp->locallab.spots.at(index).visireti = exp->get_visible();

        pp->locallab.spots.at(index).dehaz = dehaz->getIntValue();
        pp->locallab.spots.at(index).depth = depth->getIntValue();
        pp->locallab.spots.at(index).lumonly = lumonly->get_active();
        pp->locallab.spots.at(index).str = str->getValue();
        pp->locallab.spots.at(index).loglin = loglin->get_active();
        pp->locallab.spots.at(index).sensih = sensih->getIntValue();

        if (retinexMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).retinexMethod = "low";
        } else if (retinexMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).retinexMethod = "uni";
        } else if (retinexMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).retinexMethod = "high";
        }

        pp->locallab.spots.at(index).fftwreti = fftwreti->get_active();
        pp->locallab.spots.at(index).equilret = equilret->get_active();
        pp->locallab.spots.at(index).neigh = neigh->getValue();
        pp->locallab.spots.at(index).vart = vart->getValue();
        pp->locallab.spots.at(index).scalereti = scalereti->getValue();
        pp->locallab.spots.at(index).limd = limd->getValue();
        pp->locallab.spots.at(index).offs = offs->getValue();
        pp->locallab.spots.at(index).chrrt = chrrt->getValue();
        pp->locallab.spots.at(index).darkness = darkness->getValue();
        pp->locallab.spots.at(index).lightnessreti = lightnessreti->getValue();
        pp->locallab.spots.at(index).cliptm = cliptm->getValue();
        pp->locallab.spots.at(index).softradiusret = softradiusret->getValue();
        pp->locallab.spots.at(index).localTtranscurve = cTtransshape->getCurve();
        pp->locallab.spots.at(index).localTgaincurve = cTgainshape->getCurve();
        pp->locallab.spots.at(index).enaretiMask = enaretiMask->get_active();
        pp->locallab.spots.at(index).enaretiMasktmap = enaretiMasktmap->get_active();
        pp->locallab.spots.at(index).CCmaskreticurve = CCmaskretishape->getCurve();
        pp->locallab.spots.at(index).LLmaskreticurve = LLmaskretishape->getCurve();
        pp->locallab.spots.at(index).HHmaskreticurve = HHmaskretishape->getCurve();
        pp->locallab.spots.at(index).blendmaskreti = blendmaskreti->getIntValue();
        pp->locallab.spots.at(index).radmaskreti = radmaskreti->getValue();
        pp->locallab.spots.at(index).lapmaskreti = lapmaskreti->getValue();
        pp->locallab.spots.at(index).chromaskreti = chromaskreti->getValue();
        pp->locallab.spots.at(index).gammaskreti = gammaskreti->getValue();
        pp->locallab.spots.at(index).slomaskreti = slomaskreti->getValue();
        pp->locallab.spots.at(index).Lmaskreticurve = Lmaskretishape->getCurve();
        pp->locallab.spots.at(index).inversret = inversret->get_active();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabRetinex::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

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

void LocallabRetinex::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskretishape->updateLocallabBackground(normChromar);
        LLmaskretishape->updateLocallabBackground(normLumar);
        HHmaskretishape->updateLocallabBackground(normHuer);

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
    // Update Retinex GUI according to inversret button state
    updateRetinexGUI2();

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
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 0.4, 2.5, 0.01, 0.75))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 100))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 0))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sharblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARBLUR"), 0.2, 2.0, 0.05, 0.2))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    inverssha(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    sharFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHARFRAME")))),
    showmasksharMethod(Gtk::manage(new MyComboBoxText()))
{
    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    // Parameter Sharpening specific widgets
    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPSHARP_TOOLTIP"));
    }

    sharcontrast->setAdjusterListener(this);

    sharradius->setAdjusterListener(this);

    sharamount->setAdjusterListener(this);

    shardamping->setAdjusterListener(this);

    shariter->setAdjusterListener(this);

    sharblur->setAdjusterListener(this);

    if (showtooltip) {
        sensisha->set_tooltip_text(M("TP_LOCALLAB_SENSIS_TOOLTIP"));
    }

    sensisha->setAdjusterListener(this);

    inversshaConn = inverssha->signal_toggled().connect(sigc::mem_fun(*this, &LocallabSharp::inversshaChanged));

    sharFrame->set_label_align(0.025, 0.5);

    showmasksharMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasksharMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmasksharMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmasksharMethod->set_active(0);

    if (showtooltip) {
        showmasksharMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmasksharMethodConn = showmasksharMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabSharp::showmasksharMethodChanged));

    // Add Sharpening specific widgets to GUI
    pack_start(*sharcontrast);
    pack_start(*sharradius);
    pack_start(*sharamount);

    if (complexsoft < 2) {
        pack_start(*shardamping);
        pack_start(*shariter);
        pack_start(*sharblur);
    }

    pack_start(*sensisha);
    pack_start(*inverssha);
    ToolParamBlock* const sharfBox = Gtk::manage(new ToolParamBlock());
    sharfBox->pack_start(*showmasksharMethod);
    sharFrame->add(*sharfBox);
    pack_start(*sharFrame);
}

void LocallabSharp::resetMaskView()
{
    showmasksharMethodConn.block(true);
    showmasksharMethod->set_active(0);
    showmasksharMethodConn.block(false);
}

void LocallabSharp::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    sharMask = showmasksharMethod->get_active_row_number();
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
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visisharp);
        exp->setEnabled(pp->locallab.spots.at(index).expsharp);

        sharcontrast->setValue((double)pp->locallab.spots.at(index).sharcontrast);
        sharradius->setValue(pp->locallab.spots.at(index).sharradius);
        sharamount->setValue((double)pp->locallab.spots.at(index).sharamount);

        if (complexsoft < 2) {
            shardamping->setValue((double)pp->locallab.spots.at(index).shardamping);
            shariter->setValue((double)pp->locallab.spots.at(index).shariter);
            sharblur->setValue(pp->locallab.spots.at(index).sharblur);
        } else {
            shardamping->setValue(0.);
            shariter->setValue(30.);
            sharblur->setValue(0.2);
        }

        sensisha->setValue((double)pp->locallab.spots.at(index).sensisha);
        inverssha->set_active(pp->locallab.spots.at(index).inverssha);
    }

    // Enable all listeners
    enableListener();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSharp::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expsharp = exp->getEnabled();
        pp->locallab.spots.at(index).visisharp = exp->get_visible();

        pp->locallab.spots.at(index).sharcontrast = sharcontrast->getIntValue();
        pp->locallab.spots.at(index).sharradius = sharradius->getValue();
        pp->locallab.spots.at(index).sharamount = sharamount->getIntValue();
        pp->locallab.spots.at(index).shardamping = shardamping->getIntValue();
        pp->locallab.spots.at(index).shariter = shariter->getIntValue();
        pp->locallab.spots.at(index).sharblur = sharblur->getValue();
        pp->locallab.spots.at(index).sensisha = sensisha->getIntValue();
        pp->locallab.spots.at(index).inverssha = inverssha->get_active();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSharp::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

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

    // Local constrast specific widgets
    localcontMethod(Gtk::manage(new MyComboBoxText())),
    lcradius(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_RADIUS"), 10, 100, 1, 80))),
    lcamount(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_AMOUNT"), 0, 1.0, 0.01, 0))),
    lcdarkness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_DARKNESS"), 0, 3.0, 0.01, 1.0))),
    lclightness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_LIGHTNESS"), 0, 3.0, 0.01, 1.0))),
    LocalcurveEditorwav(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAV"))),
    wavshape(static_cast<FlatCurveEditor*>(LocalcurveEditorwav->addCurve(CT_Flat, "", nullptr, false, false))),
    levelwav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LEVELWAV"), 1, 9, 1, 4))),
    csThreshold(Gtk::manage(new ThresholdAdjuster(M("TP_LOCALLAB_CSTHRESHOLD"), 0, 9, 0, 0, 6, 6, 0, false))),
    expresidpyr(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::HBox())))),
    residcont(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCONT"), -100, 100, 1, 0))),
    residchro(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCHRO"), -100., 100., 1., 0.))),
    shresFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHRESFRA")))),
    residsha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDSHA"), -100., 100., 1., 0.))),
    residshathr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDSHATHR"), 0., 100., 1., 30.))),
    residhi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDHI"), -100., 100., 1., 0.))),
    residhithr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDHITHR"), 0., 100., 1., 70.))),
    sensilc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 30))),
    clariFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CLARIFRA")))),
    clarilres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARILRES"), -20., 100., 0.5, 0.))),
    claricres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARICRES"), -20., 100., 0.5, 0.))),
    clarisoft(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), -10.0, 1000.0, 0.5, 1.))),
    origlc(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ORIGLC")))),
    expcontrastpyr(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::HBox())))),
    gradwavFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADWAVFRA")))),
    wavgradl(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_GRALWFRA")))),
    strwav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -4.0, 4.0, 0.05, 0.))),
    angwav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    edgFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_EDGSHARPFRA")))),
    wavedg(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_EDGFRA")))),
    edgsBox(Gtk::manage(new ToolParamBlock())),
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
    labmNP(Gtk::manage(new Gtk::Label(M("TP_WAVELET_NPTYPE") + ":"))),
    localneiMethod(Gtk::manage(new MyComboBoxText())),
    blurlevelFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_BLURLEVELFRA")))),
    wavblur(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_BLURLEVELFRA")))),
    levelblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LEVELBLUR"), 0., 100., 0.5, 0.))),
    sigmabl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    chromablu(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMABLU"), 0.01, 5., 0.01, 1.))),
    LocalcurveEditorwavlev(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVLEV"))),
    wavshapelev(static_cast<FlatCurveEditor*>(LocalcurveEditorwavlev->addCurve(CT_Flat, "", nullptr, false, false))),
    residblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDBLUR"), 0., 100., 0.5, 0.))),
    blurlc(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_BLURLC")))),
    expcontrastpyr2(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::HBox())))),
    contFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CONTFRA")))),
    wavcont(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CONTFRA")))),
    sigma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    offset(Gtk::manage(new Adjuster(M("TP_LOCALLAB_OFFSETWAV"), 0.33, 1.66, 0.01, 1., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    chromalev(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMALEV"), 0.01, 5., 0.01, 1.))),
    LocalcurveEditorwavcon(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVCON"))),
    wavshapecon(static_cast<FlatCurveEditor*>(LocalcurveEditorwavcon->addCurve(CT_Flat, "", nullptr, false, false))),
    compreFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_COMPREFRA")))),
    wavcompre(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_COMPREFRA")))),
    LocalcurveEditorwavcompre(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVCOMPRE"))),
    wavshapecompre(static_cast<FlatCurveEditor*>(LocalcurveEditorwavcompre->addCurve(CT_Flat, "", nullptr, false, false))),
    sigmadr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    threswav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_THRESWAV"), 0.9, 2., 0.01, 1.4))),
    residcomp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCOMP"), -1., 1., 0.01, 0.))),
    compFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_COMPFRA")))),
    wavcomp(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_COMPFRA")))),
    fatdet(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATDETAIL"), -100., 300., 1., 0.))),
    fatanch(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATANCHOR"), 1., 100., 1., 50., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
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
    radmasklc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
    chromasklc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    mask2lcCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmasklcshape(static_cast<DiagonalCurveEditor*>(mask2lcCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
{
    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    const LocallabParams::LocallabSpot defSpot;

    // Parameter Local contrast specific widgets
    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPCONTRAST_TOOLTIP"));
    }

    localcontMethod->append(M("TP_LOCALLAB_LOCCONT"));
    localcontMethod->append(M("TP_LOCALLAB_WAVE"));
    localcontMethod->set_active(0);
    localcontMethodConn = localcontMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabContrast::localcontMethodChanged));

    if (complexsoft == 2) {
        lcradius->setLimits(20, 100, 1, 80);
    }

    lcradius->setAdjusterListener(this);

    lcamount->setAdjusterListener(this);

    lcdarkness->setAdjusterListener(this);

    lclightness->setAdjusterListener(this);

    LocalcurveEditorwav->setCurveListener(this);

    wavshape->setIdentityValue(0.);
    wavshape->setResetCurve(FlatCurveType(defSpot.locwavcurve.at(0)), defSpot.locwavcurve);
    wavshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    LocalcurveEditorwav->curveListComplete();

    if (showtooltip) {
        levelwav->set_tooltip_markup(M("TP_LOCALLAB_LEVELWAV_TOOLTIP"));
    }

    levelwav->setAdjusterListener(this);

    csThreshold->setAdjusterListener(this);

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

    if (showtooltip) {
        clariFrame->set_tooltip_markup(M("TP_LOCALLAB_CLARI_TOOLTIP"));
    }

    clariFrame->set_label_align(0.025, 0.5);

    clarilres->setAdjusterListener(this);

    claricres->setAdjusterListener(this);

    clarisoft->setLogScale(10, -10);
    clarisoft->setAdjusterListener(this);

    origlcConn = origlc->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::origlcChanged));

    Gtk::HBox* const LCTitleHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const LCLabel = Gtk::manage(new Gtk::Label());
    LCLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYR")) + Glib::ustring("</b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYRLAB")));
    LCLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    LCTitleHBox->pack_start(*LCLabel, Gtk::PACK_EXPAND_WIDGET, 0);
    expcontrastpyr->setLabel(LCTitleHBox);
    setExpandAlignProperties(expcontrastpyr, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    gradwavFrame->set_label_align(0.025, 0.5);

    wavgradlConn = wavgradl->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavgradlChanged));

    strwav->setAdjusterListener(this);

    angwav->setAdjusterListener(this);

    edgFrame->set_label_align(0.025, 0.5);

    if (showtooltip) {
        wavedg->set_tooltip_text(M("TP_LOCALLAB_WAVEEDG_TOOLTIP"));
    }

    wavedgConn = wavedg->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavedgChanged));

    strengthw->setAdjusterListener(this);

    sigmaed->setAdjusterListener(this);

    LocalcurveEditorwavedg->setCurveListener(this);

    wavshapeedg->setIdentityValue(0.);
    wavshapeedg->setResetCurve(FlatCurveType(defSpot.locedgwavcurve.at(0)), defSpot.locedgwavcurve);

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

    blurlevelFrame->set_label_align(0.025, 0.5);

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

    contFrame->set_label_align(0.025, 0.5);

    wavcontConn = wavcont->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavcontChanged));

    sigma->setAdjusterListener(this);

    offset->setAdjusterListener(this);

    chromalev->setAdjusterListener(this);

    LocalcurveEditorwavcon->setCurveListener(this);

    wavshapecon->setIdentityValue(0.);
    wavshapecon->setResetCurve(FlatCurveType(defSpot.locconwavcurve.at(0)), defSpot.locconwavcurve);

    LocalcurveEditorwavcon->curveListComplete();

    compreFrame->set_label_align(0.025, 0.5);

    wavcompreConn = wavcompre->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavcompreChanged));

    LocalcurveEditorwavcompre->setCurveListener(this);

    wavshapecompre->setIdentityValue(0.);
    wavshapecompre->setResetCurve(FlatCurveType(defSpot.loccomprewavcurve.at(0)), defSpot.loccomprewavcurve);

    if (showtooltip) {
        wavshapecompre->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    LocalcurveEditorwavcompre->curveListComplete();

    sigmadr->setAdjusterListener(this);

    threswav->setAdjusterListener(this);

    residcomp->setAdjusterListener(this);

    compFrame->set_label_align(0.025, 0.5);

    wavcompConn = wavcomp->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::wavcompChanged));

    if (showtooltip) {
        fatdet->set_tooltip_text(M("TP_LOCALLAB_COMPFRAME_TOOLTIP"));
    }

    fatdet->setAdjusterListener(this);

    fatanch->setAdjusterListener(this);

    LocalcurveEditorwavcomp->setCurveListener(this);

    wavshapecomp->setIdentityValue(0.);
    wavshapecomp->setResetCurve(FlatCurveType(defSpot.loccompwavcurve.at(0)), defSpot.loccompwavcurve);

    LocalcurveEditorwavcomp->curveListComplete();

    fatres->setAdjusterListener(this);

    if (showtooltip) {
        fftwlc->set_tooltip_text(M("TP_LOCALLAB_LC_FFTW_TOOLTIP"));
    }

    fftwlcConn = fftwlc->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::fftwlcChanged));

    if (showtooltip) {
        expmasklc->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    setExpandAlignProperties(expmasklc, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmasklcMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasklcMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmasklcMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmasklcMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmasklcMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmasklcMethod->set_active(0);

    if (showtooltip) {
        showmasklcMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmasklcMethodConn = showmasklcMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabContrast::showmasklcMethodChanged));

    enalcMaskConn = enalcMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::enalcMaskChanged));

    masklcCurveEditorG->setCurveListener(this);

    CCmasklcshape->setIdentityValue(0.);
    CCmasklcshape->setResetCurve(FlatCurveType(defSpot.CCmasklccurve.at(0)), defSpot.CCmasklccurve);

    if (showtooltip) {
        CCmasklcshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmasklcshape->setBottomBarColorProvider(this, 1);

    LLmasklcshape->setIdentityValue(0.);
    LLmasklcshape->setResetCurve(FlatCurveType(defSpot.LLmasklccurve.at(0)), defSpot.LLmasklccurve);

    if (showtooltip) {
        LLmasklcshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmasklcshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmasklcshape->setIdentityValue(0.);
    HHmasklcshape->setResetCurve(FlatCurveType(defSpot.HHmasklccurve.at(0)), defSpot.HHmasklccurve);

    if (showtooltip) {
        HHmasklcshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmasklcshape->setCurveColorProvider(this, 2);
    HHmasklcshape->setBottomBarColorProvider(this, 2);

    masklcCurveEditorG->curveListComplete();

    blendmasklc->setAdjusterListener(this);

    radmasklc->setLogScale(10, -10);
    radmasklc->setAdjusterListener(this);

    chromasklc->setAdjusterListener(this);

    mask2lcCurveEditorG->setCurveListener(this);

    Lmasklcshape->setResetCurve(DiagonalCurveType(defSpot.Lmasklccurve.at(0)), defSpot.Lmasklccurve);

    if (showtooltip) {
        Lmasklcshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    Lmasklcshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmasklcshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2lcCurveEditorG->curveListComplete();

    // Add Local contrast specific widgets to GUI
    if (complexsoft < 2) {
        pack_start(*localcontMethod);
    }
    shresFrame->set_label_align(0.025, 0.5);

    pack_start(*lcradius);
    pack_start(*lcamount);
    pack_start(*lcdarkness);
    pack_start(*lclightness);
    pack_start(*LocalcurveEditorwav, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    // pack_start(*levelwav);
    pack_start(*csThreshold);
    ToolParamBlock* const resiBox = Gtk::manage(new ToolParamBlock());
    resiBox->pack_start(*residcont);
    resiBox->pack_start(*residchro);
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
    gradwavFrame->set_label_widget(*wavgradl);
    ToolParamBlock* const gradwavBox = Gtk::manage(new ToolParamBlock());
    gradwavBox->pack_start(*strwav);
    gradwavBox->pack_start(*angwav);
    gradwavFrame->add(*gradwavBox);
    blurcontBox->pack_start(*gradwavFrame);
    edgFrame->set_label_widget(*wavedg);
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
    ctboxNP->pack_start(*labmNP, Gtk::PACK_SHRINK, 1);
    ctboxNP->pack_start(*localneiMethod);
    edgsBoxshow->pack_start(*ctboxNP);
    edgsBox->pack_start(*edgsBoxshow);
    edgFrame->add(*edgsBox);
    blurcontBox->pack_start(*edgFrame);
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
    Gtk::VBox* const contlevBox = Gtk::manage(new Gtk::VBox());
    contlevBox->set_spacing(2);
    contFrame->set_label_widget(*wavcont);
    contlevBox->pack_start(*sigma);
    contlevBox->pack_start(*offset);
    contlevBox->pack_start(*chromalev);
    contlevBox->pack_start(*LocalcurveEditorwavcon, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    contFrame->add(*contlevBox);
    blurcontBox2->pack_start(*contFrame);
    Gtk::VBox* const compreBox = Gtk::manage(new Gtk::VBox());
    compreBox->set_spacing(2);
    compreFrame->set_label_widget(*wavcompre);
    compreBox->pack_start(*LocalcurveEditorwavcompre, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    compreBox->pack_start(*sigmadr);
    compreBox->pack_start(*threswav);
    compreBox->pack_start(*residcomp);
    compreFrame->add(*compreBox);
    blurcontBox2->pack_start(*compreFrame);
    Gtk::VBox* const compBox = Gtk::manage(new Gtk::VBox());
    compBox->set_spacing(2);
    compFrame->set_label_widget(*wavcomp);
    compBox->pack_start(*fatdet);
    compBox->pack_start(*fatanch);
    compBox->pack_start(*LocalcurveEditorwavcomp, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    Gtk::HSeparator* const separatorcomp = Gtk::manage(new  Gtk::HSeparator());
    compBox->pack_start(*separatorcomp);
    compBox->pack_start(*fatres);
    compFrame->add(*compBox);
    // blurcontBox2->pack_start(*compFrame);
    expcontrastpyr2->add(*blurcontBox2, false);
    pack_start(*expcontrastpyr2);

    if (complexsoft < 2) {
        pack_start(*fftwlc);
    }

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

void LocallabContrast::resetMaskView()
{
    showmasklcMethodConn.block(true);
    showmasklcMethod->set_active(0);
    showmasklcMethodConn.block(false);
}

void LocallabContrast::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    lcMask = showmasklcMethod->get_active_row_number();
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
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visicontrast);
        exp->setEnabled(pp->locallab.spots.at(index).expcontrast);

        if (complexsoft < 2) {
            if (pp->locallab.spots.at(index).localcontMethod == "loc") {
                localcontMethod->set_active(0);
            } else if (pp->locallab.spots.at(index).localcontMethod == "wav") {
                localcontMethod->set_active(1);
            }
        } else {
            localcontMethod->set_active(1);
        }

        lcradius->setValue((double)pp->locallab.spots.at(index).lcradius);
        lcamount->setValue(pp->locallab.spots.at(index).lcamount);
        lcdarkness->setValue(pp->locallab.spots.at(index).lcdarkness);
        lclightness->setValue(pp->locallab.spots.at(index).lclightness);
        wavshape->setCurve(pp->locallab.spots.at(index).locwavcurve);
        levelwav->setValue((double)pp->locallab.spots.at(index).levelwav);
        csThreshold->setValue<int>(pp->locallab.spots.at(index).csthreshold);
        residcont->setValue(pp->locallab.spots.at(index).residcont);
        residchro->setValue(pp->locallab.spots.at(index).residchro);
        residsha->setValue(pp->locallab.spots.at(index).residsha);
        residshathr->setValue(pp->locallab.spots.at(index).residshathr);
        residhi->setValue(pp->locallab.spots.at(index).residhi);
        residhithr->setValue(pp->locallab.spots.at(index).residhithr);
        sensilc->setValue((double)pp->locallab.spots.at(index).sensilc);
        clarilres->setValue(pp->locallab.spots.at(index).clarilres);
        claricres->setValue(pp->locallab.spots.at(index).claricres);
        clarisoft->setValue(pp->locallab.spots.at(index).clarisoft);
        origlc->set_active(pp->locallab.spots.at(index).origlc);
        wavgradl->set_active(pp->locallab.spots.at(index).wavgradl);
        strwav->setValue(pp->locallab.spots.at(index).strwav);
        angwav->setValue(pp->locallab.spots.at(index).angwav);
        wavedg->set_active(pp->locallab.spots.at(index).wavedg);
        strengthw->setValue(pp->locallab.spots.at(index).strengthw);
        sigmaed->setValue(pp->locallab.spots.at(index).sigmaed);
        wavshapeedg->setCurve(pp->locallab.spots.at(index).locedgwavcurve);
        gradw->setValue(pp->locallab.spots.at(index).gradw);
        waveshow->set_active(pp->locallab.spots.at(index).waveshow);
        radiusw->setValue(pp->locallab.spots.at(index).radiusw);
        detailw->setValue(pp->locallab.spots.at(index).detailw);

        if (pp->locallab.spots.at(index).localedgMethod == "fir") {
            localedgMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).localedgMethod == "sec") {
            localedgMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).localedgMethod == "thr") {
            localedgMethod->set_active(2);
        }

        tloww->setValue(pp->locallab.spots.at(index).tloww);
        thigw->setValue(pp->locallab.spots.at(index).thigw);
        edgw->setValue(pp->locallab.spots.at(index).edgw);
        basew->setValue(pp->locallab.spots.at(index).basew);

        if (pp->locallab.spots.at(index).localneiMethod == "none") {
            localneiMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).localneiMethod == "low") {
            localneiMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).localneiMethod == "high") {
            localneiMethod->set_active(2);
        }

        wavblur->set_active(pp->locallab.spots.at(index).wavblur);
        levelblur->setValue(pp->locallab.spots.at(index).levelblur);
        sigmabl->setValue(pp->locallab.spots.at(index).sigmabl);
        chromablu->setValue(pp->locallab.spots.at(index).chromablu);
        wavshapelev->setCurve(pp->locallab.spots.at(index).loclevwavcurve);
        residblur->setValue(pp->locallab.spots.at(index).residblur);
        blurlc->set_active(pp->locallab.spots.at(index).blurlc);
        wavcont->set_active(pp->locallab.spots.at(index).wavcont);
        sigma->setValue(pp->locallab.spots.at(index).sigma);
        offset->setValue(pp->locallab.spots.at(index).offset);
        chromalev->setValue(pp->locallab.spots.at(index).chromalev);
        wavshapecon->setCurve(pp->locallab.spots.at(index).locconwavcurve);
        wavcompre->set_active(pp->locallab.spots.at(index).wavcompre);
        wavshapecompre->setCurve(pp->locallab.spots.at(index).loccomprewavcurve);
        sigmadr->setValue(pp->locallab.spots.at(index).sigmadr);
        threswav->setValue(pp->locallab.spots.at(index).threswav);
        residcomp->setValue(pp->locallab.spots.at(index).residcomp);
        wavcomp->set_active(pp->locallab.spots.at(index).wavcomp);
        fatdet->setValue(pp->locallab.spots.at(index).fatdet);
        fatanch->setValue(pp->locallab.spots.at(index).fatanch);
        wavshapecomp->setCurve(pp->locallab.spots.at(index).loccompwavcurve);
        fatres->setValue(pp->locallab.spots.at(index).fatres);

        if (complexsoft < 2) {
            fftwlc->set_active(pp->locallab.spots.at(index).fftwlc);
        } else {
            fftwlc->set_active(false);
        }

        enalcMask->set_active(pp->locallab.spots.at(index).enalcMask);
        CCmasklcshape->setCurve(pp->locallab.spots.at(index).CCmasklccurve);
        LLmasklcshape->setCurve(pp->locallab.spots.at(index).LLmasklccurve);
        HHmasklcshape->setCurve(pp->locallab.spots.at(index).HHmasklccurve);
        blendmasklc->setValue((double)pp->locallab.spots.at(index).blendmasklc);
        radmasklc->setValue(pp->locallab.spots.at(index).radmasklc);
        chromasklc->setValue(pp->locallab.spots.at(index).chromasklc);
        Lmasklcshape->setCurve(pp->locallab.spots.at(index).Lmasklccurve);
    }

    // Enable all listeners
    enableListener();

    // Update Local contrast GUI according to localcontMethod combobox value
    updateContrastGUI1();

    // Update Local contrast GUI according to waveshow button state
    updateContrastGUI2();

    // Update Local contrast GUI according to fftwlc button state
    updateContrastGUI3();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabContrast::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expcontrast = exp->getEnabled();
        pp->locallab.spots.at(index).visicontrast = exp->get_visible();

        if (localcontMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).localcontMethod = "loc";
        } else if (localcontMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).localcontMethod = "wav";
        }

        pp->locallab.spots.at(index).lcradius = lcradius->getIntValue();
        pp->locallab.spots.at(index).lcamount = lcamount->getValue();
        pp->locallab.spots.at(index).lcdarkness = lcdarkness->getValue();
        pp->locallab.spots.at(index).lclightness = lclightness->getValue();
        pp->locallab.spots.at(index).locwavcurve = wavshape->getCurve();
        pp->locallab.spots.at(index).levelwav = levelwav->getIntValue();
        pp->locallab.spots.at(index).csthreshold = csThreshold->getValue<int>();
        pp->locallab.spots.at(index).residcont = residcont->getValue();
        pp->locallab.spots.at(index).residchro = residchro->getValue();
        pp->locallab.spots.at(index).residsha = residsha->getValue();
        pp->locallab.spots.at(index).residshathr = residshathr->getValue();
        pp->locallab.spots.at(index).residhi = residhi->getValue();
        pp->locallab.spots.at(index).residhithr = residhithr->getValue();
        pp->locallab.spots.at(index).sensilc = sensilc->getIntValue();
        pp->locallab.spots.at(index).clarilres = clarilres->getValue();
        pp->locallab.spots.at(index).claricres = claricres->getValue();
        pp->locallab.spots.at(index).clarisoft = clarisoft->getValue();
        pp->locallab.spots.at(index).origlc = origlc->get_active();
        pp->locallab.spots.at(index).wavgradl = wavgradl->get_active();
        pp->locallab.spots.at(index).strwav = strwav->getValue();
        pp->locallab.spots.at(index).angwav = angwav->getValue();
        pp->locallab.spots.at(index).wavedg = wavedg->get_active();
        pp->locallab.spots.at(index).strengthw = strengthw->getValue();
        pp->locallab.spots.at(index).sigmaed = sigmaed->getValue();
        pp->locallab.spots.at(index).locedgwavcurve = wavshapeedg->getCurve();
        pp->locallab.spots.at(index).gradw = gradw->getValue();
        pp->locallab.spots.at(index).waveshow = waveshow->get_active();
        pp->locallab.spots.at(index).radiusw = radiusw->getValue();
        pp->locallab.spots.at(index).detailw = detailw->getValue();

        if (localedgMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).localedgMethod = "fir";
        } else if (localedgMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).localedgMethod = "sec";
        } else if (localedgMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).localedgMethod = "thr";
        }

        pp->locallab.spots.at(index).tloww = tloww->getValue();
        pp->locallab.spots.at(index).thigw = thigw->getValue();
        pp->locallab.spots.at(index).edgw = edgw->getValue();
        pp->locallab.spots.at(index).basew = basew->getValue();

        if (localneiMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).localneiMethod = "none";
        } else if (localneiMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).localneiMethod = "low";
        } else if (localneiMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).localneiMethod = "high";
        }

        pp->locallab.spots.at(index).wavblur = wavblur->get_active();
        pp->locallab.spots.at(index).levelblur = levelblur->getValue();
        pp->locallab.spots.at(index).sigmabl = sigmabl->getValue();
        pp->locallab.spots.at(index).chromablu = chromablu->getValue();
        pp->locallab.spots.at(index).loclevwavcurve = wavshapelev->getCurve();
        pp->locallab.spots.at(index).residblur = residblur->getValue();
        pp->locallab.spots.at(index).blurlc = blurlc->get_active();
        pp->locallab.spots.at(index).wavcont = wavcont->get_active();
        pp->locallab.spots.at(index).sigma = sigma->getValue();
        pp->locallab.spots.at(index).offset = offset->getValue();
        pp->locallab.spots.at(index).chromalev = chromalev->getValue();
        pp->locallab.spots.at(index).locconwavcurve = wavshapecon->getCurve();
        pp->locallab.spots.at(index).wavcompre = wavcompre->get_active();
        pp->locallab.spots.at(index).loccomprewavcurve = wavshapecompre->getCurve();
        pp->locallab.spots.at(index).sigmadr = sigmadr->getValue();
        pp->locallab.spots.at(index).threswav = threswav->getValue();
        pp->locallab.spots.at(index).residcomp = residcomp->getValue();
        pp->locallab.spots.at(index).wavcomp = wavcomp->get_active();
        pp->locallab.spots.at(index).fatdet = fatdet->getValue();
        pp->locallab.spots.at(index).fatanch = fatanch->getValue();
        pp->locallab.spots.at(index).loccompwavcurve = wavshapecomp->getCurve();
        pp->locallab.spots.at(index).fatres = fatres->getValue();
        pp->locallab.spots.at(index).fftwlc = fftwlc->get_active();
        pp->locallab.spots.at(index).enalcMask = enalcMask->get_active();
        pp->locallab.spots.at(index).CCmasklccurve = CCmasklcshape->getCurve();
        pp->locallab.spots.at(index).LLmasklccurve = LLmasklcshape->getCurve();
        pp->locallab.spots.at(index).HHmasklccurve = HHmasklcshape->getCurve();
        pp->locallab.spots.at(index).blendmasklc = blendmasklc->getIntValue();
        pp->locallab.spots.at(index).radmasklc = radmasklc->getValue();
        pp->locallab.spots.at(index).chromasklc = chromasklc->getValue();
        pp->locallab.spots.at(index).Lmasklccurve = Lmasklcshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabContrast::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster and threshold adjuster widgets
        lcradius->setDefault((double)defSpot.lcradius);
        lcamount->setDefault(defSpot.lcamount);
        lcdarkness->setDefault(defSpot.lcdarkness);
        lclightness->setDefault(defSpot.lclightness);
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
        fatdet->setDefault(defSpot.fatdet);
        fatanch->setDefault(defSpot.fatanch);
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

        if (a == fatdet) {
            if (listener) {
                listener->panelChanged(Evlocallabfatdet,
                                       fatdet->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == fatanch) {
            if (listener) {
                listener->panelChanged(Evlocallabfatanch,
                                       fatanch->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

void LocallabContrast::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmasklcshape->updateLocallabBackground(normChromar);
        LLmasklcshape->updateLocallabBackground(normLumar);
        HHmasklcshape->updateLocallabBackground(normHuer);

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
    // Update Local contrast GUI according to localcontMethod combobox value
    if (localcontMethod->get_active_row_number() == 0) {
        lcradius->show();
        lcamount->show();
        lcdarkness->show();
        lclightness->show();
        LocalcurveEditorwav->hide();
        levelwav->hide();
        csThreshold->hide();
        residcont->hide();
        residchro->hide();
        residsha->hide();
        residshathr->hide();
        residhi->hide();
        residhithr->hide();
        shresFrame->hide();
        clariFrame->hide();
        strwav->hide();
        angwav->hide();
        strengthw->hide();
        sigmaed->hide();
        LocalcurveEditorwavedg->hide();
        gradw->hide();
        radiusw->hide();
        detailw->hide();
        tloww->hide();
        thigw->hide();
        edgw->hide();
        basew->hide();
        levelblur->hide();
        sigmabl->hide();
        chromablu->hide();
        LocalcurveEditorwavlev->hide();
        residblur->hide();
        sigma->hide();
        offset->hide();
        chromalev->hide();
        LocalcurveEditorwavcon->hide();
        LocalcurveEditorwavcompre->hide();
        sigmadr->hide();
        threswav->hide();
        residcomp->hide();
        fatdet->hide();
        fatanch->hide();
        LocalcurveEditorwavcomp->hide();
        fatres->hide();
        fftwlc->show();
    } else if (localcontMethod->get_active_row_number() == 1) {
        lcradius->hide();
        lcamount->hide();
        lcdarkness->hide();
        lclightness->hide();
        LocalcurveEditorwav->show();
        levelwav->show();
        csThreshold->show();
        residcont->show();
        residchro->show();
        residsha->show();
        residshathr->show();
        residhi->show();
        residhithr->show();
        shresFrame->show();
        clariFrame->show();
        strwav->show();
        angwav->show();
        strengthw->show();
        sigmaed->show();
        LocalcurveEditorwavedg->show();
        gradw->show();
        radiusw->show();
        detailw->show();
        tloww->show();
        thigw->show();
        edgw->show();
        basew->show();
        levelblur->show();
        sigmabl->show();
        chromablu->show();
        LocalcurveEditorwavlev->show();
        residblur->show();
        sigma->show();
        offset->show();
        chromalev->show();
        LocalcurveEditorwavcon->show();
        LocalcurveEditorwavcompre->show();
        sigmadr->show();
        threswav->show();
        residcomp->show();
        fatdet->show();
        fatanch->show();
        LocalcurveEditorwavcomp->show();
        fatres->show();
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
    blurcbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURCBDL"), 0., 100., 0.1, 0.))),
    residFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RESID")))),
    clarityml(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARITYML"), 0.1, 100., 0.1, 0.1))),
    contresid(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRESID"), -100, 100, 1, 0))),
    softradiuscb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), -10.0, 1000.0, 0.5, 0.))),
    sensicb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSICB"), 0, 100, 1, 15))),
    expmaskcb(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWCB")))),
    showmaskcbMethod(Gtk::manage(new MyComboBoxText())),
    enacbMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    maskcbCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskcbshape(static_cast<FlatCurveEditor*>(maskcbCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskcbshape(static_cast<FlatCurveEditor*>(maskcbCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskcbshape(static_cast<FlatCurveEditor *>(maskcbCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
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
    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    const LocallabParams::LocallabSpot defSpot;

    // Parameter CBDL specific widgets
    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPCBDL_TOOLTIP"));
    }

    for (const auto adj : multiplier) {
        adj->setAdjusterListener(this);
    }

    if (showtooltip) {
        chromacbdl->set_tooltip_text(M("TP_LOCALLAB_CHROMACB_TOOLTIP"));
    }

    chromacbdl->setAdjusterListener(this);

    threshold->setAdjusterListener(this);

    blurcbdl->setAdjusterListener(this);

    residFrame->set_label_align(0.025, 0.5);

    clarityml->setAdjusterListener(this);

    contresid->setAdjusterListener(this);

    softradiuscb->setLogScale(10, -10);
    softradiuscb->setAdjusterListener(this);

    if (showtooltip) {
        sensicb->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    }

    sensicb->setAdjusterListener(this);

    if (showtooltip) {
        expmaskcb->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcbMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmaskcbMethod->set_active(0);

    if (showtooltip) {
        showmaskcbMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskcbMethodConn = showmaskcbMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabCBDL::showmaskcbMethodChanged));

    enacbMaskConn = enacbMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabCBDL::enacbMaskChanged));

    maskcbCurveEditorG->setCurveListener(this);

    CCmaskcbshape->setIdentityValue(0.);
    CCmaskcbshape->setResetCurve(FlatCurveType(defSpot.CCmaskcbcurve.at(0)), defSpot.CCmaskcbcurve);

    if (showtooltip) {
        CCmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmaskcbshape->setBottomBarColorProvider(this, 1);

    LLmaskcbshape->setIdentityValue(0.);
    LLmaskcbshape->setResetCurve(FlatCurveType(defSpot.LLmaskcbcurve.at(0)), defSpot.LLmaskcbcurve);

    if (showtooltip) {
        LLmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmaskcbshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskcbshape->setIdentityValue(0.);
    HHmaskcbshape->setResetCurve(FlatCurveType(defSpot.HHmaskcbcurve.at(0)), defSpot.HHmaskcbcurve);

    if (showtooltip) {
        HHmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmaskcbshape->setCurveColorProvider(this, 2);
    HHmaskcbshape->setBottomBarColorProvider(this, 2);

    maskcbCurveEditorG->curveListComplete();

    blendmaskcb->setAdjusterListener(this);

    if (showtooltip) {
        radmaskcb->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    radmaskcb->setLogScale(10, -10);
    radmaskcb->setAdjusterListener(this);

    if (showtooltip) {
        lapmaskcb->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    lapmaskcb->setAdjusterListener(this);

    chromaskcb->setAdjusterListener(this);

    gammaskcb->setAdjusterListener(this);

    slomaskcb->setAdjusterListener(this);

    mask2cbCurveEditorG->setCurveListener(this);

    Lmaskcbshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskcbcurve.at(0)), defSpot.Lmaskcbcurve);

    if (showtooltip) {
        Lmaskcbshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    Lmaskcbshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskcbshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2cbCurveEditorG->curveListComplete();

    lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumacontrastMinusPressed));

    lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumaneutralPressed));

    lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumacontrastPlusPressed));

    // Add CBDL specific widgets to GUI
    Gtk::HBox* buttonBox = Gtk::manage(new Gtk::HBox(true, 10));
    buttonBox->pack_start(*lumacontrastMinusButton);
    buttonBox->pack_start(*lumaneutralButton);
    buttonBox->pack_start(*lumacontrastPlusButton);
    pack_start(*buttonBox);

    for (const auto adj : multiplier) {
        pack_start(*adj);
    }

    Gtk::HSeparator* const separator = Gtk::manage(new  Gtk::HSeparator());
    pack_start(*separator, Gtk::PACK_SHRINK, 2);
    pack_start(*chromacbdl);
    pack_start(*threshold);
   // pack_start(*blurcbdl);
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

    if (complexsoft < 1) {
        maskcbBox->pack_start(*lapmaskcb, Gtk::PACK_SHRINK, 0);
    }

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

void LocallabCBDL::resetMaskView()
{
    showmaskcbMethodConn.block(true);
    showmaskcbMethod->set_active(0);
    showmaskcbMethodConn.block(false);
}

void LocallabCBDL::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    cbMask = showmaskcbMethod->get_active_row_number();
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
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visicbdl);
        exp->setEnabled(pp->locallab.spots.at(index).expcbdl);

        for (int i = 0; i < 6; i++) {
            multiplier[i]->setValue(pp->locallab.spots.at(index).mult[i]);
        }

        chromacbdl->setValue(pp->locallab.spots.at(index).chromacbdl);
        threshold->setValue(pp->locallab.spots.at(index).threshold);
        blurcbdl->setValue(pp->locallab.spots.at(index).blurcbdl);
        clarityml->setValue(pp->locallab.spots.at(index).clarityml);
        contresid->setValue((double)pp->locallab.spots.at(index).contresid);
        softradiuscb->setValue(pp->locallab.spots.at(index).softradiuscb);
        sensicb->setValue((double)pp->locallab.spots.at(index).sensicb);
        enacbMask->set_active(pp->locallab.spots.at(index).enacbMask);
        CCmaskcbshape->setCurve(pp->locallab.spots.at(index).CCmaskcbcurve);
        LLmaskcbshape->setCurve(pp->locallab.spots.at(index).LLmaskcbcurve);
        HHmaskcbshape->setCurve(pp->locallab.spots.at(index).HHmaskcbcurve);
        blendmaskcb->setValue((double)pp->locallab.spots.at(index).blendmaskcb);
        radmaskcb->setValue(pp->locallab.spots.at(index).radmaskcb);

        if (complexsoft == 0) {
            lapmaskcb->setValue(pp->locallab.spots.at(index).lapmaskcb);
        } else {
            lapmaskcb->setValue(0.);
        }

        chromaskcb->setValue(pp->locallab.spots.at(index).chromaskcb);
        gammaskcb->setValue(pp->locallab.spots.at(index).gammaskcb);
        slomaskcb->setValue(pp->locallab.spots.at(index).slomaskcb);
        Lmaskcbshape->setCurve(pp->locallab.spots.at(index).Lmaskcbcurve);
    }

    // Enable all listeners
    enableListener();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabCBDL::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expcbdl = exp->getEnabled();
        pp->locallab.spots.at(index).visicbdl = exp->get_visible();

        for (int i = 0; i < 6; i++) {
            pp->locallab.spots.at(index).mult[i] = multiplier[i]->getValue();
        }

        pp->locallab.spots.at(index).chromacbdl = chromacbdl->getValue();
        pp->locallab.spots.at(index).threshold = threshold->getValue();
        pp->locallab.spots.at(index).blurcbdl = blurcbdl->getValue();
        pp->locallab.spots.at(index).clarityml = clarityml->getValue();
        pp->locallab.spots.at(index).contresid = contresid->getIntValue();
        pp->locallab.spots.at(index).softradiuscb = softradiuscb->getValue();
        pp->locallab.spots.at(index).sensicb = sensicb->getIntValue();
        pp->locallab.spots.at(index).enacbMask = enacbMask->get_active();
        pp->locallab.spots.at(index).LLmaskcbcurve = LLmaskcbshape->getCurve();
        pp->locallab.spots.at(index).CCmaskcbcurve = CCmaskcbshape->getCurve();
        pp->locallab.spots.at(index).HHmaskcbcurve = HHmaskcbshape->getCurve();
        pp->locallab.spots.at(index).blendmaskcb = blendmaskcb->getIntValue();
        pp->locallab.spots.at(index).radmaskcb = radmaskcb->getValue();
        pp->locallab.spots.at(index).lapmaskcb = lapmaskcb->getValue();
        pp->locallab.spots.at(index).chromaskcb = chromaskcb->getValue();
        pp->locallab.spots.at(index).gammaskcb = gammaskcb->getValue();
        pp->locallab.spots.at(index).slomaskcb = slomaskcb->getValue();
        pp->locallab.spots.at(index).Lmaskcbcurve = Lmaskcbshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabCBDL::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        for (int i = 0; i < 6; i++) {
            multiplier[i]->setDefault(defSpot.mult[i]);
        }

        chromacbdl->setDefault(defSpot.chromacbdl);
        threshold->setDefault(defSpot.threshold);
        blurcbdl->setDefault(defSpot.blurcbdl);
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

        if (a == blurcbdl) {
            if (listener) {
                listener->panelChanged(EvLocallabblurcbdl,
                                       blurcbdl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

void LocallabCBDL::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskcbshape->updateLocallabBackground(normChromar);
        LLmaskcbshape->updateLocallabBackground(normLumar);
        HHmaskcbshape->updateLocallabBackground(normHuer);

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
    LocallabTool(this, M("TP_LOCALLAB_LOG_TOOLNAME"), M("TP_LOCALLAB_LOG"), false),

    // Log encoding specific widgets
    logPFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOGPFRA")))),
    autocompute(Gtk::manage(new Gtk::ToggleButton(M("TP_LOCALLAB_LOGAUTO")))),
    blackEv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLACK_EV"), -16.0, 0.0, 0.1, -5.0))),
    whiteEv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WHITE_EV"), 0.0, 32.0, 0.1, 10.0))),
    fullimage(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FULLIMAGE")))),
    logFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOGFRA")))),
    Autogray(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_AUTOGRAY")))),
    sourceGray(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOURCE_GRAY"), 1.0, 100.0, 0.1, 10.0))),
    targetGray(Gtk::manage(new Adjuster(M("TP_LOCALLAB_TARGET_GRAY"), 5.0, 80.0, 0.1, 18.0))),
    detail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DETAIL"), 0., 1., 0.01, 0.6))),
    baselog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BASELOG"), 1.3, 8., 0.05, 2., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    sensilog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSILOG"), 0, 100, 1, 50))),
    gradlogFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADLOGFRA")))),
    strlog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -2.0, 2.0, 0.05, 0.))),
    anglog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.)))
{
    // Parameter Log encoding specific widgets
    logPFrame->set_label_align(0.025, 0.5);

    autoconn = autocompute->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::autocomputeToggled));

    blackEv->setLogScale(2, -8);
    blackEv->setAdjusterListener(this);

    whiteEv->setLogScale(16, 0);
    whiteEv->setAdjusterListener(this);

    fullimageConn = fullimage->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::fullimageChanged));

    logFrame->set_label_align(0.025, 0.5);

    AutograyConn = Autogray->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::AutograyChanged));

    sourceGray->setAdjusterListener(this);

    targetGray->setAdjusterListener(this);

    detail->setAdjusterListener(this);

    baselog->setAdjusterListener(this);

    sensilog->setAdjusterListener(this);

    gradlogFrame->set_label_align(0.025, 0.5);

    strlog->setAdjusterListener(this);

    anglog->setAdjusterListener(this);

    // Add Log encoding specific widgets to HUI
    ToolParamBlock* const logPBox = Gtk::manage(new ToolParamBlock());
    logPBox->pack_start(*autocompute);
    logPBox->pack_start(*blackEv);
    logPBox->pack_start(*whiteEv);
    logPBox->pack_start(*fullimage);
    logPFrame->add(*logPBox);
    pack_start(*logPFrame);
    ToolParamBlock* const logFBox = Gtk::manage(new ToolParamBlock());
    logFBox->pack_start(*Autogray);
    logFBox->pack_start(*sourceGray);
    logFrame->add(*logFBox);
    pack_start(*logFrame);
    pack_start(*targetGray);
    pack_start(*detail);
    pack_start(*baselog);
    pack_start(*sensilog);
    ToolParamBlock* const gradlogBox = Gtk::manage(new ToolParamBlock());
    gradlogBox->pack_start(*strlog);
    gradlogBox->pack_start(*anglog);
    gradlogFrame->add(*gradlogBox);
    pack_start(*gradlogFrame);
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
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visilog);
        exp->setEnabled(pp->locallab.spots.at(index).explog);

        autocompute->set_active(pp->locallab.spots.at(index).autocompute);
        blackEv->setValue(pp->locallab.spots.at(index).blackEv);
        whiteEv->setValue(pp->locallab.spots.at(index).whiteEv);
        fullimage->set_active(pp->locallab.spots.at(index).fullimage);
        Autogray->set_active(pp->locallab.spots.at(index).Autogray);
        sourceGray->setValue(pp->locallab.spots.at(index).sourceGray);
        targetGray->setValue(pp->locallab.spots.at(index).targetGray);
        detail->setValue(pp->locallab.spots.at(index).detail);
        baselog->setValue(pp->locallab.spots.at(index).baselog);
        sensilog->setValue((double)pp->locallab.spots.at(index).sensilog);
        strlog->setValue(pp->locallab.spots.at(index).strlog);
        anglog->setValue(pp->locallab.spots.at(index).anglog);
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
        pp->locallab.spots.at(index).explog = exp->getEnabled();
        pp->locallab.spots.at(index).visilog = exp->get_visible();

        pp->locallab.spots.at(index).autocompute = autocompute->get_active();
        pp->locallab.spots.at(index).blackEv = blackEv->getValue();
        pp->locallab.spots.at(index).whiteEv = whiteEv->getValue();
        pp->locallab.spots.at(index).fullimage = fullimage->get_active();
        pp->locallab.spots.at(index).Autogray = Autogray->get_active();
        pp->locallab.spots.at(index).sourceGray = sourceGray->getValue();
        pp->locallab.spots.at(index).targetGray = targetGray->getValue();
        pp->locallab.spots.at(index).detail = detail->getValue();
        pp->locallab.spots.at(index).baselog = baselog->getValue();
        pp->locallab.spots.at(index).sensilog = sensilog->getIntValue();
        pp->locallab.spots.at(index).strlog = strlog->getValue();
        pp->locallab.spots.at(index).anglog = anglog->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabLog::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

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
        targetGray->set_sensitive(false);
    } else {
        blackEv->set_sensitive(true);
        whiteEv->set_sensitive(true);
        sourceGray->set_sensitive(true);
        targetGray->set_sensitive(true);
    }
}
