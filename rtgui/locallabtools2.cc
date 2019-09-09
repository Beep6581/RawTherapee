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

#define MINCHRO 0.
#define MAXCHRO 150
#define MAXCHROCC 100
#define MINNEIGH 4
#define MAXNEIGH 5000
#define CENTERNEIGH 200

using namespace rtengine;
using namespace procparams;

extern Options options;

static double retiSlider2neigh(double sval)
{
    // Slider range: 0 - 5000
    double neigh;

    if (sval <= 1000) {
        // Linear below center-temp
        neigh = MINNEIGH + (sval / 1000.0) * (CENTERNEIGH - MINNEIGH);
    } else {
        const double slope = (double)(CENTERNEIGH - MINNEIGH) / (MAXNEIGH - CENTERNEIGH);
        double x = (sval - 1000) / 1000; // x 0..1
        double y = x * slope + (1.0 - slope) * pow(x, 4.0);
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
        sval = ((neigh - MINNEIGH) / (CENTERNEIGH - MINNEIGH)) * 1000.0;
    } else {
        const double slope = (double)(CENTERNEIGH - MINNEIGH) / (MAXNEIGH - CENTERNEIGH);
        const double y = (neigh - CENTERNEIGH) / (MAXNEIGH - CENTERNEIGH);
        double x = pow(y, 0.25); // rough guess of x, will be a little lower
        double k = 0.1;
        bool add = true;

        // The y=f(x) function is a mess to invert, therefore we have this trial-refinement loop instead.
        // From tests, worst case is about 20 iterations, i.e. no problem
        for (;;) {
            double y1 = x * slope + (1.0 - slope) * pow(x, 4.0);

            if (1000 * fabs(y1 - y) < 0.1) {
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

        sval = 1000.0 + x * 1000.0;
    }

    if (sval < 0) {
        sval = 0;
    }

    if (sval > 5000) {
        sval = 5000;
    }

    return sval;
}

/* ==== LocallabTone ==== */
LocallabTone::LocallabTone():
    LocallabTool(this, M("TP_LOCALLAB_TONE_TOOLNAME"), M("TP_LOCALLAB_TM"), true, MaskNormal),

    // Tone mapping specific widgets
    amount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_AMOUNT"), 50., 100.0, 0.5, 95.))),
    stren(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STREN"), -0.5, 2.0, 0.01, 0.5))),
    equiltm(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_EQUIL")))),
    gamma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAM"), 0.4, 4.0, 0.11, 1.0))),
    satur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SATUR"), -100., 100., 0.1, 0.))), // By default satur = 0 ==> use Mantiuk value
    estop(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ESTOP"), 0.1, 4.0, 0.01, 0.5))),
    scaltm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALTM"), 0.1, 10.0, 0.01, 4.0))),
    rewei(Gtk::manage(new Adjuster(M("TP_LOCALLAB_REWEI"), 0, 3, 1, 0))),
    softradiustm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    sensitm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15)))
{
    const bool showtooltip = options.showtooltip;

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

    softradiustm->setAdjusterListener(this);

    sensitm->setAdjusterListener(this);

    if (showtooltip) {
        sensitm->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    }

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
}

void LocallabTone::getMaskView(int &colorMask, int &expMask, int &shMask, int &softMask, int &tmMask, int &retiMask, int &cbMask)
{
    tmMask = showMaskMethod->get_active_row_number();
}

void LocallabTone::disableListener()
{
    LocallabTool::disableListener();

    equiltmConn.block(true);
}

void LocallabTone::enableListener()
{
    LocallabTool::enableListener();

    equiltmConn.block(false);
}

void LocallabTone::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
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
        gamma->setValue(pp->locallab.spots.at(index).gamma);
        satur->setValue(pp->locallab.spots.at(index).satur);
        estop->setValue(pp->locallab.spots.at(index).estop);
        scaltm->setValue(pp->locallab.spots.at(index).scaltm);
        rewei->setValue(pp->locallab.spots.at(index).rewei);
        softradiustm->setValue(pp->locallab.spots.at(index).softradiustm);
        sensitm->setValue(pp->locallab.spots.at(index).sensitm);
        enaMask->set_active(pp->locallab.spots.at(index).enatmMask);
        CCMaskShape->setCurve(pp->locallab.spots.at(index).CCmasktmcurve);
        LLMaskShape->setCurve(pp->locallab.spots.at(index).LLmasktmcurve);
        HHMaskShape->setCurve(pp->locallab.spots.at(index).HHmasktmcurve);
        blendMask->setValue(pp->locallab.spots.at(index).blendmasktm);
        radMask->setValue(pp->locallab.spots.at(index).radmasktm);
        chroMask->setValue(pp->locallab.spots.at(index).chromasktm);
        gamMask->setValue(pp->locallab.spots.at(index).gammasktm);
        sloMask->setValue(pp->locallab.spots.at(index).slomasktm);
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
        pp->locallab.spots.at(index).enatmMask = enaMask->get_active();
        pp->locallab.spots.at(index).LLmasktmcurve = LLMaskShape->getCurve();
        pp->locallab.spots.at(index).CCmasktmcurve = CCMaskShape->getCurve();
        pp->locallab.spots.at(index).HHmasktmcurve = HHMaskShape->getCurve();
        pp->locallab.spots.at(index).blendmasktm = blendMask->getIntValue();
        pp->locallab.spots.at(index).radmasktm = radMask->getValue();
        pp->locallab.spots.at(index).chromasktm = chroMask->getValue();
        pp->locallab.spots.at(index).gammasktm = gamMask->getValue();
        pp->locallab.spots.at(index).slomasktm = sloMask->getValue();
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
        blendMask->setDefault((double)defSpot.blendmasktm);
        radMask->setDefault(defSpot.radmasktm);
        chroMask->setDefault(defSpot.chromasktm);
        gamMask->setDefault(defSpot.gammasktm);
        sloMask->setDefault(defSpot.slomasktm);
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

        if (a == blendMask) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmasktm,
                                       blendMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radMask) {
            if (listener) {
                listener->panelChanged(Evlocallabradmasktm,
                                       radMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chroMask) {
            if (listener) {
                listener->panelChanged(Evlocallabchromasktm,
                                       chroMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamMask) {
            if (listener) {
                listener->panelChanged(Evlocallabgammasktm,
                                       gamMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sloMask) {
            if (listener) {
                listener->panelChanged(Evlocallabslomasktm,
                                       sloMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabTone::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == CCMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmasktmshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmasktmshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmasktmshape,
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

void LocallabTone::enaMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaMask->get_active()) {
                listener->panelChanged(EvLocallabEnatmMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnatmMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabTone::showMaskMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
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

/* ==== LocallabRetinex ==== */
LocallabRetinex::LocallabRetinex():
    LocallabTool(this, M("TP_LOCALLAB_RET_TOOLNAME"), M("TP_LOCALLAB_RETI"), true, MaskWithTrMap),

    // Retinex specific widgets
    retinexMethod(Gtk::manage(new MyComboBoxText())),
    fftwreti(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FFTW")))),
    equilret(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_EQUIL")))),
    str(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STR"), 0., 100., 0.1, 0.0))),
    chrrt(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHRRT"), 0.0, 100.0, 0.1, 0.0))),
    neigh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NEIGH"), MINNEIGH, MAXNEIGH, 0.5, CENTERNEIGH, nullptr, nullptr, &retiSlider2neigh, &retiNeigh2Slider))),
    vart(Gtk::manage(new Adjuster(M("TP_LOCALLAB_VART"), 4.0, 500., 0.1, 70.))),
    scalereti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALERETI"), 1.0, 10.0, 1., 3.))),
    limd(Gtk::manage(new Adjuster(M("TP_LOCALLAB_THRESRETI"), 1.2, 100.0, 0.1, 8.))),
    darkness(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DARKRETI"), 0.01, 3.0, 0.01, 1.))),
    lightnessreti(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LIGHTRETI"), 0.01, 3.0, 0.01, 1.))),
    dehaz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DEHAZ"), 0, 100, 1, 0))),
    softradiusret(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GUIDFILTER"), 0.0, 100.0, 0.1, 0.))),
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIH"), 0, 100, 1, 30))),
    LocalcurveEditorgainT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAIN"))),
    inversret(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS"))))
{
    const bool showtooltip = options.showtooltip;

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

    fftwretiConn  = fftwreti->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::fftwretiChanged));

    equilretConn  = equilret->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::equilretChanged));

    str->setAdjusterListener(this);

    chrrt->setAdjusterListener(this);

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

    darkness->setAdjusterListener(this);

    lightnessreti->setAdjusterListener(this);

    dehaz->setAdjusterListener(this);

    if (showtooltip) {
        softradiusret->set_tooltip_text(M("TP_LOCALLAB_GUIDFILTER_TOOLTIP"));
    }

    softradiusret->setAdjusterListener(this);

    if (showtooltip) {
        sensih->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    }

    sensih->setAdjusterListener(this);

    LocalcurveEditorgainT->setCurveListener(this);

    cTgainshape = static_cast<FlatCurveEditor*>(LocalcurveEditorgainT->addCurve(CT_Flat, "", nullptr, false, false));
    cTgainshape->setIdentityValue(0.);
    cTgainshape->setResetCurve(FlatCurveType(LocallabParams::DEF_RET_CURVE.at(0)), LocallabParams::DEF_RET_CURVE);

    if (showtooltip) {
        cTgainshape->setTooltip(M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));
    }

    LocalcurveEditorgainT->curveListComplete();

    inversretConn  = inversret->signal_toggled().connect(sigc::mem_fun(*this, &LocallabRetinex::inversretChanged));

    pack_start(*retinexMethod);
    pack_start(*fftwreti);
    pack_start(*equilret);
    pack_start(*str);
    pack_start(*chrrt);
    pack_start(*neigh);
    pack_start(*vart);
    pack_start(*scalereti);
    pack_start(*limd);
    pack_start(*darkness);
    pack_start(*lightnessreti);
    pack_start(*dehaz);
    pack_start(*softradiusret);
    pack_start(*sensih);
    pack_start(*LocalcurveEditorgainT, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    pack_start(*inversret);
}

LocallabRetinex::~LocallabRetinex()
{
    delete LocalcurveEditorgainT;
}

void LocallabRetinex::getMaskView(int &colorMask, int &expMask, int &shMask, int &softMask, int &tmMask, int &retiMask, int &cbMask)
{
    retiMask = showMaskMethod->get_active_row_number();
}

void LocallabRetinex::disableListener()
{
    LocallabTool::disableListener();

    retinexMethodConn.block(true);
    fftwretiConn.block(true);
    equilretConn.block(true);
    inversretConn.block(true);
}

void LocallabRetinex::enableListener()
{
    LocallabTool::enableListener();

    retinexMethodConn.block(false);
    fftwretiConn.block(false);
    equilretConn.block(false);
    inversretConn.block(false);
}

void LocallabRetinex::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visireti);

        exp->setEnabled(pp->locallab.spots.at(index).expreti);

        if (pp->locallab.spots.at(index).retinexMethod == "low") {
            retinexMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).retinexMethod == "uni") {
            retinexMethod->set_active(1);
        } else {
            retinexMethod->set_active(2);
        }

        fftwreti->set_active(pp->locallab.spots.at(index).fftwreti);
        equilret->set_active(pp->locallab.spots.at(index).equilret);
        str->setValue(pp->locallab.spots.at(index).str);
        chrrt->setValue(pp->locallab.spots.at(index).chrrt);
        neigh->setValue(pp->locallab.spots.at(index).neigh);
        vart->setValue(pp->locallab.spots.at(index).vart);
        scalereti->setValue(pp->locallab.spots.at(index).scalereti);
        limd->setValue(pp->locallab.spots.at(index).limd);
        darkness->setValue(pp->locallab.spots.at(index).darkness);
        lightnessreti->setValue(pp->locallab.spots.at(index).lightnessreti);
        dehaz->setValue(pp->locallab.spots.at(index).dehaz);
        softradiusret->setValue(pp->locallab.spots.at(index).softradiusret);
        sensih->setValue(pp->locallab.spots.at(index).sensih);
        cTgainshape->setCurve(pp->locallab.spots.at(index).localTgaincurve);
        inversret->set_active(pp->locallab.spots.at(index).inversret);
        enaMask->set_active(pp->locallab.spots.at(index).enaretiMask);
        enaMaskTrMap->set_active(pp->locallab.spots.at(index).enaretiMasktmap);
        CCMaskShape->setCurve(pp->locallab.spots.at(index).CCmaskreticurve);
        LLMaskShape->setCurve(pp->locallab.spots.at(index).LLmaskreticurve);
        HHMaskShape->setCurve(pp->locallab.spots.at(index).HHmaskreticurve);
        blendMask->setValue(pp->locallab.spots.at(index).blendmaskreti);
        radMask->setValue(pp->locallab.spots.at(index).radmaskreti);
        chroMask->setValue(pp->locallab.spots.at(index).chromaskreti);
        gamMask->setValue(pp->locallab.spots.at(index).gammaskreti);
        sloMask->setValue(pp->locallab.spots.at(index).slomaskreti);
    }

    // Enable all listeners
    enableListener();

    // Update Retinex GUI according to scalereti adjuster value
    updateRetinexGUI();

    // Update Retinex GUI according to inversret button state
    updateRetinexGUI2();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabRetinex::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expreti = exp->getEnabled();

        pp->locallab.spots.at(index).visireti = exp->get_visible();

        if (retinexMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).retinexMethod = "low";
        } else if (retinexMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).retinexMethod = "uni";
        } else if (retinexMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).retinexMethod = "high";
        }

        pp->locallab.spots.at(index).fftwreti = fftwreti->get_active();
        pp->locallab.spots.at(index).equilret = equilret->get_active();
        pp->locallab.spots.at(index).str = str->getValue();
        pp->locallab.spots.at(index).chrrt = chrrt->getValue();
        pp->locallab.spots.at(index).neigh = neigh->getValue();
        pp->locallab.spots.at(index).vart = vart->getValue();
        pp->locallab.spots.at(index).scalereti = scalereti->getValue();
        pp->locallab.spots.at(index).limd = limd->getValue();
        pp->locallab.spots.at(index).darkness = darkness->getValue();
        pp->locallab.spots.at(index).lightnessreti = lightnessreti->getValue();
        pp->locallab.spots.at(index).dehaz = dehaz->getIntValue();
        pp->locallab.spots.at(index).softradiusret = softradiusret->getValue();
        pp->locallab.spots.at(index).sensih = sensih->getIntValue();
        pp->locallab.spots.at(index).localTgaincurve = cTgainshape->getCurve();
        pp->locallab.spots.at(index).inversret = inversret->get_active();
        pp->locallab.spots.at(index).enaretiMask = enaMask->get_active();
        pp->locallab.spots.at(index).enaretiMasktmap = enaMaskTrMap->get_active();
        pp->locallab.spots.at(index).LLmaskreticurve = LLMaskShape->getCurve();
        pp->locallab.spots.at(index).CCmaskreticurve = CCMaskShape->getCurve();
        pp->locallab.spots.at(index).HHmaskreticurve = HHMaskShape->getCurve();
        pp->locallab.spots.at(index).blendmaskreti = blendMask->getIntValue();
        pp->locallab.spots.at(index).radmaskreti = radMask->getValue();
        pp->locallab.spots.at(index).chromaskreti = chroMask->getValue();
        pp->locallab.spots.at(index).gammaskreti = gamMask->getValue();
        pp->locallab.spots.at(index).slomaskreti = sloMask->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabRetinex::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        str->setDefault(defSpot.str);
        chrrt->setDefault(defSpot.chrrt);
        neigh->setDefault(defSpot.neigh);
        vart->setDefault(defSpot.vart);
        scalereti->setDefault(defSpot.scalereti);
        limd->setDefault(defSpot.limd);
        darkness->setDefault(defSpot.darkness);
        lightnessreti->setDefault(defSpot.lightnessreti);
        dehaz->setDefault((double)defSpot.dehaz);
        softradiusret->setDefault(defSpot.softradiusret);
        sensih->setDefault((double)defSpot.sensih);
        blendMask->setDefault((double)defSpot.blendmaskreti);
        radMask->setDefault(defSpot.radmaskreti);
        chroMask->setDefault(defSpot.chromaskreti);
        gamMask->setDefault(defSpot.gammaskreti);
        sloMask->setDefault(defSpot.slomaskreti);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabRetinex::adjusterChanged(Adjuster* a, double newval)
{
    // Update Retinex GUI according to scalereti adjuster value
    if (a == scalereti) {
        updateRetinexGUI();
    }

    if (isLocActivated && exp->getEnabled()) {
        if (a == str) {
            if (listener) {
                listener->panelChanged(Evlocallabstr,
                                       str->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chrrt) {
            if (listener) {
                listener->panelChanged(Evlocallabchrrt,
                                       chrrt->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

        if (a == dehaz) {
            if (listener) {
                listener->panelChanged(Evlocallabdehaz,
                                       dehaz->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == softradiusret) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusret,
                                       softradiusret->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensih) {
            if (listener) {
                listener->panelChanged(Evlocallabsensih,
                                       sensih->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendMask) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskreti,
                                       blendMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radMask) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskreti,
                                       radMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chroMask) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskreti,
                                       chroMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamMask) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskreti,
                                       gamMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sloMask) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskreti,
                                       sloMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == cTgainshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCTgainCurve,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == CCMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskretishape,
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

void LocallabRetinex::enaMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaMask->get_active()) {
                listener->panelChanged(EvLocallabEnaretiMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaretiMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::enaMaskTrMapChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaMaskTrMap->get_active()) {
                listener->panelChanged(EvLocallabEnaretiMasktmap,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaretiMasktmap,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabRetinex::showMaskMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
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

void LocallabRetinex::updateRetinexGUI()
{
    // Update Retinex GUI according to scalereti adjuster value
    if (scalereti->getValue() == 1) {
        // limd->hide();
        LocalcurveEditorgainT->hide();
        retinexMethod->hide();
    } else {
        // limd->show();
        LocalcurveEditorgainT->show();
        retinexMethod->show();
    }
}

void LocallabRetinex::updateRetinexGUI2()
{
    // Update Retinex GUI according to inversret button state
    if (inversret->get_active()) {
        maskExp->hide();
    } else {
        maskExp->show();
    }
}

/* ==== LocallabSharp ==== */
LocallabSharp::LocallabSharp():
    LocallabTool(this, M("TP_LOCALLAB_SHARP_TOOLNAME"), M("TP_LOCALLAB_SHARP"), true, MaskNone),

    // Sharpening specific widgets
    sharcontrast(Gtk::manage(new Adjuster(M("TP_SHARPENING_CONTRAST"), 0, 200, 1, 20))),
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 0.4, 2.5, 0.01, 0.75))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 100))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 0))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sharblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARBLUR"), 0.2, 2.0, 0.05, 0.2))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    inverssha(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS"))))
{
    const bool showtooltip = options.showtooltip;

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

    inversshaConn  = inverssha->signal_toggled().connect(sigc::mem_fun(*this, &LocallabSharp::inversshaChanged));

    pack_start(*sharcontrast);
    pack_start(*sharradius);
    pack_start(*sharamount);
    pack_start(*shardamping);
    pack_start(*shariter);
    pack_start(*sharblur);
    pack_start(*sensisha);
    pack_start(*inverssha);
}

void LocallabSharp::disableListener()
{
    LocallabTool::disableListener();

    inversshaConn.block(true);
}

void LocallabSharp::enableListener()
{
    LocallabTool::enableListener();

    inversshaConn.block(false);
}

void LocallabSharp::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visisharp);

        exp->setEnabled(pp->locallab.spots.at(index).expsharp);
        sharcontrast->setValue(pp->locallab.spots.at(index).sharcontrast);
        sharradius->setValue(pp->locallab.spots.at(index).sharradius);
        sharamount->setValue(pp->locallab.spots.at(index).sharamount);
        shardamping->setValue(pp->locallab.spots.at(index).shardamping);
        shariter->setValue(pp->locallab.spots.at(index).shariter);
        sharblur->setValue(pp->locallab.spots.at(index).sharblur);
        sensisha->setValue(pp->locallab.spots.at(index).sensisha);
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

/* ==== LocallabContrast ==== */
LocallabContrast::LocallabContrast():
    LocallabTool(this, M("TP_LOCALLAB_LC_TOOLNAME"), M("TP_LOCALLAB_LOC_CONTRAST"), false, MaskNone),

    // Local constrast specific widgets
    localcontMethod(Gtk::manage(new MyComboBoxText())),
    lcradius(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_RADIUS"), 20, 400, 1, 80))),
    lcamount(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_AMOUNT"), 0, 1.0, 0.01, 0))),
    lcdarkness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_DARKNESS"), 0, 3.0, 0.01, 1.0))),
    lclightness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_LIGHTNESS"), 0, 3.0, 0.01, 1.0))),
    LocalcurveEditorwav(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAV"))),
    levelwav(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LEVELWAV"), 3, 9, 1, 4))),
    residcont(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCONT"), -100, 100, 1, 0))),
    sensilc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    fftwlc(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FFTW"))))
{
    const bool showtooltip = options.showtooltip;

    localcontMethod->append(M("TP_LOCALLAB_LOCCONT"));
    localcontMethod->append(M("TP_LOCALLAB_WAVE"));
    localcontMethod->set_active(0);

    if (showtooltip) {
        // localcontMethod->set_tooltip_markup(M("TP_LOCALLAB_LOCMETHOD_TOOLTIP"));
    }

    localcontMethodConn = localcontMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabContrast::localcontMethodChanged));

    lcradius->setAdjusterListener(this);

    lcamount->setAdjusterListener(this);

    lcdarkness->setAdjusterListener(this);

    lclightness->setAdjusterListener(this);

    LocalcurveEditorwav->setCurveListener(this);

    wavshape = static_cast<FlatCurveEditor*>(LocalcurveEditorwav->addCurve(CT_Flat, "", nullptr, false, false));
    wavshape->setIdentityValue(0.);
    wavshape->setResetCurve(FlatCurveType(LocallabParams::DEF_LC_CURVE.at(0)), LocallabParams::DEF_LC_CURVE);

    if (showtooltip) {
        wavshape->setTooltip(M("TP_RETINEX_WAV_TOOLTIP"));
    }

    LocalcurveEditorwav->curveListComplete();

    if (showtooltip) {
        levelwav->set_tooltip_markup(M("TP_LOCALLAB_LEVELWAV_TOOLTIP"));
    }

    levelwav->setAdjusterListener(this);

    residcont->setAdjusterListener(this);

    sensilc->setAdjusterListener(this);

    if (showtooltip) {
        fftwlc->set_tooltip_text(M("TP_LOCALLAB_LC_FFTW_TOOLTIP"));
    }

    fftwlcConn  = fftwlc->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::fftwlcChanged));

    pack_start(*localcontMethod);
    pack_start(*lcradius);
    pack_start(*lcamount);
    pack_start(*lcdarkness);
    pack_start(*lclightness);
    pack_start(*LocalcurveEditorwav, Gtk::PACK_SHRINK, 4);
    pack_start(*levelwav);
    pack_start(*residcont);
    pack_start(*sensilc);
    pack_start(*fftwlc);
}

LocallabContrast::~LocallabContrast()
{
    delete LocalcurveEditorwav;
}

void LocallabContrast::disableListener()
{
    LocallabTool::disableListener();

    localcontMethodConn.block(true);
    fftwlcConn.block(true);
}

void LocallabContrast::enableListener()
{
    LocallabTool::enableListener();

    localcontMethodConn.block(false);
    fftwlcConn.block(false);
}

void LocallabContrast::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visicontrast);

        exp->setEnabled(pp->locallab.spots.at(index).expcontrast);

        if (pp->locallab.spots.at(index).localcontMethod == "loc") {
            localcontMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).localcontMethod == "wav") {
            localcontMethod->set_active(1);
        }

        lcradius->setValue(pp->locallab.spots.at(index).lcradius);
        lcamount->setValue(pp->locallab.spots.at(index).lcamount);
        lcdarkness->setValue(pp->locallab.spots.at(index).lcdarkness);
        lclightness->setValue(pp->locallab.spots.at(index).lclightness);
        wavshape->setCurve(pp->locallab.spots.at(index).locwavcurve);
        levelwav->setValue(pp->locallab.spots.at(index).levelwav);
        residcont->setValue(pp->locallab.spots.at(index).residcont);
        sensilc->setValue(pp->locallab.spots.at(index).sensilc);
        fftwlc->set_active(pp->locallab.spots.at(index).fftwlc);
    }

    // Enable all listeners
    enableListener();

    // Update Local constrast GUI according to localcontMethod combobox value
    updateContrastGUI();

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
        pp->locallab.spots.at(index).residcont = residcont->getValue();
        pp->locallab.spots.at(index).sensilc = sensilc->getIntValue();
        pp->locallab.spots.at(index).fftwlc = fftwlc->get_active();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabContrast::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        lcradius->setDefault((double)defSpot.lcradius);
        lcamount->setDefault(defSpot.lcamount);
        lcdarkness->setDefault(defSpot.lcdarkness);
        lclightness->setDefault(defSpot.lclightness);
        levelwav->setDefault(defSpot.levelwav);
        residcont->setDefault(defSpot.residcont);
        sensilc->setDefault((double)defSpot.sensilc);
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

        if (a == sensilc) {
            if (listener) {
                listener->panelChanged(Evlocallabsensilc,
                                       sensilc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
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

void LocallabContrast::localcontMethodChanged()
{
    // Update Local constrast GUI according to localcontMethod combobox value
    updateContrastGUI();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallablocalcontMethod,
                                   localcontMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabContrast::fftwlcChanged()
{
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

void LocallabContrast::updateContrastGUI()
{
    // Update Local constrast GUI according to localcontMethod combobox value
    if (localcontMethod->get_active_row_number() == 0) {
        lcradius->show();
        lcamount->show();
        lcdarkness->show();
        lclightness->show();
        LocalcurveEditorwav->hide();
        levelwav->hide();
        residcont->hide();
        fftwlc->show();
    } else if (localcontMethod->get_active_row_number() == 1) {
        lcradius->hide();
        lcamount->hide();
        lcdarkness->hide();
        lclightness->hide();
        LocalcurveEditorwav->show();
        levelwav->show();
        residcont->show();
        fftwlc->hide();
    }
}

/* ==== LocallabCBDL ==== */
LocallabCBDL::LocallabCBDL():
    LocallabTool(this, M("TP_LOCALLAB_CBDL_TOOLNAME"), M("TP_LOCALLAB_CBDL"), true, MaskNormal),

    // CBDL specific widgets
    chromacbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMACBDL"), 0., 1.5, 0.01, 0.))),
    threshold(Gtk::manage(new Adjuster(M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 1., 0.01, 0.2))),
    blurcbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURCBDL"), 0., 100., 0.1, 0.))),
    clarityml(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARITYML"), 0.1, 100., 0.1, 0.1))),
    contresid(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRESID"), -100, 100, 1, 0))),
    softradiuscb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    sensicb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSICB"), 0, 100, 1, 15))),

    lumacontrastMinusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")))),
    lumaneutralButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")))),
    lumacontrastPlusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS"))))
{
    const bool showtooltip = options.showtooltip;

    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPCBDL_TOOLTIP"));
    }

    for (int i = 0; i < 6; i++) {
        Glib::ustring ss;
        ss = Glib::ustring::format(i);

        if (i == 0) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMAFINEST"));
        } else if (i == 5) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        }

        multiplier[i] = Gtk::manage(new Adjuster(std::move(ss), 0.0, 4.0, 0.01, 1.0));
        multiplier[i]->setAdjusterListener(this);
    }

    if (showtooltip) {
        chromacbdl->set_tooltip_text(M("TP_LOCALLAB_CHROMACB_TOOLTIP"));
    }

    chromacbdl->setAdjusterListener(this);

    threshold->setAdjusterListener(this);

    blurcbdl->setAdjusterListener(this);

    clarityml->setAdjusterListener(this);

    contresid->setAdjusterListener(this);

    softradiuscb->setAdjusterListener(this);

    if (showtooltip) {
        sensicb->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    }

    sensicb->setAdjusterListener(this);

    lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumacontrastMinusPressed));

    lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumaneutralPressed));

    lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect(sigc::mem_fun(*this, &LocallabCBDL::lumacontrastPlusPressed));

    Gtk::HBox* buttonBox = Gtk::manage(new Gtk::HBox(true, 10));
    buttonBox->pack_start(*lumacontrastMinusButton);
    buttonBox->pack_start(*lumaneutralButton);
    buttonBox->pack_start(*lumacontrastPlusButton);
    pack_start(*buttonBox);

    for (int i = 0; i < 6; i++) {
        pack_start(*multiplier[i]);
    }

    Gtk::HSeparator* const separator = Gtk::manage(new  Gtk::HSeparator());
    pack_start(*separator, Gtk::PACK_SHRINK, 2);
    pack_start(*chromacbdl);
    pack_start(*threshold);
    pack_start(*blurcbdl);
    Gtk::Frame* const residFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RESID")));
    residFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const residBox = Gtk::manage(new ToolParamBlock());
    residBox->pack_start(*clarityml);
    residBox->pack_start(*contresid);
    residFrame->add(*residBox);
    pack_start(*residFrame);
    pack_start(*softradiuscb);
    pack_start(*sensicb);
}

void LocallabCBDL::getMaskView(int &colorMask, int &expMask, int &shMask, int &softMask, int &tmMask, int &retiMask, int &cbMask)
{
    cbMask = showMaskMethod->get_active_row_number();
}

void LocallabCBDL::disableListener()
{
    LocallabTool::disableListener();

    lumacontrastMinusPressedConn.block(true);
    lumaneutralPressedConn.block(true);
    lumacontrastPlusPressedConn.block(true);
}

void LocallabCBDL::enableListener()
{
    LocallabTool::enableListener();

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
        contresid->setValue(pp->locallab.spots.at(index).contresid);
        softradiuscb->setValue(pp->locallab.spots.at(index).softradiuscb);
        sensicb->setValue(pp->locallab.spots.at(index).sensicb);
        enaMask->set_active(pp->locallab.spots.at(index).enacbMask);
        CCMaskShape->setCurve(pp->locallab.spots.at(index).CCmaskcbcurve);
        LLMaskShape->setCurve(pp->locallab.spots.at(index).LLmaskcbcurve);
        HHMaskShape->setCurve(pp->locallab.spots.at(index).HHmaskcbcurve);
        blendMask->setValue(pp->locallab.spots.at(index).blendmaskcb);
        radMask->setValue(pp->locallab.spots.at(index).radmaskcb);
        chroMask->setValue(pp->locallab.spots.at(index).chromaskcb);
        gamMask->setValue(pp->locallab.spots.at(index).gammaskcb);
        sloMask->setValue(pp->locallab.spots.at(index).slomaskcb);
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
        pp->locallab.spots.at(index).enacbMask = enaMask->get_active();
        pp->locallab.spots.at(index).LLmaskcbcurve = LLMaskShape->getCurve();
        pp->locallab.spots.at(index).CCmaskcbcurve = CCMaskShape->getCurve();
        pp->locallab.spots.at(index).HHmaskcbcurve = HHMaskShape->getCurve();
        pp->locallab.spots.at(index).blendmaskcb = blendMask->getIntValue();
        pp->locallab.spots.at(index).radmaskcb = radMask->getValue();
        pp->locallab.spots.at(index).chromaskcb = chroMask->getValue();
        pp->locallab.spots.at(index).gammaskcb = gamMask->getValue();
        pp->locallab.spots.at(index).slomaskcb = sloMask->getValue();
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
        blendMask->setDefault((double)defSpot.blendmaskcb);
        radMask->setDefault(defSpot.radmaskcb);
        chroMask->setDefault(defSpot.chromaskcb);
        gamMask->setDefault(defSpot.gammaskcb);
        sloMask->setDefault(defSpot.slomaskcb);
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

        if (a == blendMask) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskcb,
                                       blendMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radMask) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskcb,
                                       radMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chroMask) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskcb,
                                       chroMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamMask) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskcb,
                                       gamMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sloMask) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskcb,
                                       sloMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabCBDL::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == CCMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskcbshape,
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

void LocallabCBDL::enaMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaMask->get_active()) {
                listener->panelChanged(EvLocallabEnacbMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnacbMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabCBDL::showMaskMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
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

/* ==== LocallabDenoise ==== */
LocallabDenoise::LocallabDenoise():
    LocallabTool(this, M("TP_LOCALLAB_DEN_TOOLNAME"), M("TP_LOCALLAB_DENOIS"), true, MaskNone),

    // Denoise specific widgets
    noiselumf0(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINEZERO"), MINCHRO, MAXCHRO, 1, 0))),
    noiselumf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINE"), MINCHRO, MAXCHRO, 1, 0))),
    noiselumf2(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINETWO"), MINCHRO, MAXCHRO, 1, 0))),
    noiselumc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMCOARSE"), MINCHRO, MAXCHROCC, 1, 0))),
    noiselumdetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMDETAIL"), 0, 100, 1, 0))),
    noiselequal(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELEQUAL"), -2, 10, 1, 7, Gtk::manage(new RTImage("circle-white-small.png")), Gtk::manage(new RTImage("circle-black-small.png"))))),
    noisechrof(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROFINE"), MINCHRO, MAXCHRO, 1, 0))),
    noisechroc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROCOARSE"), MINCHRO, MAXCHROCC, 1, 0))),
    noisechrodetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHRODETAIL"), 0, 100, 1, 0))),
    adjblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ADJ"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-red-small.png"))))),
    bilateral(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BILATERAL"), 0, 100, 1, 0))),
    sensiden(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIDEN"), 0, 100, 1, 20)))
{
    const bool showtooltip = options.showtooltip;

    noiselumf0->setAdjusterListener(this);

    noiselumf->setAdjusterListener(this);

    noiselumf2->setAdjusterListener(this);

    if (showtooltip) {
        noiselumc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    }

    noiselumc->setAdjusterListener(this);

    noiselumdetail->setAdjusterListener(this);

    noiselequal->setAdjusterListener(this);

    noisechrof->setAdjusterListener(this);

    if (showtooltip) {
        noisechroc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    }

    noisechroc->setAdjusterListener(this);

    noisechrodetail->setAdjusterListener(this);

    adjblur->setAdjusterListener(this);

    bilateral->setAdjusterListener(this);

    sensiden->setAdjusterListener(this);

    Gtk::Frame* const wavFrame = Gtk::manage(new Gtk::Frame());
    ToolParamBlock* const wavBox = Gtk::manage(new ToolParamBlock());
    wavBox->pack_start(*noiselumf0);
    wavBox->pack_start(*noiselumf);
    wavBox->pack_start(*noiselumf2);
    wavBox->pack_start(*noiselumc);
    wavBox->pack_start(*noiselumdetail);
    wavBox->pack_start(*noiselequal);
    wavBox->pack_start(*noisechrof);
    wavBox->pack_start(*noisechroc);
    // wavBox->pack_start(*noisechrodetail); // Uncomment this line to use the noisechrodetail adjuster
    wavBox->pack_start(*adjblur);
    wavFrame->add(*wavBox);
    pack_start(*wavFrame);
    pack_start(*bilateral);
    pack_start(*sensiden);
}

void LocallabDenoise::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visidenoi);

        exp->setEnabled(pp->locallab.spots.at(index).expdenoi);
        noiselumf0->setValue(pp->locallab.spots.at(index).noiselumf0);
        noiselumf->setValue(pp->locallab.spots.at(index).noiselumf);
        noiselumf2->setValue(pp->locallab.spots.at(index).noiselumf2);
        noiselumc->setValue(pp->locallab.spots.at(index).noiselumc);
        noiselumdetail->setValue(pp->locallab.spots.at(index).noiselumdetail);
        noiselequal->setValue(pp->locallab.spots.at(index).noiselequal);
        noisechrof->setValue(pp->locallab.spots.at(index).noisechrof);
        noisechroc->setValue(pp->locallab.spots.at(index).noisechroc);
        noisechrodetail->setValue(pp->locallab.spots.at(index).noisechrodetail);
        adjblur->setValue(pp->locallab.spots.at(index).adjblur);
        bilateral->setValue(pp->locallab.spots.at(index).bilateral);
        sensiden->setValue(pp->locallab.spots.at(index).sensiden);
    }

    // Enable all listeners
    enableListener();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabDenoise::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expdenoi = exp->getEnabled();

        pp->locallab.spots.at(index).visidenoi = exp->get_visible();

        pp->locallab.spots.at(index).noiselumf0 = noiselumf0->getIntValue();
        pp->locallab.spots.at(index).noiselumf = noiselumf->getIntValue();
        pp->locallab.spots.at(index).noiselumf2 = noiselumf2->getIntValue();
        pp->locallab.spots.at(index).noiselumc = noiselumc->getIntValue();
        pp->locallab.spots.at(index).noiselumdetail = noiselumdetail->getIntValue();
        pp->locallab.spots.at(index).noiselequal = noiselequal->getIntValue();
        pp->locallab.spots.at(index).noisechrof = noisechrof->getIntValue();
        pp->locallab.spots.at(index).noisechroc = noisechroc->getIntValue();
        pp->locallab.spots.at(index).noisechrodetail = noisechrodetail->getIntValue();
        pp->locallab.spots.at(index).adjblur = adjblur->getIntValue();
        pp->locallab.spots.at(index).bilateral = bilateral->getIntValue();
        pp->locallab.spots.at(index).sensiden = sensiden->getIntValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabDenoise::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        noiselumf0->setDefault((double)defSpot.noiselumf0);
        noiselumf->setDefault((double)defSpot.noiselumf);
        noiselumf2->setDefault((double)defSpot.noiselumf2);
        noiselumc->setDefault((double)defSpot.noiselumc);
        noiselumdetail->setDefault((double)defSpot.noiselumdetail);
        noiselequal->setDefault((double)defSpot.noiselequal);
        noisechrof->setDefault((double)defSpot.noisechrof);
        noisechroc->setDefault((double)defSpot.noisechroc);
        noisechrodetail->setDefault((double)defSpot.noisechrodetail);
        adjblur->setDefault((double)defSpot.adjblur);
        bilateral->setDefault((double)defSpot.bilateral);
        sensiden->setDefault((double)defSpot.sensiden);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabDenoise::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == noiselumf0) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf0,
                                       noiselumf0->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumf) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf,
                                       noiselumf->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumf2) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf2,
                                       noiselumf2->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumc) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumc,
                                       noiselumc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumdetail) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumdetail,
                                       noiselumdetail->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselequal) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselequal,
                                       noiselequal->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noisechrof) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechrof,
                                       noisechrof->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noisechroc) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechroc,
                                       noisechroc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noisechrodetail) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechrodetail,
                                       noisechrodetail->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == adjblur) {
            if (listener) {
                listener->panelChanged(Evlocallabadjblur,
                                       adjblur->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == bilateral) {
            if (listener) {
                listener->panelChanged(Evlocallabbilateral,
                                       bilateral->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensiden) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiden,
                                       sensiden->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabDenoise::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenadenoi,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenadenoi,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}
