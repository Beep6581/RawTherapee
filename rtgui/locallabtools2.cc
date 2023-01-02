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
#include "../rtengine/color.h"

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
    repartm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGREPART"), 1.0, 100.0, 1., 100.0))),
    amount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_AMOUNT"), 50., 100.0, 0.5, 95.))),
    stren(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STREN"), -0.5, 2.0, 0.01, 0.5))),
    equiltm(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_EQUIL")))),
    gamma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAM"), 0.4, 4.0, 0.11, 1.0))),
    satur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SATUR"), -100., 100., 0.1, 0.))), // By default satur = 0 ==> use Mantiuk value
    estop(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ESTOP"), 0.1, 4., 0.01, 1.4))),
    scaltm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALTM"), 0.1, 10.0, 0.01, 1.0))),
    rewei(Gtk::manage(new Adjuster(M("TP_LOCALLAB_REWEI"), 0, 3, 1, 0))),
    softradiustm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),//unused here, but used for normalize_mean_dt 
    sensitm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
    exprecovt(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_DENOI2_EXP")))),
    maskusablet(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUSABLE")))),
    maskunusablet(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUNUSABLE")))),
    recothrest(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKRECOTHRES"), 0., 2., 0.01, 1.))),
    lowthrest(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHRLOW"), 1., 80., 0.5, 12.))),
    higthrest(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHR"), 20., 99., 0.5, 85.))),
    decayt(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKDDECAY"), 0.5, 4., 0.1, 2.))),
    expmasktm(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWT")))),
    showmasktmMethod(Gtk::manage(new MyComboBoxText())),
    enatmMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    enatmMaskaft(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_AFTER_MASK")))),
 //   masktmCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    masktmCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    CCmasktmshape(static_cast<FlatCurveEditor*>(masktmCurveEditorG->addCurve(CT_Flat, "C", nullptr, false, false))),
    LLmasktmshape(static_cast<FlatCurveEditor*>(masktmCurveEditorG->addCurve(CT_Flat, "L", nullptr, false, false))),
    HHmasktmshape(static_cast<FlatCurveEditor *>(masktmCurveEditorG->addCurve(CT_Flat, "LC(h)", nullptr, false, true))),
    blendmasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    lapmasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    radmasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.05, 5.0, 0.01, 1.))),
    slomasktm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2tmCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmasktmshape(static_cast<DiagonalCurveEditor*>(mask2tmCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
    const LocallabParams::LocallabSpot defSpot;

    // Parameter Tone Mapping specific widgets
    amount->setAdjusterListener(this);

    repartm->setAdjusterListener(this);

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

    recothrest->setAdjusterListener(this);
    lowthrest->setAdjusterListener(this);
    higthrest->setAdjusterListener(this);
    decayt->setAdjusterListener(this);
    setExpandAlignProperties(exprecovt, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

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
    Gtk::Separator* const separatortm = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));

    // Add Tone Mapping specific widgets to GUI
    // pack_start(*amount); // To use if we change transit_shapedetect parameters
    pack_start(*sensitm);
    pack_start(*repartm);
    pack_start(*separatortm);
    pack_start(*stren);
    pack_start(*equiltm);
    pack_start(*gamma);
    pack_start(*satur);
    pack_start(*estop);
    pack_start(*scaltm);
    pack_start(*rewei);
    // pack_start(*softradiustm); //unused here, but used for normalize_mean_dt 
//    pack_start(*sensitm);
    ToolParamBlock* const tmBox3 = Gtk::manage(new ToolParamBlock());
    tmBox3->pack_start(*maskusablet, Gtk::PACK_SHRINK, 0);
    tmBox3->pack_start(*maskunusablet, Gtk::PACK_SHRINK, 0);
    tmBox3->pack_start(*recothrest);
    tmBox3->pack_start(*lowthrest);
    tmBox3->pack_start(*higthrest);
    tmBox3->pack_start(*decayt);
   // colBox3->pack_start(*invmaskc);
    exprecovt->add(*tmBox3, false);
    pack_start(*exprecovt, false, false);
    
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

void LocallabTone::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask)
{
    tmMask = showmasktmMethod->get_active_row_number();
}

void LocallabTone::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_TONEMAP_TOOLTIP"));
        recothrest->set_tooltip_text(M("TP_LOCALLAB_RECOTHRES02_TOOLTIP"));
        exprecovt->set_tooltip_markup(M("TP_LOCALLAB_MASKRESTM_TOOLTIP"));
        equiltm->set_tooltip_text(M("TP_LOCALLAB_EQUILTM_TOOLTIP"));
        repartm->set_tooltip_text(M("TP_LOCALLAB_REPARTM_TOOLTIP"));
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
        radmasktm->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        mask2tmCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        Lmasktmshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        masktmCurveEditorG->set_tooltip_markup(M("TP_LOCALLAB_MASKCURVE_TOOLTIP"));
        gammasktm->set_tooltip_text(M("TP_LOCALLAB_GAMMASK_TOOLTIP"));
        chromasktm->set_tooltip_text(M("TP_LOCALLAB_CHROMASK_TOOLTIP"));
        slomasktm->set_tooltip_text(M("TP_LOCALLAB_SLOMASK_TOOLTIP"));
        lapmasktm->set_tooltip_text(M("TP_LOCALLAB_LAPRAD1_TOOLTIP"));
        decayt->set_tooltip_text(M("TP_LOCALLAB_MASKDECAY_TOOLTIP"));
        lowthrest->set_tooltip_text(M("TP_LOCALLAB_MASKLOWTHRESTM_TOOLTIP"));
        higthrest->set_tooltip_text(M("TP_LOCALLAB_MASKHIGTHRESTM_TOOLTIP"));

    } else {
        exp->set_tooltip_text("");
        equiltm->set_tooltip_text("");
        recothrest->set_tooltip_text("");
        repartm->set_tooltip_text("");
        gamma->set_tooltip_text("");
        estop->set_tooltip_text("");
        scaltm->set_tooltip_text("");
        rewei->set_tooltip_text("");
        sensitm->set_tooltip_text("");
        expmasktm->set_tooltip_markup("");
        CCmasktmshape->setTooltip("");
        LLmasktmshape->setTooltip("");
        HHmasktmshape->setTooltip("");
        blendmasktm->set_tooltip_text("");
        radmasktm->set_tooltip_text("");
        mask2tmCurveEditorG->set_tooltip_text("");
        Lmasktmshape->setTooltip("");
        mask2tmCurveEditorG->set_tooltip_text("");
        Lmasktmshape->setTooltip("");
        masktmCurveEditorG->set_tooltip_markup("");
        gammasktm->set_tooltip_text("");
        chromasktm->set_tooltip_text("");
        slomasktm->set_tooltip_text("");
        lapmasktm->set_tooltip_text("");
        exprecovt->set_tooltip_markup("");
        decayt->set_tooltip_text("");
        lowthrest->set_tooltip_text("");
        higthrest->set_tooltip_text("");
    }
}

void LocallabTone::setDefaultExpanderVisibility()
{
    exprecovt->set_expanded(false);
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

        exp->set_visible(spot.visitonemap);
        exp->setEnabled(spot.exptonemap);
        complexity->set_active(spot.complextonemap);

        amount->setValue(spot.amount);
        stren->setValue(spot.stren);
        repartm->setValue(spot.repartm);
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
        recothrest->setValue((double)spot.recothrest);
        lowthrest->setValue((double)spot.lowthrest);
        higthrest->setValue((double)spot.higthrest);
        decayt->setValue((double)spot.decayt);
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
        spot.repartm = repartm->getValue();
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
        spot.recothrest = recothrest->getValue();
        spot.lowthrest = lowthrest->getValue();
        spot.higthrest = higthrest->getValue();
        spot.decayt = decayt->getValue();
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
        repartm->setDefault(defSpot.repartm);
        rewei->setDefault((double)defSpot.rewei);
        softradiustm->setDefault(defSpot.softradiustm);
        sensitm->setDefault((double)defSpot.sensitm);
        blendmasktm->setDefault((double)defSpot.blendmasktm);
        lapmasktm->setDefault(defSpot.lapmasktm);
        radmasktm->setDefault(defSpot.radmasktm);
        chromasktm->setDefault(defSpot.chromasktm);
        gammasktm->setDefault(defSpot.gammasktm);
        slomasktm->setDefault(defSpot.slomasktm);
        recothrest->setDefault((double)defSpot.recothrest);
        lowthrest->setDefault((double)defSpot.lowthrest);
        higthrest->setDefault((double)defSpot.higthrest);
        decayt->setDefault((double)defSpot.decayt);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabTone::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled() && listener) {
        const auto spName = " (" + escapeHtmlChars(getSpotName()) + ")";

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
        } else if (a == repartm) {
            listener->panelChanged(Evlocallabrepartm, repartm->getTextValue() + spName);
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
        } else if (a == recothrest) {
            listener->panelChanged(Evlocallabrecothrest, recothrest->getTextValue() + spName);
        } else if (a == lowthrest) {
            listener->panelChanged(Evlocallablowthrest, lowthrest->getTextValue() + spName);
        } else if (a == higthrest) {
            listener->panelChanged(Evlocallabhigthrest, higthrest->getTextValue() + spName);
        } else if (a == decayt) {
            listener->panelChanged(Evlocallabdecayt, decayt->getTextValue() + spName);
        }
    }
}

void LocallabTone::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled() && listener) {
        const auto spName = M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")";

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
                               + " (" + escapeHtmlChars(getSpotName()) + ")");
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
//    CCmasktmshape->setCurve(defSpot.CCmasktmcurve);
//    LLmasktmshape->setCurve(defSpot.LLmasktmcurve);
//    HHmasktmshape->setCurve(defSpot.HHmasktmcurve);
//    blendmasktm->setValue((double)defSpot.blendmasktm);
//    radmasktm->setValue(defSpot.radmasktm);
//    chromasktm->setValue(defSpot.chromasktm);
//    Lmasktmshape->setCurve(defSpot.Lmasktmcurve);
    recothrest->setValue(defSpot.recothrest);
    lowthrest->setValue(defSpot.lowthrest);
    higthrest->setValue(defSpot.higthrest);
    decayt->setValue(defSpot.decayt);

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
            exprecovt->hide();
            decayt->hide();
            maskusablet->hide();
            maskunusablet->hide();

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
            exprecovt->show();
            decayt->hide();
            if (enatmMask->get_active()) {
                maskusablet->show();
                maskunusablet->hide();
                
            } else {
                maskusablet->hide();
                maskunusablet->show();
            }

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
            exprecovt->show();
            decayt->show();
            if (enatmMask->get_active()) {
                maskusablet->show();
                maskunusablet->hide();
                
            } else {
                maskusablet->hide();
                maskunusablet->show();
            }
    }
}

void LocallabTone::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz)
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabequiltm,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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

    if(exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
    }
}

void LocallabTone::enatmMaskChanged()
{
    if (enatmMask->get_active()) {
        maskusablet->show();
        maskunusablet->hide();
    } else {
        maskusablet->hide();
        maskunusablet->show();
    }
    
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enatmMask->get_active()) {
                listener->panelChanged(EvLocallabEnatmMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnatmMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnatmMaskaft,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
    dehazeSaturation(Gtk::manage(new Adjuster(M("TP_DEHAZE_SATURATION"), 0, 100, 1, 50))),
    retiFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RETIFRA")))),
    str(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STR"), 0., 100., 0.2, 0.))),
    loglin(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_LOGLIN")))),
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
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
    exprecovr(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_DENOI2_EXP")))),
    maskusabler(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUSABLE")))),
    maskunusabler(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUNUSABLE")))),
    recothresr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKRECOTHRES"), 1., 2., 0.01, 1.))),
    lowthresr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHRLOW2"), 1., 80., 0.5, 12.))),
    higthresr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHR2"), 20., 99., 0.5, 85.))),
    decayr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKDDECAY"), 0.5, 4., 0.1, 2.))),
    expmaskreti(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWR")))),
    showmaskretiMethod(Gtk::manage(new MyComboBoxText())),
    enaretiMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    enaretiMasktmap(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TM_MASK")))),
 //   maskretiCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    maskretiCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    CCmaskretishape(static_cast<FlatCurveEditor*>(maskretiCurveEditorG->addCurve(CT_Flat, "C", nullptr, false, false))),
    LLmaskretishape(static_cast<FlatCurveEditor*>(maskretiCurveEditorG->addCurve(CT_Flat, "L", nullptr, false, false))),
    HHmaskretishape(static_cast<FlatCurveEditor *>(maskretiCurveEditorG->addCurve(CT_Flat, "LC(h)", nullptr, false, true))),
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
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
    const LocallabParams::LocallabSpot defSpot;

    // Parameter Retinex specific widgets
    dehaz->setAdjusterListener(this);

    dehazeSaturation->setAdjusterListener(this);
    depth->setAdjusterListener(this);

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

    recothresr->setAdjusterListener(this);
    lowthresr->setAdjusterListener(this);
    higthresr->setAdjusterListener(this);
    decayr->setAdjusterListener(this);
    setExpandAlignProperties(exprecovr, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

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
    pack_start(*sensih);
    ToolParamBlock* const auxBox = Gtk::manage(new ToolParamBlock());
//    Gtk::Frame* const dehaFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_DEHAFRA")));
    dehaFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const dehaBox = Gtk::manage(new ToolParamBlock());
    dehaBox->pack_start(*dehaz);
    dehaBox->pack_start(*depth);
    dehaBox->pack_start(*dehazeSaturation);
    dehaFrame->add(*dehaBox);
    auxBox->add(*dehaFrame);
    ToolParamBlock* const deharetiBox = Gtk::manage(new ToolParamBlock());
    deharetiBox->pack_start(*str);
    deharetiBox->pack_start(*loglin);
    retiFrame->add(*deharetiBox);
    auxBox->add(*retiFrame);
 //   ToolParamBlock* const scopeBox = Gtk::manage(new ToolParamBlock());
 //   scopeBox->pack_start(*sensih);
 //   auxBox->add(*scopeBox);
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
  //  toolretiBox->pack_start(*chrrt);
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
    ToolParamBlock* const reBox3 = Gtk::manage(new ToolParamBlock());
    reBox3->pack_start(*maskusabler, Gtk::PACK_SHRINK, 0);
    reBox3->pack_start(*maskunusabler, Gtk::PACK_SHRINK, 0);
    reBox3->pack_start(*recothresr);
    reBox3->pack_start(*lowthresr);
    reBox3->pack_start(*higthresr);
    reBox3->pack_start(*decayr);
   // colBox3->pack_start(*invmaskc);
    exprecovr->add(*reBox3, false);

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
    retiBox->pack_start(*exprecovr, false, false);
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

void LocallabRetinex::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask)
{
    retiMask = showmaskretiMethod->get_active_row_number();
}

void LocallabRetinex::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        dehaFrame->set_tooltip_text(M("TP_LOCALLAB_DEHAZFRAME_TOOLTIP"));
        dehaz->set_tooltip_text(M("TP_LOCALLAB_DEHAZ_TOOLTIP"));
        retiFrame->set_tooltip_text(M("TP_LOCALLAB_RETIFRAME_TOOLTIP"));
        exprecovr->set_tooltip_markup(M("TP_LOCALLAB_MASKRESRETI_TOOLTIP"));
        loglin->set_tooltip_text(M("TP_LOCALLAB_RETI_LOGLIN_TOOLTIP"));
        sensih->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        fftwreti->set_tooltip_text(M("TP_LOCALLAB_LC_FFTW_TOOLTIP"));
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
        transLabels2->set_tooltip_markup(M("TP_LOCALLAB_TLABEL_TOOLTIP"));
        cTgainshape->setTooltip(M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));
        expmaskreti->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        enaretiMasktmap->set_tooltip_markup(M("TP_LOCALLAB_ENARETIMASKTMAP_TOOLTIP"));
        CCmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskretishape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        blendmaskreti->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        radmaskreti->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        mask2retiCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        Lmaskretishape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        maskretiCurveEditorG->set_tooltip_markup(M("TP_LOCALLAB_MASKCURVE_TOOLTIP"));
        gammaskreti->set_tooltip_text(M("TP_LOCALLAB_GAMMASK_TOOLTIP"));
        chromaskreti->set_tooltip_text(M("TP_LOCALLAB_CHROMASK_TOOLTIP"));
        slomaskreti->set_tooltip_text(M("TP_LOCALLAB_SLOMASK_TOOLTIP"));
        lapmaskreti->set_tooltip_text(M("TP_LOCALLAB_LAPRAD1_TOOLTIP"));
        decayr->set_tooltip_text(M("TP_LOCALLAB_MASKDECAY_TOOLTIP"));
        lowthresr->set_tooltip_text(M("TP_LOCALLAB_MASKLOWTHRESRETI_TOOLTIP"));
        higthresr->set_tooltip_text(M("TP_LOCALLAB_MASKHIGTHRESRETI_TOOLTIP"));

    } else {
        dehaFrame->set_tooltip_text("");
        dehaz->set_tooltip_text("");
        retiFrame->set_tooltip_text("");
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
        mMLabels->set_tooltip_markup("");
        transLabels->set_tooltip_markup("");
        transLabels2->set_tooltip_markup("");
        cTgainshape->setTooltip("");
        expmaskreti->set_tooltip_markup("");
        enaretiMasktmap->set_tooltip_markup("");
        CCmaskretishape->setTooltip("");
        LLmaskretishape->setTooltip("");
        HHmaskretishape->setTooltip("");
        blendmaskreti->set_tooltip_text("");
        radmaskreti->set_tooltip_text("");
        mask2retiCurveEditorG->set_tooltip_text("");
        Lmaskretishape->setTooltip("");
        maskretiCurveEditorG->set_tooltip_markup("");
        gammaskreti->set_tooltip_text("");
        chromaskreti->set_tooltip_text("");
        slomaskreti->set_tooltip_text("");
        lapmaskreti->set_tooltip_text("");
        exprecovr->set_tooltip_markup("");
        decayr->set_tooltip_text("");
        lowthresr->set_tooltip_text("");
        higthresr->set_tooltip_text("");
    }
}

void LocallabRetinex::setDefaultExpanderVisibility()
{
    exprecovr->set_expanded(false);
    expretitools->set_expanded(false);
    expmaskreti->set_expanded(false);
}

void LocallabRetinex::disableListener()
{
    LocallabTool::disableListener();

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

        exp->set_visible(spot.visireti);
        exp->setEnabled(spot.expreti);
        complexity->set_active(spot.complexreti);

        dehaz->setValue((double)spot.dehaz);
        depth->setValue((double)spot.depth);
        dehazeSaturation->setValue((double)spot.dehazeSaturation);
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
        chrrt->setValue(0.);
       // chrrt->setValue(spot.chrrt);
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
        recothresr->setValue((double)spot.recothresr);
        lowthresr->setValue((double)spot.lowthresr);
        higthresr->setValue((double)spot.higthresr);
        decayr->setValue((double)spot.decayr);
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
        spot.dehazeSaturation = dehazeSaturation->getIntValue();
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
        spot.recothresr = recothresr->getValue();
        spot.lowthresr = lowthresr->getValue();
        spot.higthresr = higthresr->getValue();
        spot.decayr = decayr->getValue();
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
        dehazeSaturation->setDefault((double)defSpot.dehazeSaturation);
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
        recothresr->setDefault((double)defSpot.recothresr);
        lowthresr->setDefault((double)defSpot.lowthresr);
        higthresr->setDefault((double)defSpot.higthresr);
        decayr->setDefault((double)defSpot.decayr);
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
                                       dehaz->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == dehazeSaturation) {
            if (listener) {
                listener->panelChanged(EvlocallabdehazeSaturation,
                                       dehazeSaturation->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == depth) {
            if (listener) {
                listener->panelChanged(Evlocallabdepth,
                                       depth->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == str) {
            if (listener) {
                listener->panelChanged(Evlocallabstr,
                                       str->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sensih) {
            if (listener) {
                listener->panelChanged(Evlocallabsensih,
                                       sensih->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == neigh) {
            if (listener) {
                listener->panelChanged(Evlocallabneigh,
                                       neigh->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == vart) {
            if (listener) {
                listener->panelChanged(Evlocallabvart,
                                       vart->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == scalereti) {
            if (listener) {
                listener->panelChanged(Evlocallabscalereti,
                                       scalereti->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == limd) {
            if (listener) {
                listener->panelChanged(Evlocallablimd,
                                       limd->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == offs) {
            if (listener) {
                listener->panelChanged(Evlocallaboffs,
                                       offs->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chrrt) {
            if (listener) {
                listener->panelChanged(Evlocallabchrrt,
                                       chrrt->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == darkness) {
            if (listener) {
                listener->panelChanged(Evlocallabdarkness,
                                       darkness->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lightnessreti) {
            if (listener) {
                listener->panelChanged(Evlocallablightnessreti,
                                       lightnessreti->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == cliptm) {
            if (listener) {
                listener->panelChanged(Evlocallabcliptm,
                                       cliptm->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == softradiusret) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusret,
                                       softradiusret->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == recothresr) {
            
            if (listener) {
                listener->panelChanged(Evlocallabrecothresr,
                                       recothresr->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lowthresr) {
            if (listener) {
                listener->panelChanged(Evlocallablowthresr,
                                       lowthresr->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == higthresr) {
            if (listener) {
                listener->panelChanged(Evlocallabhigthresr,
                                       higthresr->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == decayr) {
            if (listener) {
                listener->panelChanged(Evlocallabdecayr,
                                       decayr->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == blendmaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskreti,
                                       blendmaskreti->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == radmaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskreti,
                                       radmaskreti->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lapmaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskreti,
                                       lapmaskreti->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskreti,
                                       chromaskreti->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == gammaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskreti,
                                       gammaskreti->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == slomaskreti) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskreti,
                                       slomaskreti->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == cTgainshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCTgainCurve,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == CCmaskretishape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == LLmaskretishape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == HHmaskretishape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == Lmaskretishape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskretishape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocenareti,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
    recothresr->setValue(defSpot.recothresr);
    lowthresr->setValue(defSpot.lowthresr);
    higthresr->setValue(defSpot.higthresr);
    decayr->setValue(defSpot.decayr);

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

    recothresr->setValue(defSpot.recothresr);
    lowthresr->setValue(defSpot.lowthresr);
    higthresr->setValue(defSpot.higthresr);
    decayr->setValue(defSpot.decayr);
    enableListener();
    
}

void LocallabRetinex::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            retiFrame->hide();
            retitoolFrame->hide();
            exprecovr->hide();
            decayr->hide();
            maskusabler->hide();
            maskunusabler->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            retiFrame->hide();
            retitoolFrame->hide();
            // Specific Simple mode widgets are shown in Normal mode
            exprecovr->hide();
            decayr->hide();
            maskusabler->hide();
            maskunusabler->hide();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            retiFrame->show();
            retitoolFrame->show();
            exprecovr->show();
            decayr->show();
            if (enaretiMask->get_active()) {
                maskusabler->show();
                maskunusabler->hide();
                
            } else {
                maskusabler->hide();
                maskunusabler->show();
            }
    }
}

void LocallabRetinex::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz)
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

void LocallabRetinex::loglinChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (loglin->get_active()) {
                listener->panelChanged(Evlocallabloglin,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabloglin,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabRetinex::retinexMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabretinexMethod,
                                   retinexMethod->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}

void LocallabRetinex::fftwretiChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (fftwreti->get_active()) {
                listener->panelChanged(Evlocallabfftwreti,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabfftwreti,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabequilret,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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

    if(exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
    }
}

void LocallabRetinex::enaretiMaskChanged()
{
    if (enaretiMask->get_active()) {
        maskusabler->show();
        maskunusabler->hide();

    } else {
        maskusabler->hide();
        maskunusabler->show();
    }
    
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaretiMask->get_active()) {
                listener->panelChanged(EvLocallabEnaretiMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaretiMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaretiMasktmap,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabinversret,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
    shargam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMC"), 0.5, 3.0, 0.05, 1.))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 100))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 0))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 0.4, 2.5, 0.01, 0.75))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 40))),
    inverssha(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    sharFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHARFRAME")))),
    showmasksharMethod(Gtk::manage(new MyComboBoxText()))
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
    // Parameter Sharpening specific widgets
    sharcontrast->setAdjusterListener(this);

    sharradius->setAdjusterListener(this);

    sharamount->setAdjusterListener(this);

    shardamping->setAdjusterListener(this);

    shariter->setAdjusterListener(this);

    sharblur->setAdjusterListener(this);

    shargam->setAdjusterListener(this);

    sensisha->setAdjusterListener(this);

    inversshaConn = inverssha->signal_toggled().connect(sigc::mem_fun(*this, &LocallabSharp::inversshaChanged));

    showmasksharMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasksharMethod->append(M("TP_LOCALLAB_SHOWMODIF2"));
    showmasksharMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmasksharMethod->set_active(0);
    showmasksharMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmasksharMethodConn = showmasksharMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabSharp::showmasksharMethodChanged));

    // Add Sharpening specific widgets to GUI
    pack_start(*sensisha);
    pack_start(*sharcontrast);
    pack_start(*sharblur);
    pack_start(*shargam);
    pack_start(*sharradius);
    pack_start(*sharamount);
    pack_start(*shardamping);
    pack_start(*shariter);
//    pack_start(*sensisha);
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

void LocallabSharp::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask)
{
    sharMask = showmasksharMethod->get_active_row_number();
}

void LocallabSharp::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPSHARP_TOOLTIP"));
        sensisha->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        shargam->set_tooltip_text(M("TP_LOCALLAB_GAMCOL_TOOLTIP"));

    } else {
        exp->set_tooltip_text("");
        sensisha->set_tooltip_text("");
        shargam->set_tooltip_text("");
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

        exp->set_visible(spot.visisharp);
        exp->setEnabled(spot.expsharp);
        complexity->set_active(spot.complexsharp);

        sharcontrast->setValue((double)spot.sharcontrast);
        sharradius->setValue(spot.sharradius);
        sharamount->setValue((double)spot.sharamount);
        shardamping->setValue((double)spot.shardamping);
        shariter->setValue((double)spot.shariter);
        sharblur->setValue(spot.sharblur);
        shargam->setValue(spot.shargam);
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
        spot.shargam = shargam->getValue();
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
        shargam->setDefault(defSpot.shargam);
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
                                       sharcontrast->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sharradius) {
            if (listener) {
                listener->panelChanged(Evlocallabsharradius,
                                       sharradius->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sharamount) {
            if (listener) {
                listener->panelChanged(Evlocallabsharamount,
                                       sharamount->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == shardamping) {
            if (listener) {
                listener->panelChanged(Evlocallabshardamping,
                                       shardamping->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == shariter) {
            if (listener) {
                listener->panelChanged(Evlocallabshariter,
                                       shariter->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sharblur) {
            if (listener) {
                listener->panelChanged(Evlocallabsharblur,
                                       sharblur->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == shargam) {
            if (listener) {
                listener->panelChanged(Evlocallabshargam,
                                       shargam->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sensisha) {
            if (listener) {
                listener->panelChanged(Evlocallabsensis,
                                       sensisha->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocenasharp,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
    shargam->setValue(defSpot.shargam);

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
            shargam->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            sharcontrast->hide();
            sharblur->hide();
            shargam->hide();
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
            shargam->show();
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabinverssha,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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

    if(exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
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
    expresidpyr(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::Box())))),
    residcont(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCONT"), -100, 100, 1, 0))),
    residchro(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDCHRO"), -100., 100., 1., 0.))),
    residsha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDSHA"), -100., 100., 1., 0.))),
    residshathr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDSHATHR"), 0., 100., 1., 30.))),
    residhi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDHI"), -100., 100., 1., 0.))),
    residhithr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RESIDHITHR"), 0., 100., 1., 70.))),
    gamlc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMW"), 0.5, 3., 0.01, 1.))),
    residgam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMSH"), 0.25, 15.0, 0.01, 2.4))),
    residslop(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOSH"), 0.0, 500.0, 0.01, 12.92))),
    sensilc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
    reparw(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGREPART"), 1.0, 100.0, 1., 100.0))),
    clariFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CLARIFRA")))),
    clarilres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARILRES"), -20., 100., 0.5, 0.))),
    claricres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARICRES"), -20., 100., 0.5, 0.))),
    clarisoft(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.5, 1.))),
    origlc(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ORIGLC")))),
    expcontrastpyr(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::Box())))),
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
    expcontrastpyr2(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::Box())))),
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
    //fatres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATRES"), 0., 100., 1., 0.))),
    fftwlc(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FFTW")))),
    exprecovw(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_DENOI2_EXP")))),
    maskusablew(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUSABLE")))),
    maskunusablew(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUNUSABLE")))),
    recothresw(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKRECOTHRES"), 0., 2., 0.01, 1.))),
    lowthresw(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHRLOW"), 1., 80., 0.5, 12.))),
    higthresw(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHR"), 20., 99., 0.5, 85.))),
    decayw(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKDDECAY"), 0.5, 4., 0.1, 2.))),

    expmasklc(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWLC")))),
    showmasklcMethod(Gtk::manage(new MyComboBoxText())),
    enalcMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
//    masklcCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    masklcCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    CCmasklcshape(static_cast<FlatCurveEditor*>(masklcCurveEditorG->addCurve(CT_Flat, "C", nullptr, false, false))),
    LLmasklcshape(static_cast<FlatCurveEditor*>(masklcCurveEditorG->addCurve(CT_Flat, "L", nullptr, false, false))),
    HHmasklcshape(static_cast<FlatCurveEditor *>(masklcCurveEditorG->addCurve(CT_Flat, "LC(h)", nullptr, false, true))),
    blendmasklc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmasklc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromasklc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    mask2lcCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmasklcshape(static_cast<DiagonalCurveEditor*>(mask2lcCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
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
//    wavshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    LocalcurveEditorwav->curveListComplete();

    csThreshold->setAdjusterListener(this);

    levelwav->setAdjusterListener(this);

    Gtk::Box* const LresTitleHBox = Gtk::manage(new Gtk::Box());
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

    gamlc->setAdjusterListener(this);


    residgam->setAdjusterListener(this);

    residslop->setAdjusterListener(this);
    residslop->setLogScale(16, 0);

    sensilc->setAdjusterListener(this);

    reparw->setAdjusterListener(this);

    clariFrame->set_label_align(0.025, 0.5);

    clarilres->setAdjusterListener(this);

    claricres->setAdjusterListener(this);

    clarisoft->setLogScale(10, 0);
    clarisoft->setAdjusterListener(this);

    origlcConn = origlc->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::origlcChanged));

    Gtk::Box *TittleVBox;
    TittleVBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    TittleVBox->set_spacing(2);

    Gtk::Box* const LCTitleHBox = Gtk::manage(new Gtk::Box());
    Gtk::Box* const LCTitleHBox11 = Gtk::manage(new Gtk::Box());
    Gtk::Label* const LCLabel = Gtk::manage(new Gtk::Label());
    Gtk::Label* const LCLabel11 = Gtk::manage(new Gtk::Label());

    LCLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYR")) + Glib::ustring("</b>"));
    LCLabel11->set_markup(escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYRLAB")));
    LCLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    LCLabel11->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    LCTitleHBox->pack_start(*LCLabel, Gtk::PACK_EXPAND_WIDGET, 0);
    LCTitleHBox11->pack_start(*LCLabel11, Gtk::PACK_EXPAND_WIDGET, 0);
    TittleVBox->pack_start(*LCTitleHBox, Gtk::PACK_SHRINK);
    TittleVBox->pack_start(*LCTitleHBox11, Gtk::PACK_SHRINK);
    expcontrastpyr->setLabel(TittleVBox);
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
//    wavshapeedg->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

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

    Gtk::Box *TittleVBox2;
    TittleVBox2 = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    TittleVBox2->set_spacing(2);

    Gtk::Box* const LCTitleHBox2 = Gtk::manage(new Gtk::Box());
    Gtk::Box* const LCTitleHBox22 = Gtk::manage(new Gtk::Box());
    Gtk::Label* const LCLabel2 = Gtk::manage(new Gtk::Label());
    Gtk::Label* const LCLabel22 = Gtk::manage(new Gtk::Label());

    LCLabel2->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYR2")) + Glib::ustring("</b>"));
    LCLabel22->set_markup(escapeHtmlChars(M("TP_LOCALLAB_LOC_CONTRASTPYR2LAB")));
    LCLabel2->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    LCLabel22->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    LCTitleHBox2->pack_start(*LCLabel2, Gtk::PACK_EXPAND_WIDGET, 0);
    LCTitleHBox22->pack_start(*LCLabel22, Gtk::PACK_EXPAND_WIDGET, 0);
    TittleVBox2->pack_start(*LCTitleHBox2, Gtk::PACK_SHRINK);
    TittleVBox2->pack_start(*LCTitleHBox22, Gtk::PACK_SHRINK);
    expcontrastpyr2->setLabel(TittleVBox2);
    
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
//    wavshapecomp->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    LocalcurveEditorwavcomp->curveListComplete();

    //fatres->setAdjusterListener(this);

    fftwlcConn = fftwlc->signal_toggled().connect(sigc::mem_fun(*this, &LocallabContrast::fftwlcChanged));

    recothresw->setAdjusterListener(this);
    lowthresw->setAdjusterListener(this);
    higthresw->setAdjusterListener(this);
    decayw->setAdjusterListener(this);
    setExpandAlignProperties(exprecovw, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);


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
//    LLmasklcshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

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
    pack_start(*sensilc);
    pack_start(*reparw);
    pack_start(*localcontMethod);
    pack_start(*lcradius);
    pack_start(*lcamount);
    pack_start(*lcdarkness);
    pack_start(*lclightness);
    pack_start(*csThreshold);
    ToolParamBlock* const coBox = Gtk::manage(new ToolParamBlock());
    coBox->pack_start(*sigmalc);
    coBox->pack_start(*LocalcurveEditorwav, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
   // coBox->pack_start(*csThreshold);
    contFrame->add(*coBox);
    pack_start(*contFrame);
    // pack_start(*levelwav);
    ToolParamBlock* const resiBox = Gtk::manage(new ToolParamBlock());
    resiBox->pack_start(*residcont);
    resiBox->pack_start(*residchro);
    Gtk::Frame* const shresFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHRESFRA")));
    shresFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const shresBox = Gtk::manage(new ToolParamBlock());
    Gtk::Separator* const separatorsh = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    shresBox->pack_start(*residsha);
    shresBox->pack_start(*residshathr);
    shresBox->pack_start(*residhi);
    shresBox->pack_start(*residhithr);
    shresBox->pack_start(*separatorsh);
    shresBox->pack_start(*residgam);
    shresBox->pack_start(*residslop);
    shresFrame->add(*shresBox);
    resiBox->pack_start(*shresFrame);
    expresidpyr->add(*resiBox, false);
    pack_start(*expresidpyr);
//    pack_start(*sensilc);
    Gtk::Separator* const separatorcontr = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    pack_start(*separatorcontr);
    ToolParamBlock* const clariBox = Gtk::manage(new ToolParamBlock());
    clariBox->pack_start(*clarilres);
    clariBox->pack_start(*claricres);
    clariBox->pack_start(*clarisoft);
    clariBox->pack_start(*origlc);
    clariFrame->add(*clariBox);
    pack_start(*clariFrame);
    ToolParamBlock* const blurcontBox = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const gradwavFrame = Gtk::manage(new Gtk::Frame());
    gradwavFrame->set_label_align(0.025, 0.5);
    gradwavFrame->set_label_widget(*wavgradl);
    ToolParamBlock* const gradwavBox = Gtk::manage(new ToolParamBlock());
    gradwavBox->pack_start(*sigmalc2);
    gradwavBox->pack_start(*strwav);
    gradwavBox->pack_start(*angwav);
    gradwavFrame->add(*gradwavBox);
    blurcontBox->pack_start(*gradwavFrame);
    Gtk::Frame* const edgFrame = Gtk::manage(new Gtk::Frame());
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
    Gtk::Box* const edbox = Gtk::manage(new Gtk::Box());
    Gtk::Label* const labmedgr = Gtk::manage(new Gtk::Label(M("TP_WAVELET_MEDGREINF") + ":"));
    edbox->pack_start(*labmedgr, Gtk::PACK_SHRINK, 1);
    edbox->pack_start(*localedgMethod);
    edgsBoxshow->pack_start(*edbox);
    Gtk::Separator* const separatoredg2 = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    edgsBoxshow->pack_start(*separatoredg2);
    edgsBoxshow->pack_start(*tloww);
    edgsBoxshow->pack_start(*thigw);
    Gtk::Separator* const separatoredg = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    edgsBoxshow->pack_start(*separatoredg);
    edgsBoxshow->pack_start(*edgw);
    edgsBoxshow->pack_start(*basew);
    Gtk::Box* const ctboxNP = Gtk::manage(new Gtk::Box());
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
    Gtk::Box* const blurlevcontBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    blurlevcontBox->set_spacing(2);
    blurlevcontBox->pack_start(*levelblur);
    blurlevcontBox->pack_start(*sigmabl);
    blurlevcontBox->pack_start(*chromablu);
    blurlevcontBox->pack_start(*LocalcurveEditorwavlev, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    Gtk::Separator* const separatorblu = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    blurlevcontBox->pack_start(*separatorblu);
    blurlevcontBox->pack_start(*residblur);
    blurlevcontBox->pack_start(*blurlc);
    blurlevelFrame->add(*blurlevcontBox);
    blurcontBox->pack_start(*blurlevelFrame);
    expcontrastpyr->add(*blurcontBox, false);
    pack_start(*gamlc);    
    pack_start(*expcontrastpyr);
    ToolParamBlock* const blurcontBox2 = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const contFrame2 = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CONTFRA")));
    contFrame2->set_label_align(0.025, 0.5);
    contFrame2->set_label_widget(*wavcont);
    Gtk::Box* const contlevBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
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
    Gtk::Box* const compreBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
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
    Gtk::Box* const compBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    compBox->set_spacing(2);
    compBox->pack_start(*sigmadc);
    compBox->pack_start(*deltad);
    compBox->pack_start(*LocalcurveEditorwavcomp, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    // Gtk::Separator* const separatorcomp = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    // compBox->pack_start(*separatorcomp);
    // compBox->pack_start(*fatres);
    compFrame->add(*compBox);
    blurcontBox2->pack_start(*compFrame);
    expcontrastpyr2->add(*blurcontBox2, false);
    pack_start(*expcontrastpyr2);
    pack_start(*fftwlc);
    ToolParamBlock* const wwBox3 = Gtk::manage(new ToolParamBlock());
    wwBox3->pack_start(*maskusablew, Gtk::PACK_SHRINK, 0);
    wwBox3->pack_start(*maskunusablew, Gtk::PACK_SHRINK, 0);
    wwBox3->pack_start(*recothresw);
    wwBox3->pack_start(*lowthresw);
    wwBox3->pack_start(*higthresw);
    wwBox3->pack_start(*decayw);
   // colBox3->pack_start(*invmaskc);
    exprecovw->add(*wwBox3, false);
    pack_start(*exprecovw, false, false);

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

void LocallabContrast::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask)
{
    lcMask = showmasklcMethod->get_active_row_number();
}

void LocallabContrast::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        contFrame->set_tooltip_text(M("TP_LOCALLAB_EXPCONTRAST_TOOLTIP"));
        recothresw->set_tooltip_text(M("TP_LOCALLAB_RECOTHRES02_TOOLTIP"));
        LocalcurveEditorwav->set_tooltip_markup(M("TP_LOCALLAB_WAT_LEVELLOCCONTRAST_TOOLTIP"));
        csThreshold->set_tooltip_markup(M("TP_LOCALLAB_WAT_THRESHOLDWAV_TOOLTIP"));
        levelwav->set_tooltip_markup(M("TP_LOCALLAB_LEVELWAV_TOOLTIP"));
        clariFrame->set_tooltip_markup(M("TP_LOCALLAB_CLARI_TOOLTIP"));
        clarisoft->set_tooltip_markup(M("TP_LOCALLAB_CLARISOFT_TOOLTIP"));
        exprecovw->set_tooltip_markup(M("TP_LOCALLAB_MASKRESWAV_TOOLTIP"));
        gamlc->set_tooltip_text(M("TP_LOCALLAB_GAMC_TOOLTIP"));
        wavshape->setTooltip(M("TP_LOCALLAB_WAT_WAVSHAPE_TOOLTIP"));
        clarilres->set_tooltip_text(M("TP_LOCALLAB_WAT_CLARIL_TOOLTIP"));
        claricres->set_tooltip_text(M("TP_LOCALLAB_WAT_CLARIC_TOOLTIP"));
        sigmalc->set_tooltip_text(M("TP_LOCALLAB_WAT_SIGMALC_TOOLTIP"));
        sigmalc2->set_tooltip_text(M("TP_LOCALLAB_WAT_SIGMALC_TOOLTIP"));
        sigmaed->set_tooltip_text(M("TP_LOCALLAB_WAT_SIGMALC_TOOLTIP"));
        sigmabl->set_tooltip_text(M("TP_LOCALLAB_WAT_SIGMALC_TOOLTIP"));
        sigma->set_tooltip_text(M("TP_LOCALLAB_WAT_SIGMALC_TOOLTIP"));
        sigmadc->set_tooltip_text(M("TP_LOCALLAB_WAT_SIGMALC_TOOLTIP"));
        sigmadr->set_tooltip_text(M("TP_LOCALLAB_WAT_SIGMALC_TOOLTIP"));
        origlc->set_tooltip_text(M("TP_LOCALLAB_WAT_ORIGLC_TOOLTIP"));
        strwav->set_tooltip_text(M("TP_LOCALLAB_WAT_STRWAV_TOOLTIP"));
        angwav->set_tooltip_text(M("TP_LOCALLAB_WAT_STRWAV_TOOLTIP"));
        strengthw->set_tooltip_text(M("TP_LOCALLAB_WAT_STRENGTHW_TOOLTIP"));
        LocalcurveEditorwavedg->set_tooltip_markup(M("TP_LOCALLAB_WAT_LOCCONTRASTEDG_TOOLTIP"));
        wavshapeedg->setTooltip(M("TP_LOCALLAB_WAT_LOCCONTRASTEDG_TOOLTIP"));
        gradw->set_tooltip_text(M("TP_LOCALLAB_WAT_GRADW_TOOLTIP"));
        waveshow->set_tooltip_text(M("TP_LOCALLAB_WAT_WAVESHOW_TOOLTIP"));
        LocalcurveEditorwavlev->set_tooltip_markup(M("TP_LOCALLAB_WAT_WAVBLURCURV_TOOLTIP"));
        wavshapelev->setTooltip(M("TP_LOCALLAB_WAT_WAVBLURCURV_TOOLTIP"));
        levelblur->set_tooltip_text(M("TP_LOCALLAB_WAT_WAVLEVELBLUR_TOOLTIP"));
        residblur->set_tooltip_text(M("TP_LOCALLAB_WAT_RESIDBLUR_TOOLTIP"));
        blurlc->set_tooltip_text(M("TP_LOCALLAB_WAT_BLURLC_TOOLTIP"));
        offset->set_tooltip_text(M("TP_LOCALLAB_WAT_CONTOFFSET_TOOLTIP"));
        chromalev->set_tooltip_text(M("TP_LOCALLAB_WAT_CONTCHROMALEV_TOOLTIP"));
        LocalcurveEditorwavcon->set_tooltip_markup(M("TP_LOCALLAB_WAT_WAVCBDL_TOOLTIP"));
        wavshapecon->setTooltip(M("TP_LOCALLAB_WAT_WAVCBDL_TOOLTIP"));
        LocalcurveEditorwavcompre->set_tooltip_markup(M("TP_LOCALLAB_WAT_WAVTM_TOOLTIP"));
        wavshapecompre->setTooltip(M("TP_LOCALLAB_WAT_WAVTM_TOOLTIP"));
        deltad->set_tooltip_text(M("TP_LOCALLAB_WAT_DELTABAL_TOOLTIP"));
        LocalcurveEditorwavcomp->set_tooltip_markup(M("TP_LOCALLAB_WAT_WAVDELTABAL_TOOLTIP"));
        wavshapecomp->setTooltip(M("TP_LOCALLAB_WAT_WAVDELTABAL_TOOLTIP"));
        threswav->set_tooltip_text(M("TP_LOCALLAB_WAT_BALTHRES_TOOLTIP"));
        residcomp->set_tooltip_text(M("TP_LOCALLAB_WAT_RESIDCOMP_TOOLTIP"));

        
        expresidpyr->set_tooltip_text(M("TP_LOCALLAB_WAT_EXPRESID_TOOLTIP"));
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
        masklcCurveEditorG->set_tooltip_markup(M("TP_LOCALLAB_MASKCURVE_TOOLTIP"));
        chromasklc->set_tooltip_text(M("TP_LOCALLAB_CHROMASK_TOOLTIP"));
        sensilc->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        reparw->set_tooltip_text(M("TP_LOCALLAB_REPARW_TOOLTIP"));
        decayw->set_tooltip_text(M("TP_LOCALLAB_MASKDECAY_TOOLTIP"));
        lowthresw->set_tooltip_text(M("TP_LOCALLAB_MASKLOWTHRESWAV_TOOLTIP"));
        higthresw->set_tooltip_text(M("TP_LOCALLAB_MASKHIGTHRESWAV_TOOLTIP"));
    } else {
        contFrame->set_tooltip_text("");
        recothresw->set_tooltip_text("");
        LocalcurveEditorwav->set_tooltip_markup("");
        csThreshold->set_tooltip_markup("");
        expresidpyr->set_tooltip_text("");
        levelwav->set_tooltip_markup("");
        clariFrame->set_tooltip_markup("");
        clarisoft->set_tooltip_markup("");
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
        expmasklc->set_tooltip_markup("");
        CCmasklcshape->setTooltip("");
        LLmasklcshape->setTooltip("");
        HHmasklcshape->setTooltip("");
        blendmasklc->set_tooltip_text("");
        mask2lcCurveEditorG->set_tooltip_text("");
        Lmasklcshape->setTooltip("");
        masklcCurveEditorG->set_tooltip_markup("");
        chromasklc->set_tooltip_text("");
        sensilc->set_tooltip_text("");
        reparw->set_tooltip_text("");
        gamlc->set_tooltip_text("");

        wavshape->setTooltip("");
        clarilres->set_tooltip_text("");
        claricres->set_tooltip_text("");
        sigmalc->set_tooltip_text("");
        sigmalc2->set_tooltip_text("");
        sigmaed->set_tooltip_text("");
        sigmabl->set_tooltip_text("");
        sigma->set_tooltip_text("");
        sigmadc->set_tooltip_text("");
        sigmadr->set_tooltip_text("");
        origlc->set_tooltip_text("");
        strwav->set_tooltip_text("");
        angwav->set_tooltip_text("");
        strengthw->set_tooltip_text("");
        LocalcurveEditorwavedg->set_tooltip_markup("");
        wavshapeedg->setTooltip("");
        gradw->set_tooltip_text("");
        waveshow->set_tooltip_text("");
        LocalcurveEditorwavlev->set_tooltip_markup("");
        wavshapelev->setTooltip("");
        residblur->set_tooltip_text("");
        blurlc->set_tooltip_text("");
        levelblur->set_tooltip_text("");
        offset->set_tooltip_text("");
        chromalev->set_tooltip_text("");
        LocalcurveEditorwavcon->set_tooltip_markup("");
        wavshapecon->setTooltip("");
        LocalcurveEditorwavcompre->set_tooltip_markup("");
        wavshapecompre->setTooltip("");
        deltad->set_tooltip_text("");
        LocalcurveEditorwavcomp->set_tooltip_markup("");
        wavshapecomp->setTooltip("");
        threswav->set_tooltip_text("");
        residcomp->set_tooltip_text("");
        exprecovw->set_tooltip_markup("");
        decayw->set_tooltip_text("");
        lowthresw->set_tooltip_text("");
        higthresw->set_tooltip_text("");

    }
}

void LocallabContrast::setDefaultExpanderVisibility()
{
    expresidpyr->set_expanded(false);
    expcontrastpyr->set_expanded(false);
    expcontrastpyr2->set_expanded(false);
    expmasklc->set_expanded(false);
    exprecovw->set_expanded(false);
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
        gamlc->setValue((double)spot.gamlc);
        residgam->setValue(spot.residgam);
        residslop->setValue(spot.residslop);
        sensilc->setValue((double)spot.sensilc);
        reparw->setValue(spot.reparw);
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
        //fatres->setValue(spot.fatres);
        enalcMask->set_active(spot.enalcMask);
        CCmasklcshape->setCurve(spot.CCmasklccurve);
        LLmasklcshape->setCurve(spot.LLmasklccurve);
        HHmasklcshape->setCurve(spot.HHmasklccurve);
        blendmasklc->setValue((double)spot.blendmasklc);
        radmasklc->setValue(spot.radmasklc);
        chromasklc->setValue(spot.chromasklc);
        Lmasklcshape->setCurve(spot.Lmasklccurve);
        recothresw->setValue((double)spot.recothresw);
        lowthresw->setValue((double)spot.lowthresw);
        higthresw->setValue((double)spot.higthresw);
        decayw->setValue((double)spot.decayw);
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
        spot.gamlc = gamlc->getValue();
        spot.residgam = residgam->getValue();
        spot.residslop = residslop->getValue();
        spot.sensilc = sensilc->getIntValue();
        spot.reparw = reparw->getValue();
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
        //spot.fatres = fatres->getValue();
        spot.fftwlc = fftwlc->get_active();
        spot.enalcMask = enalcMask->get_active();
        spot.CCmasklccurve = CCmasklcshape->getCurve();
        spot.LLmasklccurve = LLmasklcshape->getCurve();
        spot.HHmasklccurve = HHmasklcshape->getCurve();
        spot.blendmasklc = blendmasklc->getIntValue();
        spot.radmasklc = radmasklc->getValue();
        spot.chromasklc = chromasklc->getValue();
        spot.Lmasklccurve = Lmasklcshape->getCurve();
        spot.recothresw = recothresw->getValue();
        spot.lowthresw = lowthresw->getValue();
        spot.higthresw = higthresw->getValue();
        spot.decayw = decayw->getValue();
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
        gamlc->setDefault((double)defSpot.gamlc);
        residgam->setDefault(defSpot.residgam);
        residslop->setDefault(defSpot.residslop);
        sensilc->setDefault((double)defSpot.sensilc);
        reparw->setDefault(defSpot.reparw);
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
        //fatres->setDefault(defSpot.fatres);
        blendmasklc->setDefault((double)defSpot.blendmasklc);
        radmasklc->setDefault(defSpot.radmasklc);
        chromasklc->setDefault(defSpot.chromasklc);
        recothresw->setDefault((double)defSpot.recothresw);
        lowthresw->setDefault((double)defSpot.lowthresw);
        higthresw->setDefault((double)defSpot.higthresw);
        decayw->setDefault((double)defSpot.decayw);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabContrast::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == lcradius) {
            if (listener) {
                listener->panelChanged(Evlocallablcradius,
                                       lcradius->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lcamount) {
            if (listener) {
                listener->panelChanged(Evlocallablcamount,
                                       lcamount->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lcdarkness) {
            if (listener) {
                listener->panelChanged(Evlocallablcdarkness,
                                       lcdarkness->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lclightness) {
            if (listener) {
                listener->panelChanged(Evlocallablclightness,
                                       lclightness->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sigmalc) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmalc,
                                       sigmalc->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == levelwav) {
            if (listener) {
                listener->panelChanged(Evlocallablevelwav,
                                       levelwav->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residcont) {
            if (listener) {
                listener->panelChanged(Evlocallabresidcont,
                                       residcont->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residchro) {
            if (listener) {
                listener->panelChanged(Evlocallabresidchro,
                                       residchro->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residsha) {
            if (listener) {
                listener->panelChanged(Evlocallabresidsha,
                                       residsha->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residshathr) {
            if (listener) {
                listener->panelChanged(Evlocallabresidshathr,
                                       residshathr->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residhi) {
            if (listener) {
                listener->panelChanged(Evlocallabresidhi,
                                       residhi->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residhithr) {
            if (listener) {
                listener->panelChanged(Evlocallabresidhithr,
                                       residhithr->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == gamlc) {
            if (listener) {
                listener->panelChanged(Evlocallabgamlc,
                                       gamlc->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residgam) {
            if (listener) {
                listener->panelChanged(Evlocallabresidgam,
                                       residgam->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residslop) {
            if (listener) {
                listener->panelChanged(Evlocallabresidslop,
                                       residslop->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sensilc) {
            if (listener) {
                listener->panelChanged(Evlocallabsensilc,
                                       sensilc->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == reparw) {
            if (listener) {
                listener->panelChanged(Evlocallabreparw,
                                       reparw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == clarilres) {
            if (listener) {
                listener->panelChanged(Evlocallabclarilres,
                                       clarilres->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == claricres) {
            if (listener) {
                listener->panelChanged(Evlocallabclaricres,
                                       claricres->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == clarisoft) {
            if (listener) {
                listener->panelChanged(Evlocallabclarisoft,
                                       clarisoft->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sigmalc2) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmalc2,
                                       sigmalc2->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == strwav) {
            if (listener) {
                listener->panelChanged(Evlocallabstrwav,
                                       strwav->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == angwav) {
            if (listener) {
                listener->panelChanged(Evlocallabangwav,
                                       angwav->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == strengthw) {
            if (listener) {
                listener->panelChanged(Evlocallabstrengthw,
                                       strengthw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sigmaed) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmaed,
                                       sigmaed->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == gradw) {
            if (listener) {
                listener->panelChanged(Evlocallabgradw,
                                       gradw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == radiusw) {
            if (listener) {
                listener->panelChanged(Evlocallabradiusw,
                                       radiusw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == detailw) {
            if (listener) {
                listener->panelChanged(Evlocallabdetailw,
                                       detailw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == tloww) {
            if (listener) {
                listener->panelChanged(Evlocallabtloww,
                                       tloww->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == thigw) {
            if (listener) {
                listener->panelChanged(Evlocallabthigw,
                                       thigw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == edgw) {
            if (listener) {
                listener->panelChanged(Evlocallabedgw,
                                       edgw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == basew) {
            if (listener) {
                listener->panelChanged(Evlocallabbasew,
                                       basew->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == levelblur) {
            if (listener) {
                listener->panelChanged(Evlocallablevelblur,
                                       levelblur->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sigmabl) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmabl,
                                       sigmabl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromablu) {
            if (listener) {
                listener->panelChanged(Evlocallabchromablu,
                                       chromablu->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residblur) {
            if (listener) {
                listener->panelChanged(Evlocallabresidblur,
                                       residblur->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sigma) {
            if (listener) {
                listener->panelChanged(Evlocallabsigma,
                                       sigma->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == offset) {
            if (listener) {
                listener->panelChanged(Evlocallaboffset,
                                       offset->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromalev) {
            if (listener) {
                listener->panelChanged(Evlocallabchromalev,
                                       chromalev->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sigmadr) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmadr,
                                       sigmadr->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }


        if (a == threswav) {
            if (listener) {
                listener->panelChanged(Evlocallabthreswav,
                                       threswav->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == residcomp) {
            if (listener) {
                listener->panelChanged(Evlocallabresidcomp,
                                       residcomp->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sigmadc) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmadc,
                                       sigmadc->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == deltad) {
            if (listener) {
                listener->panelChanged(Evlocallabdeltad,
                                       deltad->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        //if (a == fatres) {
        //    if (listener) {
        //        listener->panelChanged(Evlocallabfatres,
        //                               fatres->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
        //    }
        //}

        if (a == recothresw) {
            
            if (listener) {
                listener->panelChanged(Evlocallabrecothresw,
                                       recothresw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lowthresw) {
            if (listener) {
                listener->panelChanged(Evlocallablowthresw,
                                       lowthresw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == higthresw) {
            if (listener) {
                listener->panelChanged(Evlocallabhigthresw,
                                       higthresw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == decayw) {
            if (listener) {
                listener->panelChanged(Evlocallabdecayw,
                                       decayw->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == blendmasklc) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmasklc,
                                       blendmasklc->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == radmasklc) {
            if (listener) {
                listener->panelChanged(Evlocallabradmasklc,
                                       radmasklc->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromasklc) {
            if (listener) {
                listener->panelChanged(Evlocallabchromasklc,
                                       chromasklc->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabContrast::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabcsThreshold,
                                   csThreshold->getHistoryString() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}

void LocallabContrast::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == wavshape) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurve,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == wavshapeedg) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurveedg,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == wavshapelev) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvelev,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == wavshapecon) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvecon,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == wavshapecompre) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvecompre,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == wavshapecomp) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvecomp,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == CCmasklcshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmasklcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == LLmasklcshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmasklcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == HHmasklcshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmasklcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == Lmasklcshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmasklcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocenacontrast,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabContrast::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();
    gamlc->setValue(defSpot.gamlc);

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
    //fatres->setValue(defSpot.fatres);
    fftwlc->set_active(defSpot.fftwlc);
    decayw->setValue(defSpot.decayw);

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
    gamlc->setValue(defSpot.gamlc);

    // Set hidden specific GUI widgets in Simple mode to default spot values
    if (defSpot.localcontMethod == "loc") {
        localcontMethod->set_active(0);
    } else if (defSpot.localcontMethod == "wav") {
        localcontMethod->set_active(1);
    }

    showmasklcMethod->set_active(0);
    enalcMask->set_active(defSpot.enalcMask);
//    CCmasklcshape->setCurve(defSpot.CCmasklccurve);
//    LLmasklcshape->setCurve(defSpot.LLmasklccurve);
//    HHmasklcshape->setCurve(defSpot.HHmasklccurve);
//    blendmasklc->setValue((double)defSpot.blendmasklc);
//    radmasklc->setValue(defSpot.radmasklc);
//    chromasklc->setValue(defSpot.chromasklc);
//    Lmasklcshape->setCurve(defSpot.Lmasklccurve);

    // Enable all listeners
    recothresw->setValue(defSpot.recothresw);
    lowthresw->setValue(defSpot.lowthresw);
    higthresw->setValue(defSpot.higthresw);
    decayw->setValue(defSpot.decayw);
    
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
            exprecovw->hide();
            decayw->hide();
            maskusablew->hide();
            maskunusablew->hide();
            gamlc->hide();

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
            exprecovw->show();
            decayw->hide();
            if (enalcMask->get_active()) {
                maskusablew->show();
                maskunusablew->hide();
                
            } else {
                maskusablew->hide();
                maskunusablew->show();
            }
            gamlc->hide();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            localcontMethod->show();
            origlc->show();

            if (localcontMethod->get_active_row_number() != 0) { // Keep widgets hidden when localcontMethod is equal to 0
                expcontrastpyr->show();
                expcontrastpyr2->show();
                gamlc->show();
            }

            if (localcontMethod->get_active_row_number() != 1) { // Keep widget hidden when localcontMethod is equal to 1
                fftwlc->show();
            }

            expmasklc->show();
            exprecovw->show();
            decayw->show();
            if (enalcMask->get_active()) {
                maskusablew->show();
                maskunusablew->hide();
                
            } else {
                maskusablew->hide();
                maskunusablew->show();
            }
            
    }
}

void LocallabContrast::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz)
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
                                   localcontMethod->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}

void LocallabContrast::origlcChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (origlc->get_active()) {
                listener->panelChanged(Evlocallaboriglc,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallaboriglc,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabwavgradl,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabwavedg,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabContrast::localedgMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallablocaledgMethod,
                                   localedgMethod->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabwaveshow,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabContrast::localneiMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallablocalneiMethod,
                                   localneiMethod->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}

void LocallabContrast::wavblurChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (wavblur->get_active()) {
                listener->panelChanged(Evlocallabwavblur,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabwavblur,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabblurlc,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabwavcont,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabwavcompre,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabwavcomp,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabfftwlc,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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

    if(exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
    }
}

void LocallabContrast::enalcMaskChanged()
{
    if (enalcMask->get_active()) {
        maskusablew->show();
        maskunusablew->hide();

    } else {
        maskusablew->hide();
        maskunusablew->show();
    }
    
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enalcMask->get_active()) {
                listener->panelChanged(EvLocallabEnalcMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnalcMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
        csThreshold->hide();
        levelwav->hide();
        expresidpyr->hide();
        clariFrame->hide();
        expcontrastpyr->hide();
        expcontrastpyr2->hide();
        gamlc->hide();
        if (mode == Expert) { // Keep widget hidden in Normal and Simple mode
            fftwlc->show();
        }
    } else if (localcontMethod->get_active_row_number() == 1) {
        lcradius->hide();
        lcamount->hide();
        lcdarkness->hide();
        lclightness->hide();
        contFrame->show();
        csThreshold->show();
        levelwav->show();
        expresidpyr->show();
        clariFrame->show();

        if (mode == Expert) { // Keep widget hidden in Normal and Simple mode
            expcontrastpyr->show();
            expcontrastpyr2->show();
            gamlc->show();
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
    sensicb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
    exprecovcb(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_DENOI2_EXP")))),
    maskusablecb(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUSABLE")))),
    maskunusablecb(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUNUSABLE")))),
    recothrescb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKRECOTHRES"), 1., 2., 0.01, 1.))),
    lowthrescb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHRLOW"), 1., 80., 0.5, 12.))),
    higthrescb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHR"), 20., 99., 0.5, 85.))),
    decaycb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKDDECAY"), 0.5, 4., 0.1, 2.))),
    expmaskcb(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWCB")))),
    showmaskcbMethod(Gtk::manage(new MyComboBoxText())),
    enacbMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
  //  maskcbCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    maskcbCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    CCmaskcbshape(static_cast<FlatCurveEditor*>(maskcbCurveEditorG->addCurve(CT_Flat, "C", nullptr, false, false))),
    LLmaskcbshape(static_cast<FlatCurveEditor*>(maskcbCurveEditorG->addCurve(CT_Flat, "L", nullptr, false, false))),
    HHmaskcbshape(static_cast<FlatCurveEditor *>(maskcbCurveEditorG->addCurve(CT_Flat, "LC(h)", nullptr, false, true))),
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
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
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

    recothrescb->setAdjusterListener(this);
    lowthrescb->setAdjusterListener(this);
    higthrescb->setAdjusterListener(this);
    decaycb->setAdjusterListener(this);
    setExpandAlignProperties(exprecovcb, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

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
    pack_start(*sensicb);

    // Add CBDL specific widgets to GUI
    ToolParamBlock* const levBox = Gtk::manage(new ToolParamBlock());
    Gtk::Box* buttonBox = Gtk::manage(new Gtk::Box());
    buttonBox->set_spacing(2);
    buttonBox->set_homogeneous(true);
    buttonBox->pack_start(*lumacontrastMinusButton);
    buttonBox->pack_start(*lumaneutralButton);
    buttonBox->pack_start(*lumacontrastPlusButton);
    levBox->pack_start(*buttonBox);

    for (const auto adj : multiplier) {
        levBox->pack_start(*adj);
    }

    Gtk::Separator* const separator = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
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
//    pack_start(*sensicb);
    ToolParamBlock* const cbBox3 = Gtk::manage(new ToolParamBlock());
    cbBox3->pack_start(*maskusablecb, Gtk::PACK_SHRINK, 0);
    cbBox3->pack_start(*maskunusablecb, Gtk::PACK_SHRINK, 0);
    cbBox3->pack_start(*recothrescb);
    cbBox3->pack_start(*lowthrescb);
    cbBox3->pack_start(*higthrescb);
    cbBox3->pack_start(*decaycb);
   // colBox3->pack_start(*invmaskc);
    exprecovcb->add(*cbBox3, false);
    pack_start(*exprecovcb, false, false);

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

void LocallabCBDL::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask)
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

        exprecovcb->set_tooltip_markup(M("TP_LOCALLAB_MASKRESCB_TOOLTIP"));
        chromacbdl->set_tooltip_text(M("TP_LOCALLAB_CHROMACB_TOOLTIP"));
        threshold->set_tooltip_text(M("TP_LOCALLAB_CBDL_THRES_TOOLTIP"));
        clarityml->set_tooltip_text(M("TP_LOCALLAB_CBDLCLARI_TOOLTIP"));
        sensicb->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        expmaskcb->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        blendmaskcb->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        radmaskcb->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        mask2cbCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        Lmaskcbshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        maskcbCurveEditorG->set_tooltip_markup(M("TP_LOCALLAB_MASKCURVE_TOOLTIP"));
        gammaskcb->set_tooltip_text(M("TP_LOCALLAB_GAMMASK_TOOLTIP"));
        chromaskcb->set_tooltip_text(M("TP_LOCALLAB_CHROMASK_TOOLTIP"));
        slomaskcb->set_tooltip_text(M("TP_LOCALLAB_SLOMASK_TOOLTIP"));
        lapmaskcb->set_tooltip_text(M("TP_LOCALLAB_LAPRAD1_TOOLTIP"));
        decaycb->set_tooltip_text(M("TP_LOCALLAB_MASKDECAY_TOOLTIP"));
        lowthrescb->set_tooltip_text(M("TP_LOCALLAB_MASKLOWTHRESCB_TOOLTIP"));
        higthrescb->set_tooltip_text(M("TP_LOCALLAB_MASKHIGTHRESCB_TOOLTIP"));
    } else {
        levFrame->set_tooltip_text("");

        for (const auto adj : multiplier) {
            adj->set_tooltip_text("");
        }

        chromacbdl->set_tooltip_text("");
        threshold->set_tooltip_text("");
        clarityml->set_tooltip_text("");
        sensicb->set_tooltip_text("");
        expmaskcb->set_tooltip_markup("");
        CCmaskcbshape->setTooltip("");
        LLmaskcbshape->setTooltip("");
        HHmaskcbshape->setTooltip("");
        blendmaskcb->set_tooltip_text("");
        radmaskcb->set_tooltip_text("");
        mask2cbCurveEditorG->set_tooltip_text("");
        Lmaskcbshape->setTooltip("");
        maskcbCurveEditorG->set_tooltip_markup("");
        gammaskcb->set_tooltip_text("");
        chromaskcb->set_tooltip_text("");
        slomaskcb->set_tooltip_text("");
        lapmaskcb->set_tooltip_text("");
        exprecovcb->set_tooltip_markup("");
        decaycb->set_tooltip_text("");
        lowthrescb->set_tooltip_text("");
        higthrescb->set_tooltip_text("");
    }
}

void LocallabCBDL::setDefaultExpanderVisibility()
{
    exprecovcb->set_expanded(false);
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
        recothrescb->setValue((double)spot.recothrescb);
        lowthrescb->setValue((double)spot.lowthrescb);
        higthrescb->setValue((double)spot.higthrescb);
        decaycb->setValue((double)spot.decaycb);
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
        spot.recothrescb = recothrescb->getValue();
        spot.lowthrescb = lowthrescb->getValue();
        spot.higthrescb = higthrescb->getValue();
        spot.decaycb = decaycb->getValue();
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
        recothrescb->setDefault((double)defSpot.recothrescb);
        lowthrescb->setDefault((double)defSpot.lowthrescb);
        higthrescb->setDefault((double)defSpot.higthrescb);
        decaycb->setDefault((double)defSpot.decaycb);
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
                                       + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromacbdl) {
            if (listener) {
                listener->panelChanged(Evlocallabchromacbdl,
                                       chromacbdl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == threshold) {
            if (listener) {
                listener->panelChanged(EvlocallabThresho,
                                       threshold->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == clarityml) {
            if (listener) {
                listener->panelChanged(EvLocallabclarityml,
                                       clarityml->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == contresid) {
            if (listener) {
                listener->panelChanged(EvLocallabcontresid,
                                       contresid->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == softradiuscb) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiuscb,
                                       softradiuscb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sensicb) {
            if (listener) {
                listener->panelChanged(Evlocallabsensicb,
                                       sensicb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == recothrescb) {
            
            if (listener) {
                listener->panelChanged(Evlocallabrecothrescb,
                                       recothrescb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lowthrescb) {
            if (listener) {
                listener->panelChanged(Evlocallablowthrescb,
                                       lowthrescb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == higthrescb) {
            if (listener) {
                listener->panelChanged(Evlocallabhigthrescb,
                                       higthrescb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == decaycb) {
            if (listener) {
                listener->panelChanged(Evlocallabdecaycb,
                                       decaycb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == blendmaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskcb,
                                       blendmaskcb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == radmaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskcb,
                                       radmaskcb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lapmaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskcb,
                                       lapmaskcb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskcb,
                                       chromaskcb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == gammaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskcb,
                                       gammaskcb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == slomaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskcb,
                                       slomaskcb->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == LLmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == HHmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == Lmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskcbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocenacbdl,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
    decaycb->setValue(defSpot.decaycb);

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
//    CCmaskcbshape->setCurve(defSpot.CCmaskcbcurve);
//    LLmaskcbshape->setCurve(defSpot.LLmaskcbcurve);
//    HHmaskcbshape->setCurve(defSpot.HHmaskcbcurve);
//    blendmaskcb->setValue((double)defSpot.blendmaskcb);
//    radmaskcb->setValue(defSpot.radmaskcb);
//    chromaskcb->setValue(defSpot.chromaskcb);
//    gammaskcb->setValue(defSpot.gammaskcb);
//    slomaskcb->setValue(defSpot.slomaskcb);
//    Lmaskcbshape->setCurve(defSpot.Lmaskcbcurve);
    recothrescb->setValue(defSpot.recothrescb);
    lowthrescb->setValue(defSpot.lowthrescb);
    higthrescb->setValue(defSpot.higthrescb);
    decaycb->setValue(defSpot.decaycb);

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
            exprecovcb->hide();
            decaycb->hide();
            maskusablecb->hide();
            maskunusablecb->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            lapmaskcb->hide();
            // Specific Simple mode widgets are shown in Normal mode
            softradiuscb->show();
            expmaskcb->show();
            exprecovcb->show();
            decaycb->hide();
            if (enacbMask->get_active()) {
                maskusablecb->show();
                maskunusablecb->hide();
                
            } else {
                maskusablecb->hide();
                maskunusablecb->show();
            }

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            softradiuscb->show();
            expmaskcb->show();
            lapmaskcb->show();
            exprecovcb->show();
            decaycb->show();
            if (enacbMask->get_active()) {
                maskusablecb->show();
                maskunusablecb->hide();
                
            } else {
                maskusablecb->hide();
                maskunusablecb->show();
            }
    }
}

void LocallabCBDL::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz)
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

    if(exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
    }
}

void LocallabCBDL::enacbMaskChanged()
{
    if (enacbMask->get_active()) {
        maskusablecb->show();
        maskunusablecb->hide();
    } else {
        maskusablecb->hide();
        maskunusablecb->show();
    }
    
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enacbMask->get_active()) {
                listener->panelChanged(EvLocallabEnacbMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnacbMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
    repar(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGREPART"), 1.0, 100.0, 1., 100.0))),
    ciecam(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CIEC")))),
    autocompute(Gtk::manage(new Gtk::ToggleButton(M("TP_LOCALLAB_LOGAUTO")))),
    logPFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOGPFRA")))),
    blackEv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLACK_EV"), -16.0, 0.0, 0.1, -5.0))),
    whiteEv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WHITE_EV"), 0., 32.0, 0.1, 10.0))),
    fullimage(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FULLIMAGE")))),
    logFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOGFRA")))),
    Autogray(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_AUTOGRAY")))),
    sourceGray(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOURCE_GRAY"), 1.0, 100.0, 0.1, 10.0))),
    sourceabs(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOURCE_ABS"), 0.01, 16384.0, 0.01, 2000.0))),
    sursour(Gtk::manage (new MyComboBoxText ())),
    surHBox(Gtk::manage(new Gtk::Box())),
    log1Frame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOG1FRA")))),
    log2Frame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOG2FRA")))),
    targetGray(Gtk::manage(new Adjuster(M("TP_LOCALLAB_TARGET_GRAY"), 4.0, 80.0, 0.1, 18.0))),
    detail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DETAIL"), 0., 1., 0.01, 0.6))),
    catad(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CATAD"), -100., 100., 0.5, 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-orange-small.png"))))),
    lightl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGLIGHTL"), -100., 100., 0.5, 0.))),
    lightq(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGLIGHTQ"), -100., 100., 0.5, 0.))),
    contl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCONTL"), -100., 100., 0.5, 0.))),
    contq(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCONQL"), -100., 100., 0.5, 0.))),
    contthres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCONTHRES"), -1., 1., 0.01, 0.))),
    colorfl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCOLORFL"), -100., 100., 0.5, 0.))),
    saturl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SATURV"), -100., 100., 0.5, 0.))),
    chroml(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROML"), -100., 100., 0.5, 0.))),
    expL(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_LOGEXP")))),
    //CurveEditorL(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_LOGCONTQ"))),
    //LshapeL(static_cast<DiagonalCurveEditor*>(CurveEditorL->addCurve(CT_Diagonal, "Q(Q)"))),
    targabs(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOURCE_ABS"), 0.01, 16384.0, 0.01, 16.0))),
    surround(Gtk::manage (new MyComboBoxText ())),
    surrHBox(Gtk::manage(new Gtk::Box())),
    baselog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BASELOG"), 1.3, 3., 0.05, 2.))),//, Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    exprecovl(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_DENOI2_EXP")))),
    maskusablel(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUSABLE")))),
    maskunusablel(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUNUSABLE")))),
    recothresl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKRECOTHRES"), 0., 2., 0.01, 1.))),
    lowthresl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHRLOW"), 1., 80., 0.5, 12.))),
    higthresl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHR"), 20., 99., 0.5, 85.))),
    decayl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKDDECAY"), 0.5, 4., 0.1, 2.))),

    sensilog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
    gradlogFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADLOGFRA")))),
    strlog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -2.0, 2.0, 0.05, 0.))),
    anglog(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    expmaskL(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWC")))),
    showmaskLMethod(Gtk::manage(new MyComboBoxText())),
    enaLMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
 //   maskCurveEditorL(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASKCOL"))),
    maskCurveEditorL(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    CCmaskshapeL(static_cast<FlatCurveEditor*>(maskCurveEditorL->addCurve(CT_Flat, "C", nullptr, false, false))),
    LLmaskshapeL(static_cast<FlatCurveEditor*>(maskCurveEditorL->addCurve(CT_Flat, "L", nullptr, false, false))),
    HHmaskshapeL(static_cast<FlatCurveEditor *>(maskCurveEditorL->addCurve(CT_Flat, "LC(h)", nullptr, false, true))),
    blendmaskL(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskL(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskL(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    mask2CurveEditorL(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    LmaskshapeL(static_cast<DiagonalCurveEditor*>(mask2CurveEditorL->addCurve(CT_Diagonal, "L(L)")))
  
    
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
    // Parameter Log encoding specific widgets
    autoconn = autocompute->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::autocomputeToggled));
    const LocallabParams::LocallabSpot defSpot;
    repar->setAdjusterListener(this);

    blackEv->setLogScale(2, -8);
    blackEv->setAdjusterListener(this);

    whiteEv->setLogScale(16, 0);
    whiteEv->setAdjusterListener(this);
    ciecamconn = ciecam->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::ciecamChanged));

    fullimageConn = fullimage->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::fullimageChanged));

    AutograyConn = Autogray->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::AutograyChanged));

    sourceGray->setAdjusterListener(this);
    sourceGray->setLogScale(10, 18, true);

    sourceabs->setLogScale(500, 0);

    sourceabs->setAdjusterListener(this);

    targetGray->setAdjusterListener(this);
    targetGray->setLogScale(10, 18, true);

    detail->setAdjusterListener(this);

    catad->setAdjusterListener(this);

    setExpandAlignProperties(expL, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    saturl->setAdjusterListener(this);

    chroml->setAdjusterListener(this);

    lightl->setAdjusterListener(this);

    lightq->setAdjusterListener(this);
    contl->setAdjusterListener(this);
    contthres->setAdjusterListener(this);

    contq->setAdjusterListener(this);
    colorfl->setAdjusterListener(this);

    //CurveEditorL->setCurveListener(this);

    //LshapeL->setResetCurve(DiagonalCurveType(defSpot.LcurveL.at(0)), defSpot.LcurveL);
    //LshapeL->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    //LshapeL->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    //CurveEditorL->curveListComplete();


    targabs->setLogScale(500, 0);

    targabs->setAdjusterListener(this);

    baselog->setAdjusterListener(this);
    recothresl->setAdjusterListener(this);
    lowthresl->setAdjusterListener(this);
    higthresl->setAdjusterListener(this);
    decayl->setAdjusterListener(this);
    setExpandAlignProperties(exprecovl, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    sensilog->setAdjusterListener(this);

    strlog->setAdjusterListener(this);

    anglog->setAdjusterListener(this);
    
    surHBox->set_spacing (2);
    surHBox->set_tooltip_markup (M ("TP_LOCALLAB_LOGSURSOUR_TOOLTIP"));
    Gtk::Label* surLabel = Gtk::manage (new Gtk::Label (M ("TP_COLORAPP_SURROUND") + ":"));
    surHBox->pack_start (*surLabel, Gtk::PACK_SHRINK);
    sursour->append (M ("TP_COLORAPP_SURROUND_AVER"));
    sursour->append (M ("TP_COLORAPP_SURROUND_DIM"));
    sursour->append (M ("TP_COLORAPP_SURROUND_DARK"));
    sursour->set_active (0);
    surHBox->pack_start (*sursour);
    sursourconn = sursour->signal_changed().connect ( sigc::mem_fun (*this, &LocallabLog::sursourChanged) );



    surrHBox->set_spacing (2);
    surrHBox->set_tooltip_markup (M ("TP_COLORAPP_SURROUND_TOOLTIP"));
    Gtk::Label* surrLabel = Gtk::manage (new Gtk::Label (M ("TP_COLORAPP_SURROUND") + ":"));
    surrHBox->pack_start (*surrLabel, Gtk::PACK_SHRINK);
    surround->append (M ("TP_COLORAPP_SURROUND_AVER"));
    surround->append (M ("TP_COLORAPP_SURROUND_DIM"));
    surround->append (M ("TP_COLORAPP_SURROUND_DARK"));
    surround->append (M ("TP_COLORAPP_SURROUND_EXDARK"));
    surround->set_active (0);
    surrHBox->pack_start (*surround);
    surroundconn = surround->signal_changed().connect ( sigc::mem_fun (*this, &LocallabLog::surroundChanged) );

    setExpandAlignProperties(expmaskL, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskLMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskLMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskLMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskLMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskLMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskLMethod->set_active(0);
    showmaskLMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskLMethodConn  = showmaskLMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabLog::showmaskLMethodChanged));


    enaLMaskConn = enaLMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabLog::enaLMaskChanged));

    maskCurveEditorL->setCurveListener(this);

    CCmaskshapeL->setIdentityValue(0.);
    CCmaskshapeL->setResetCurve(FlatCurveType(defSpot.CCmaskcurveL.at(0)), defSpot.CCmaskcurveL);
    CCmaskshapeL->setBottomBarColorProvider(this, 1);

    LLmaskshapeL->setIdentityValue(0.);
    LLmaskshapeL->setResetCurve(FlatCurveType(defSpot.LLmaskcurveL.at(0)), defSpot.LLmaskcurveL);
    LLmaskshapeL->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskshapeL->setIdentityValue(0.);
    HHmaskshapeL->setResetCurve(FlatCurveType(defSpot.HHmaskcurveL.at(0)), defSpot.HHmaskcurveL);
    HHmaskshapeL->setCurveColorProvider(this, 2);
    HHmaskshapeL->setBottomBarColorProvider(this, 2);

    maskCurveEditorL->curveListComplete();

    blendmaskL->setAdjusterListener(this);
    radmaskL->setAdjusterListener(this);
    chromaskL->setAdjusterListener(this);

    mask2CurveEditorL->setCurveListener(this);

    LmaskshapeL->setResetCurve(DiagonalCurveType(defSpot.LmaskcurveL.at(0)), defSpot.LmaskcurveL);
    LmaskshapeL->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    LmaskshapeL->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2CurveEditorL->curveListComplete();

    // Add Log encoding specific widgets to GUI
    pack_start(*sensilog);
    pack_start(*repar);
    pack_start(*ciecam);
    logPFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const logPBox = Gtk::manage(new ToolParamBlock());
    logPBox->pack_start(*autocompute);
    logPBox->pack_start(*blackEv);
    logPBox->pack_start(*whiteEv);
    logPBox->pack_start(*fullimage);
    logPFrame->add(*logPBox);
    pack_start(*logPFrame);
//    Gtk::Frame* const logFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOGFRA")));
    logFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const logFBox = Gtk::manage(new ToolParamBlock());
    logFBox->pack_start(*Autogray);
    logFBox->pack_start(*sourceGray);
    logFBox->pack_start(*sourceabs);
    logFBox->pack_start (*surHBox);
//    logFBox->pack_start(*baselog);
    logFrame->add(*logFBox);
    pack_start(*logFrame);
    log1Frame->set_label_align(0.025, 0.5);
    ToolParamBlock* const logP1Box = Gtk::manage(new ToolParamBlock());
    logP1Box->pack_start(*detail);
    logP1Box->pack_start(*contl);
    logP1Box->pack_start(*contthres);
    logP1Box->pack_start(*saturl);
    Gtk::Separator* const separatorchro = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    ToolParamBlock* const logP11Box = Gtk::manage(new ToolParamBlock());
    logP11Box->pack_start(*lightl);
    logP11Box->pack_start(*lightq);
    logP11Box->pack_start(*contq);
    logP11Box->pack_start(*separatorchro);
    logP11Box->pack_start(*chroml);
    logP11Box->pack_start(*colorfl);
    expL->add(*logP11Box, false);
    logP1Box->pack_start(*expL, false, false);
    
//    logP1Box->pack_start(*CurveEditorL, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    log1Frame->add(*logP1Box);
    pack_start(*log1Frame);
    log2Frame->set_label_align(0.025, 0.5);    
    ToolParamBlock* const logP2Box = Gtk::manage(new ToolParamBlock());
    logP2Box->pack_start(*targetGray);
    logP2Box->pack_start(*targabs);
    logP2Box->pack_start(*catad);
    logP2Box->pack_start (*surrHBox);
    ToolParamBlock* const logBox3 = Gtk::manage(new ToolParamBlock());
    logBox3->pack_start(*maskusablel, Gtk::PACK_SHRINK, 0);
    logBox3->pack_start(*maskunusablel, Gtk::PACK_SHRINK, 0);
    logBox3->pack_start(*recothresl);
    logBox3->pack_start(*lowthresl);
    logBox3->pack_start(*higthresl);
    logBox3->pack_start(*decayl);
   // colBox3->pack_start(*invmaskc);
    exprecovl->add(*logBox3, false);
    
    ToolParamBlock* const logP3Box = Gtk::manage(new ToolParamBlock());
    logP3Box->pack_start(*showmaskLMethod, Gtk::PACK_SHRINK, 4);
    logP3Box->pack_start(*enaLMask, Gtk::PACK_SHRINK, 0);
    logP3Box->pack_start(*maskCurveEditorL, Gtk::PACK_SHRINK, 4);
    logP3Box->pack_start(*blendmaskL);
    logP3Box->pack_start(*radmaskL);
    logP3Box->pack_start(*chromaskL);
    logP3Box->pack_start(*mask2CurveEditorL, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    expmaskL->add(*logP3Box, false);


    log2Frame->add(*logP2Box);
    pack_start(*log2Frame);
    pack_start(*exprecovl, false, false);
    
//    pack_start(*baselog);
//    pack_start(*sensilog);
    pack_start(*expmaskL, false, false);
    
 //   Gtk::Frame* const gradlogFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADLOGFRA")));
    gradlogFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const gradlogBox = Gtk::manage(new ToolParamBlock());
    gradlogBox->pack_start(*strlog);
    gradlogBox->pack_start(*anglog);
    gradlogFrame->add(*gradlogBox);
    pack_start(*gradlogFrame);
}

LocallabLog::~LocallabLog()
{
    delete maskCurveEditorL;
    delete mask2CurveEditorL;
    //delete CurveEditorL;
    
}

void LocallabLog::setDefaultExpanderVisibility()
{
    exprecovl->set_expanded(false);
    expmaskL->set_expanded(false);
    expL->set_expanded(false);

}

void LocallabLog::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_LOGENCOD_TOOLTIP"));
        repar->set_tooltip_text(M("TP_LOCALLAB_LOGREPART_TOOLTIP"));
        recothresl->set_tooltip_text(M("TP_LOCALLAB_RECOTHRES02_TOOLTIP"));
        logPFrame->set_tooltip_text(M("TP_LOCALLAB_LOGFRAME_TOOLTIP"));
        logFrame->set_tooltip_text(M("TP_LOCALLAB_LOGSCENE_TOOLTIP"));
        log1Frame->set_tooltip_text(M("TP_LOCALLAB_LOGIMAGE_TOOLTIP"));
        log2Frame->set_tooltip_text(M("TP_LOCALLAB_LOGVIEWING_TOOLTIP"));
        autocompute->set_tooltip_text(M("TP_LOCALLAB_LOGAUTO_TOOLTIP"));
        Autogray->set_tooltip_text(M("TP_LOCALLAB_LOGAUTOGRAY_TOOLTIP"));
    //    blackEv->set_tooltip_text(M("TP_LOCALLAB_LOGBLACKWHEV_TOOLTIP"));
    //    whiteEv->set_tooltip_text(M("TP_LOCALLAB_LOGBLACKWHEV_TOOLTIP"));
        exprecovl->set_tooltip_markup(M("TP_LOCALLAB_MASKRELOG_TOOLTIP"));
        blackEv->set_tooltip_text("");
        whiteEv->set_tooltip_text("");
        sourceGray->set_tooltip_text(M("TP_LOCALLAB_JZLOGYBOUT_TOOLTIP"));
        sourceabs->set_tooltip_text(M("TP_COLORAPP_ADAPSCEN_TOOLTIP"));
        targabs->set_tooltip_text(M("TP_COLORAPP_VIEWING_ABSOLUTELUMINANCE_TOOLTIP"));
        targetGray->set_tooltip_text(M("TP_COLORAPP_YBOUT_TOOLTIP"));
        baselog->set_tooltip_text(M("TP_LOCALLAB_LOGBASE_TOOLTIP"));
        strlog->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
        anglog->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
        contl->set_tooltip_text(M("TP_LOCALLAB_LOGCONTL_TOOLTIP"));
        contq->set_tooltip_text(M("TP_LOCALLAB_LOGCONTQ_TOOLTIP"));
        contthres->set_tooltip_text(M("TP_LOCALLAB_LOGCONTTHRES_TOOLTIP"));
        colorfl->set_tooltip_text(M("TP_LOCALLAB_LOGCOLORF_TOOLTIP"));
        lightl->set_tooltip_text(M("TP_LOCALLAB_LOGLIGHTL_TOOLTIP"));        
        lightq->set_tooltip_text(M("TP_LOCALLAB_LOGLIGHTQ_TOOLTIP"));        
        saturl->set_tooltip_text(M("TP_LOCALLAB_LOGSATURL_TOOLTIP"));
        chroml->set_tooltip_text(M("TP_COLORAPP_CHROMA_TOOLTIP"));
        detail->set_tooltip_text(M("TP_LOCALLAB_LOGDETAIL_TOOLTIP"));
        catad->set_tooltip_text(M("TP_LOCALLAB_LOGCATAD_TOOLTIP"));
        sensilog->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        fullimage->set_tooltip_text(M("TP_LOCALLAB_FULLIMAGELOG_TOOLTIP"));
        ciecam->set_tooltip_text(M("TP_LOCALLAB_CIECAMLOG_TOOLTIP"));
        expmaskL->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmaskshapeL->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskshapeL->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskshapeL->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        blendmaskL->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        radmaskL->set_tooltip_text(M("TP_LOCALLAB_LAPRAD2_TOOLTIP"));
        chromaskL->set_tooltip_text(M("TP_LOCALLAB_CHROMASK_TOOLTIP"));
//        mask2CurveEditorL->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        LmaskshapeL->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        decayl->set_tooltip_text(M("TP_LOCALLAB_MASKDECAY_TOOLTIP"));
        lowthresl->set_tooltip_text(M("TP_LOCALLAB_MASKLOWTHRESL_TOOLTIP"));
        higthresl->set_tooltip_text(M("TP_LOCALLAB_MASKHIGTHRESL_TOOLTIP"));


    } else {
        exp->set_tooltip_text("");
        repar->set_tooltip_text("");
        recothresl->set_tooltip_text("");
        logPFrame->set_tooltip_text("");
        logFrame->set_tooltip_text("");
        log1Frame->set_tooltip_text("");
        log2Frame->set_tooltip_text("");
        exprecovl->set_tooltip_markup("");
        autocompute->set_tooltip_text("");
        blackEv->set_tooltip_text("");
        whiteEv->set_tooltip_text("");
        sourceGray->set_tooltip_text("");
        sourceabs->set_tooltip_text("");
        targabs->set_tooltip_text("");
        targetGray->set_tooltip_text("");
        baselog->set_tooltip_text("");
        strlog->set_tooltip_text("");
        anglog->set_tooltip_text("");
        detail->set_tooltip_text("");
        Autogray->set_tooltip_text("");
        sensilog->set_tooltip_text("");
        fullimage->set_tooltip_text("");
        ciecam->set_tooltip_text("");
        contl->set_tooltip_text("");
        lightl->set_tooltip_text("");
        lightq->set_tooltip_text("");
        contq->set_tooltip_text("");
        contthres->set_tooltip_text("");
        colorfl->set_tooltip_text("");
        saturl->set_tooltip_text("");
        chroml->set_tooltip_text("");
        catad->set_tooltip_text("");
        expmaskL->set_tooltip_markup("");
        CCmaskshapeL->setTooltip("");
        LLmaskshapeL->setTooltip("");
        HHmaskshapeL->setTooltip("");
        blendmaskL->set_tooltip_text("");
        radmaskL->set_tooltip_text("");
        chromaskL->set_tooltip_text("");
        mask2CurveEditorL->set_tooltip_text("");
        LmaskshapeL->setTooltip("");
        decayl->set_tooltip_text("");
        lowthresl->set_tooltip_text("");
        higthresl->set_tooltip_text("");

    }
}

void LocallabLog::disableListener()
{
    LocallabTool::disableListener();

    autoconn.block(true);
    fullimageConn.block(true);
    ciecamconn.block(true);
    enaLMaskConn.block(true);
    surroundconn.block (true);
    sursourconn.block (true);
    AutograyConn.block(true);
    showmaskLMethodConn.block(true);
}

void LocallabLog::enableListener()
{
    LocallabTool::enableListener();

    autoconn.block(false);
    fullimageConn.block(false);
    ciecamconn.block(false);
    enaLMaskConn.block(false);
    surroundconn.block (false);
    sursourconn.block (false);
    AutograyConn.block(false);
    showmaskLMethodConn.block(false);
}

bool LocallabLog::isMaskViewActive()
{
    return ((showmaskLMethod->get_active_row_number() != 0));
}

void LocallabLog::resetMaskView()
{
    showmaskLMethodConn.block(true);

    showmaskLMethod->set_active(0);

    showmaskLMethodConn.block(false);
}

void LocallabLog::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask)
{
    logMask = showmaskLMethod->get_active_row_number();
}


void LocallabLog::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        exp->set_visible(spot.visilog);
        exp->setEnabled(spot.explog);
        complexity->set_active(spot.complexlog);

        autocompute->set_active(spot.autocompute);
        blackEv->setValue(spot.blackEv);
        repar->setValue(spot.repar);

        whiteEv->setValue(spot.whiteEv);
/*        if(whiteEv->getValue() < 1.5){
            whiteEv->setValue(1.5);
        }
*/
        if (spot.sursour == "Average") {
            sursour->set_active (0);
        } else if (spot.sursour == "Dim") {
            sursour->set_active (1);
        } else if (spot.sursour == "Dark") {
            sursour->set_active (2);
        }


        if (spot.surround == "Average") {
            surround->set_active (0);
        } else if (spot.surround == "Dim") {
            surround->set_active (1);
        } else if (spot.surround == "Dark") {
            surround->set_active (2);
        } else if (spot.surround == "ExtremelyDark") {
            surround->set_active (3);
        }
        recothresl->setValue((double)spot.recothresl);
        lowthresl->setValue((double)spot.lowthresl);
        higthresl->setValue((double)spot.higthresl);
        decayl->setValue((double)spot.decayl);

        ciecam->set_active(spot.ciecam);
        fullimage->set_active(spot.fullimage);
        Autogray->set_active(spot.Autogray);
        sourceGray->setValue(spot.sourceGray);
        sourceabs->setValue(spot.sourceabs);
        catad->setValue(spot.catad);
        saturl->setValue(spot.saturl);
        chroml->setValue(spot.chroml);
        lightl->setValue(spot.lightl);
        lightq->setValue(spot.lightq);
        contl->setValue(spot.contl);
        contthres->setValue(spot.contthres);
        contq->setValue(spot.contq);
        colorfl->setValue(spot.colorfl);
        //LshapeL->setCurve(spot.LcurveL);
        targabs->setValue(spot.targabs);
        targetGray->setValue(spot.targetGray);
        detail->setValue(spot.detail);
        baselog->setValue(spot.baselog);
        sensilog->setValue((double)spot.sensilog);
        strlog->setValue(spot.strlog);
        anglog->setValue(spot.anglog);
        CCmaskshapeL->setCurve(spot.CCmaskcurveL);
        LLmaskshapeL->setCurve(spot.LLmaskcurveL);
        HHmaskshapeL->setCurve(spot.HHmaskcurveL);
        enaLMask->set_active(spot.enaLMask);
        blendmaskL->setValue(spot.blendmaskL);
        radmaskL->setValue(spot.radmaskL);
        chromaskL->setValue(spot.chromaskL);
        LmaskshapeL->setCurve(spot.LmaskcurveL);

        
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Update Log Encoding GUI according to autocompute button state
    updateLogGUI();
    updateLogGUI2();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabLog::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.explog = exp->getEnabled();
        spot.visilog = exp->get_visible();
        spot.complexlog = complexity->get_active_row_number();

        spot.autocompute = autocompute->get_active();
        spot.repar = repar->getValue();
        spot.blackEv = blackEv->getValue();
        spot.whiteEv = whiteEv->getValue();
        spot.fullimage = fullimage->get_active();
        spot.ciecam = ciecam->get_active();
        spot.Autogray = Autogray->get_active();
        spot.sourceGray = sourceGray->getValue();
        spot.sourceabs = sourceabs->getValue();
        spot.targabs = targabs->getValue();
        spot.targetGray = targetGray->getValue();
        spot.catad = catad->getValue();
        spot.saturl = saturl->getValue();
        spot.chroml = chroml->getValue();
        spot.lightl = lightl->getValue();
        spot.lightq = lightq->getValue();
        spot.contl = contl->getValue();
        spot.contthres = contthres->getValue();
        spot.contq = contq->getValue();
        spot.colorfl = colorfl->getValue();
        //spot.LcurveL = LshapeL->getCurve();
        spot.detail = detail->getValue();
        spot.baselog = baselog->getValue();
        spot.sensilog = sensilog->getIntValue();
        spot.strlog = strlog->getValue();
        spot.anglog = anglog->getValue();
        spot.CCmaskcurveL = CCmaskshapeL->getCurve();
        spot.LLmaskcurveL = LLmaskshapeL->getCurve();
        spot.HHmaskcurveL = HHmaskshapeL->getCurve();
        spot.enaLMask = enaLMask->get_active();
        spot.blendmaskL = blendmaskL->getValue();
        spot.radmaskL = radmaskL->getValue();
        spot.chromaskL = chromaskL->getValue();
        spot.LmaskcurveL = LmaskshapeL->getCurve();

        spot.recothresl = recothresl->getValue();
        spot.lowthresl = lowthresl->getValue();
        spot.higthresl = higthresl->getValue();
        spot.decayl = decayl->getValue();

        if (sursour->get_active_row_number() == 0) {
            spot.sursour = "Average";
        } else if (sursour->get_active_row_number() == 1) {
            spot.sursour = "Dim";
        } else if (sursour->get_active_row_number() == 2) {
            spot.sursour = "Dark";
        }

        if (surround->get_active_row_number() == 0) {
            spot.surround = "Average";
        } else if (surround->get_active_row_number() == 1) {
            spot.surround = "Dim";
        } else if (surround->get_active_row_number() == 2) {
            spot.surround = "Dark";
        } else if (surround->get_active_row_number() == 3) {
            spot.surround = "ExtremelyDark";
        }
        
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabLog::enaLMaskChanged()
{
    if (enaLMask->get_active()) {
        maskusablel->show();
        maskunusablel->hide();

    } else {
        maskusablel->hide();
        maskunusablel->show();
    }
    
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaLMask->get_active()) {
                listener->panelChanged(EvLocallabEnaLMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaLMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}



void LocallabLog::updateGUIToMode(const modeType new_type)
{

    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            ciecam->hide();
            ciecam->set_active(false);
            sourceabs->hide();
            targabs->hide();
            saturl->hide();
            chroml->hide();
            contl->hide();
            contthres->hide();
            lightl->hide();
            lightq->hide();
            contq->hide();
            colorfl->hide();
            catad->hide();
            surrHBox->hide();
            expL->hide();
            surHBox->hide();
            expmaskL->hide();
            gradlogFrame->hide();
            exprecovl->hide();
            maskusablel->hide();
            maskunusablel->hide();
            decayl->hide();
            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            ciecam->hide();
            ciecam->set_active(true);

            sourceabs->show();
            targabs->show();
            catad->show();
            saturl->show();
            chroml->show();
            lightl->show();
            lightq->show();
            contl->show();
            contthres->show();
            contq->show();
            colorfl->show();
            surrHBox->show();
            expL->hide();
            surHBox->hide();
            expmaskL->show();
            gradlogFrame->show();
            if (enaLMask->get_active()) {
                maskusablel->show();
                maskunusablel->hide();
                
            } else {
                maskusablel->hide();
                maskunusablel->show();
            }
            
            exprecovl->show();
            decayl->hide();

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            ciecam->hide();
            ciecam->set_active(true);
            sourceabs->show();
            targabs->show();
            catad->show();
            saturl->show();
            chroml->show();
            lightl->show();
            lightq->show();
            contl->show();
            contthres->show();
            contq->show();
            colorfl->show();
            surrHBox->show();
            expL->show();
            expmaskL->show();
            gradlogFrame->show();
            surHBox->show();
            if (enaLMask->get_active()) {
                maskusablel->show();
                maskunusablel->hide();
                
            } else {
                maskusablel->hide();
                maskunusablel->show();
            }
            exprecovl->show();
            decayl->show();

    }
}




void LocallabLog::convertParamToSimple()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();
    ciecam->set_active(false);
    contq->setValue(defSpot.contq);
    contthres->setValue(defSpot.contthres);
    colorfl->setValue(defSpot.colorfl);
    lightl->setValue(defSpot.lightl);
    lightq->setValue(defSpot.lightq);
    sursour->set_active(0);
    strlog->setValue(defSpot.strlog);
    anglog->setValue(defSpot.anglog);
    enaLMask->set_active(false);
    recothresl->setValue(defSpot.recothresl);
    lowthresl->setValue(defSpot.lowthresl);
    higthresl->setValue(defSpot.higthresl);
    decayl->setValue(defSpot.decayl);
    // Enable all listeners
    enableListener();
}


void LocallabLog::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();
    ciecam->set_active(true);
    contq->setValue(defSpot.contq);
    colorfl->setValue(defSpot.colorfl);
    lightl->setValue(defSpot.lightl);
    lightq->setValue(defSpot.lightq);
    sursour->set_active(0);
//    enaLMask->set_active(true);
    decayl->setValue(defSpot.decayl);
    // Enable all listeners
    enableListener();

}



void LocallabLog::showmaskLMethodChanged()
{

    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if(exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
    }
}

void LocallabLog::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == HHmaskshapeL) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskshapeL,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == LLmaskshapeL) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskshapeL,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == CCmaskshapeL) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskshapeL,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == LmaskshapeL) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskshapeL,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        //if (ce == LshapeL) {
        //    if (listener) {
        //        listener->panelChanged(EvlocallabLshapeL,
        //                               M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
        //    }
        //}

    }
}


void LocallabLog::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        // Set default value for adjuster widgets
        repar->setDefault(defSpot.repar);
        blackEv->setDefault(defSpot.blackEv);
        whiteEv->setDefault(defSpot.whiteEv);
        sourceGray->setDefault(defSpot.sourceGray);
        sourceabs->setDefault(defSpot.sourceabs);
        targabs->setDefault(defSpot.targabs);
        targetGray->setDefault(defSpot.targetGray);
        catad->setDefault(defSpot.catad);
        saturl->setDefault(defSpot.saturl);
        chroml->setDefault(defSpot.chroml);
        lightl->setDefault(defSpot.lightl);
        lightq->setDefault(defSpot.lightq);
        contl->setDefault(defSpot.contl);
        contthres->setDefault(defSpot.contthres);
        contq->setDefault(defSpot.contq);
        colorfl->setDefault(defSpot.colorfl);
        detail->setDefault(defSpot.detail);
        baselog->setDefault(defSpot.baselog);
        sensilog->setDefault((double)defSpot.sensilog);
        strlog->setDefault(defSpot.strlog);
        anglog->setDefault(defSpot.anglog);
        blendmaskL->setDefault(defSpot.blendmaskL);
        radmaskL->setDefault(defSpot.radmaskL);
        chromaskL->setDefault(defSpot.chromaskL);
        recothresl->setDefault((double)defSpot.recothresl);
        lowthresl->setDefault((double)defSpot.lowthresl);
        higthresl->setDefault((double)defSpot.higthresl);
        decayl->setDefault((double)defSpot.decayl);


    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabLog::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == repar) {
            if (listener) {
                listener->panelChanged(Evlocallabrepar,
                                       repar->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == blackEv) {
            if (listener) {
                listener->panelChanged(EvlocallabblackEv,
                                       blackEv->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == whiteEv) {
            if (listener) {
                listener->panelChanged(EvlocallabwhiteEv,
                                       whiteEv->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sourceGray) {
            if (listener) {
                listener->panelChanged(EvlocallabsourceGray,
                                       sourceGray->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == sourceabs) {
            if (listener) {
                listener->panelChanged(Evlocallabsourceabs,
                                       sourceabs->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == targabs) {
            if (listener) {
                listener->panelChanged(Evlocallabtargabs,
                                       targabs->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == targetGray) {
            if (listener) {
                listener->panelChanged(EvlocallabtargetGray,
                                       targetGray->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == catad) {
            if (listener) {
                listener->panelChanged(Evlocallabcatad,
                                       catad->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == saturl) {
            if (listener) {
                listener->panelChanged(Evlocallabsaturl,
                                       saturl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chroml) {
            if (listener) {
                listener->panelChanged(Evlocallabchroml,
                                       chroml->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lightl) {
            if (listener) {
                listener->panelChanged(Evlocallablightl,
                                       lightl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lightq) {
            if (listener) {
                listener->panelChanged(Evlocallablightq,
                                       lightq->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }


        if (a == contl) {
            if (listener) {
                listener->panelChanged(Evlocallabcontl,
                                       contl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == contthres) {
            if (listener) {
                listener->panelChanged(Evlocallabcontthres,
                                       contthres->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == contq) {
            if (listener) {
                listener->panelChanged(Evlocallabcontq,
                                       contq->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == colorfl) {
            if (listener) {
                listener->panelChanged(Evlocallabcolorfl,
                                       colorfl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == detail) {
            if (listener) {
                listener->panelChanged(Evlocallabdetail,
                                       detail->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == baselog) {
            if (listener) {
                listener->panelChanged(Evlocallabbaselog,
                                       baselog->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == recothresl) {
            
            if (listener) {
                listener->panelChanged(Evlocallabrecothresl,
                                       recothresl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lowthresl) {
            if (listener) {
                listener->panelChanged(Evlocallablowthresl,
                                       lowthresl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == higthresl) {
            if (listener) {
                listener->panelChanged(Evlocallabhigthresl,
                                       higthresl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == decayl) {
            if (listener) {
                listener->panelChanged(Evlocallabdecayl,
                                       decayl->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }


        if (a == sensilog) {
            if (listener) {
                listener->panelChanged(Evlocallabsensilog,
                                       sensilog->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == strlog) {
            if (listener) {
                listener->panelChanged(Evlocallabstrlog,
                                       strlog->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == anglog) {
            if (listener) {
                listener->panelChanged(Evlocallabanglog,
                                       anglog->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
        
        if (a == blendmaskL) {
            if (listener) {
                listener->panelChanged(EvLocallabblendmaskL,
                                       blendmaskL->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == radmaskL) {
            if (listener) {
                listener->panelChanged(EvLocallabradmaskL,
                                       radmaskL->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromaskL) {
            if (listener) {
                listener->panelChanged(EvLocallabchromaskL,
                                       chromaskL->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        
    }
}

void LocallabLog::updateAutocompute(const float blackev, const float whiteev, const float sourceg, const float sourceab, const float targetg, const float jz1)
{
    if (autocompute->get_active()) {
        idle_register.add(
        [this, blackev, whiteev, sourceg, sourceab, targetg]() -> bool {
            GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

            // Update adjuster values according to autocomputed ones
            disableListener();

            blackEv->setValue(blackev);
            whiteEv->setValue(whiteev);
            sourceGray->setValue(sourceg);
            sourceabs->setValue(sourceab);
            targetGray->setValue(targetg);

            enableListener();

        return false;
        }
        );
    }
}

void LocallabLog::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenalog,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocenalog,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabLog::sursourChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(Evlocallabsursour,
                                   sursour->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}


void LocallabLog::surroundChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(Evlocallabsurround,
                                   surround->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabAuto,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabLog::ciecamChanged()
{
    /*
    if(ciecam->get_active()){
        sourceabs->set_sensitive(true);
        targabs->set_sensitive(true);
        catad->set_sensitive(true);
        surrHBox->set_sensitive(true);
        
        sourceabs->show();
        targabs->show();
        catad->show();
        saturl->show();
        lightl->show();
        contl->show();
        contq->show();
        surrHBox->show();
    } else {
        sourceabs->hide();
        targabs->hide();
        saturl->hide();
        contl->hide();
        lightl->hide();
        contq->hide();
        catad->hide();
        surrHBox->hide();
    }
*/
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (ciecam->get_active()) {
                listener->panelChanged(Evlocallabciecam,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabciecam,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabfullimage,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabLog::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskshapeL->updateLocallabBackground(normChromar);
        LLmaskshapeL->updateLocallabBackground(normLumar);
        HHmaskshapeL->updateLocallabBackground(normHuer);
        LmaskshapeL->updateLocallabBackground(normLumar);

        return false;
    }
    );
}



void LocallabLog::AutograyChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (Autogray->get_active()) {
                listener->panelChanged(EvlocallabAutogray,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvlocallabAutogray,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void LocallabLog::updateLogGUI2()
{
    /*
    if(ciecam->get_active()){
        sourceabs->show();
        targabs->show();
        catad->show();
        saturl->show();
        contl->show();
        lightl->show();
        contq->show();
        surrHBox->show();
    } else {
        sourceabs->hide();
        targabs->hide();
        catad->hide();
        saturl->hide();
        lightl->hide();
        contl->hide();
        contq->hide();
        surrHBox->hide();
    }
    */
}


void LocallabLog::updateLogGUI()
{
    const int mode = complexity->get_active_row_number();
    
    if (autocompute->get_active()) {
        blackEv->set_sensitive(false);
        whiteEv->set_sensitive(false);
        sourceGray->set_sensitive(false);
        if (mode == Expert || mode == Normal) {
            sourceabs->set_sensitive(false);
        } else {
            sourceabs->hide();
        }
    } else {
        blackEv->set_sensitive(true);
        whiteEv->set_sensitive(true);
        sourceGray->set_sensitive(true);
        if (mode == Expert || mode == Normal){
            sourceabs->set_sensitive(true);
        } else {
            sourceabs->hide();
        }
    }
        if (mode == Expert || mode == Normal) { // Keep widget hidden in Simple mode
            exprecovl->show();
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
//    mask_CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASKCOL"))),
    mask_CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    CCmask_shape(static_cast<FlatCurveEditor*>(mask_CurveEditorG->addCurve(CT_Flat, "C", nullptr, false, false))),
    LLmask_shape(static_cast<FlatCurveEditor*>(mask_CurveEditorG->addCurve(CT_Flat, "L", nullptr, false, false))),
    HHmask_shape(static_cast<FlatCurveEditor *>(mask_CurveEditorG->addCurve(CT_Flat, "LC(h)", nullptr, false, true))),
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
    HHhmask_shape(static_cast<FlatCurveEditor *>(mask_HCurveEditorG->addCurve(CT_Flat, "h(h)", nullptr, false, true))),
    mask2CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmask_shape(static_cast<DiagonalCurveEditor*>(mask2CurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    mask2CurveEditorGwav(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVMASK"))),
    LLmask_shapewav(static_cast<FlatCurveEditor*>(mask2CurveEditorGwav->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    csThresholdmask(Gtk::manage(new ThresholdAdjuster(M("TP_LOCALLAB_CSTHRESHOLDBLUR"), 0, 9, 0, 0, 6, 5, 0, false))),
    gradFramemask(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADFRA")))),
    str_mask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -2., 2., 0.05, 0.))),
    ang_mask(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180., 180., 0.1, 0.)))
{
    
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
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
//    LLmask_shapewav->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

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

void LocallabMask::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask)
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
        mask_HCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_HHMASK_TOOLTIP"));
        mask2CurveEditorG->set_tooltip_text(M("TP_LOCALLAB_WAVMASK_TOOLTIP"));
        Lmask_shape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        mask2CurveEditorGwav->set_tooltip_text(M("TP_LOCALLAB_WAVMASK_TOOLTIP"));
        LLmask_shapewav->setTooltip(M("TP_LOCALLAB_LMASK_LEVEL_TOOLTIP"));
        mask_CurveEditorG->set_tooltip_markup(M("TP_LOCALLAB_MASKCURVE_TOOLTIP"));
        strumaskmask->set_tooltip_text(M("TP_LOCALLAB_STRUSTRMASK_TOOLTIP"));
        blurFrame->set_tooltip_text(M("TP_LOCALLAB_BLURMASK_TOOLTIP"));
        toolmask->set_tooltip_text(M("TP_LOCALLAB_TOOLMASK_TOOLTIP"));
        toolmaskFrame->set_tooltip_text(M("TP_LOCALLAB_TOOLCOLFRMASK_TOOLTIP"));
        fftmask->set_tooltip_text(M("TP_LOCALLAB_FFTMASK_TOOLTIP"));
        gammask->set_tooltip_text(M("TP_LOCALLAB_GAMMASK_TOOLTIP"));
        chromask->set_tooltip_text(M("TP_LOCALLAB_CHROMASK_TOOLTIP"));
        slopmask->set_tooltip_text(M("TP_LOCALLAB_SLOMASK_TOOLTIP"));
        shadmask->set_tooltip_text(M("TP_LOCALLAB_SHADMASK_TOOLTIP"));
        contmask->set_tooltip_text(M("TP_LOCALLAB_CONTTHMASK_TOOLTIP"));
        blurmask->set_tooltip_text(M("TP_LOCALLAB_BLURRMASK_TOOLTIP"));
        lapmask->set_tooltip_text(M("TP_LOCALLAB_LAPRAD1_TOOLTIP"));
        csThresholdmask->set_tooltip_text(M("TP_LOCALLAB_WAVEMASK_LEVEL_TOOLTIP"));
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
        mask_HCurveEditorG->set_tooltip_text("");
        mask2CurveEditorG->set_tooltip_text("");
        Lmask_shape->setTooltip("");
        mask2CurveEditorGwav->set_tooltip_text("");
        LLmask_shapewav->setTooltip("");
        mask_CurveEditorG->set_tooltip_markup("");
        strumaskmask->set_tooltip_text("");
        blurFrame->set_tooltip_text("");
        toolmask->set_tooltip_text("");
        toolmaskFrame->set_tooltip_text("");
        fftmask->set_tooltip_text("");
        gammask->set_tooltip_text("");
        chromask->set_tooltip_text("");
        slopmask->set_tooltip_text("");
        shadmask->set_tooltip_text("");
        contmask->set_tooltip_text("");
        blurmask->set_tooltip_text("");
        lapmask->set_tooltip_text("");
        csThresholdmask->set_tooltip_text("");
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
                                       sensimask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == blendmask) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmask,
                                       blendmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == blendmaskab) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskab,
                                       blendmaskab->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == softradiusmask) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusmask,
                                       softradiusmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == strumaskmask) {
            if (listener) {
                listener->panelChanged(Evlocallabstrumaskmask,
                                       strumaskmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == contmask) {
            if (listener) {
                listener->panelChanged(Evlocallabcontmask,
                                       contmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == blurmask) {
            if (listener) {
                listener->panelChanged(Evlocallabblurmask,
                                       blurmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == radmask) {
            if (listener) {
                listener->panelChanged(Evlocallabradmask,
                                       radmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lapmask) {
            if (listener) {
                listener->panelChanged(Evlocallablapmask,
                                       lapmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromask) {
            if (listener) {
                listener->panelChanged(Evlocallabchromask,
                                       chromask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == gammask) {
            if (listener) {
                listener->panelChanged(Evlocallabgammask,
                                       gammask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == slopmask) {
            if (listener) {
                listener->panelChanged(Evlocallabslopmask,
                                       slopmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == shadmask) {
            if (listener) {
                listener->panelChanged(Evlocallabshadmask,
                                       shadmask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == str_mask) {
            if (listener) {
                listener->panelChanged(Evlocallabstr_mask,
                                       str_mask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == ang_mask) {
            if (listener) {
                listener->panelChanged(Evlocallabang_mask,
                                       ang_mask->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       csThresholdmask->getHistoryString() + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == LLmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == HHmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == HHhmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHhmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == Lmask_shape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmask_shape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == LLmask_shapewav) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmask_shapewav,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                   M("TP_LOCALLAB_MODE_SIMPLE") + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    } else if (complexity->get_active_row_number() == Normal) { // New selected mode is Normal one
        // Convert tool widget parameters
        convertParamToNormal();
        // Update GUI based on new mode
        updateGUIToMode(Normal);

        if (listener && isLocActivated) {
            listener->panelChanged(EvlocallabcomplexityWithRefresh,
                                   M("TP_LOCALLAB_MODE_NORMAL") + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    } else if (complexity->get_active_row_number() == Expert) { // New selected mode is Expert one
        // Update GUI based on new mode
        updateGUIToMode(Expert);

        if (listener && isLocActivated) {
            listener->panelChanged(EvlocallabcomplexityWithRefresh,
                                   M("TP_LOCALLAB_MODE_EXPERT") + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}

void LocallabMask::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocena_mask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocena_mask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
    gammask->setValue(defSpot.gammask);
    slopmask->setValue(defSpot.slopmask);
    //Lmask_shape->setCurve(defSpot.Lmask_curve);

    // Enable all listeners
    enableListener();
}

void LocallabMask::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            // Expert and Normal mode widgets are hidden in Simple mode
            softradiusmask->show();
            toolmaskFrame->show();
            struFrame->hide();
            blurFrame->hide();
            gammask->hide();
            slopmask->hide();
            shadmask->hide();
            lapmask->hide();
            mask_HCurveEditorG->hide();
            mask2CurveEditorGwav->hide();
            csThresholdmask->hide();
            gradFramemask->hide();

            break;

        case Normal:
            // Expert mode widgets are hidden in Normal mode
            softradiusmask->show();
            struFrame->hide();
            blurFrame->hide();
            lapmask->hide();
            gammask->show();
            slopmask->show();
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
            softradiusmask->show();
            struFrame->show();
            blurFrame->show();
            toolmaskFrame->show();
            lapmask->show();
            gammask->show();
            slopmask->show();
            shadmask->show();
            mask_HCurveEditorG->show();
            mask2CurveEditorGwav->show();
            csThresholdmask->show();
            gradFramemask->show();
    }
    
}

void LocallabMask::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz)
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

    if(exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
    }
}

void LocallabMask::enamaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enamask->get_active()) {
                listener->panelChanged(EvLocallabEnaMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabtoolmask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabfftmask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
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

/*==== Locallabcie ====*/
Locallabcie::Locallabcie():
    LocallabTool(this, M("TP_LOCALLAB_CIE_TOOLNAME"), M("TP_LOCALLAB_CIE"), false),
    // ciecam specific widgets
    sensicie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
    reparcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGREPART"), 1.0, 100.0, 1., 100.0))),
    jabcie(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_JAB")))),
    modecam(Gtk::manage (new MyComboBoxText ())),
    modecie(Gtk::manage (new MyComboBoxText ())),
    jzFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_JZFRA")))),
    modeHBoxcam(Gtk::manage(new Gtk::Box())),
    modeHBoxcie(Gtk::manage(new Gtk::Box())),
    cieFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOGFRA")))),
    Autograycie(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_AUTOGRAYCIE")))),
    sourceGraycie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOURCE_GRAY"), 1.0, 100.0, 0.1, 18.0))),
    sourceabscie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOURCE_ABS"), 0.01, 16384.0, 0.01, 2000.0))),
    sursourcie(Gtk::manage (new MyComboBoxText ())),
    surHBoxcie(Gtk::manage(new Gtk::Box())),
    cie1Frame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOG1FRA")))),
    cie1lightFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CIELIGHTFRA")))),
    cie1contFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CIECONTFRA")))),
    cie1colorFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CIECOLORFRA")))),
    czlightFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CIELIGHTCONTFRA")))),
    czcolorFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CIECOLORFRA")))),
    PQFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_JZPQFRA")))),
    qtoj(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_JZQTOJ")))),
    lightlcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGLIGHTL"), -100., 100., 0.01, 0.))),
    lightjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZLIGHT"), -100., 100., 0.01, 0.))),
    contjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZCONT"), -100., 100., 0.5, 0.))),
    adapjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZADAP"), 1., 10., 0.05, 4.))),
    jz100(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZ100"), 0.10, 0.90, 0.01, 0.25))),
    pqremap(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZPQREMAP"), 100., 10000., 0.1, 120.))),
    pqremapcam16(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CAM16PQREMAP"), 100., 10000., 0.1, 100.))),
    forcejz(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_JZFORCE")))),
    expjz(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::Box())))),
    jzshFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_JZSHFRA")))),
    hljzcie(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), 0., 100., 1., 0.))),
    hlthjzcie(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_HLTONALW"), 20., 100., 1., 70.))),
    shjzcie(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_SHADOWS"), 0., 100., 1., 0.))),
    shthjzcie(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_SHTONALW"), 20., 100., 1., 40.))),
    radjzcie(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_RADIUS"), 0., 100., 1., 40.))),
    expwavjz(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_JZWAVEXP")))),
    contFramejz(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CONTWFRA")))),
    sigmalcjz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMAWAV"), 0.2, 2.5, 0.01, 1.))),
    LocalcurveEditorwavjz(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAV"))),
    wavshapejz(static_cast<FlatCurveEditor*>(LocalcurveEditorwavjz->addCurve(CT_Flat, "", nullptr, false, false))),
    csThresholdjz(Gtk::manage(new ThresholdAdjuster(M("TP_LOCALLAB_CSTHRESHOLD"), 0, 9, 0, 0, 7, 4, 0, false))),
    clariFramejz(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_CLARIFRA")))),
    clarilresjz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZCLARILRES"), -20., 100., 0.5, 0.))),
    claricresjz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZCLARICRES"), -20., 100., 0.5, 0.))),
    clarisoftjz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.5, 0.))),

    expcam16(Gtk::manage(new MyExpander(false, Gtk::manage(new Gtk::Box())))),
    lightqcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGLIGHTQ"), -100., 100., 0.05, 0.))),
    contlcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCONTL"), -100., 100., 0.5, 0.))),
    contqcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCONQL"), -100., 100., 0.5, 0.))),
    contthrescie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCONTHRES"), -1., 1., 0.01, 0.))),

    logjzFrame(Gtk::manage(new Gtk::Frame())),
    logjz(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_JZLOG")))),
    blackEvjz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLACK_EV"), -16.0, 0.0, 0.1, -5.0))),
    whiteEvjz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WHITE_EV"), 0., 32.0, 0.1, 10.0))),
    targetjz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZTARGET_EV"), 4., 80.0, 0.1, 18.0))),
    bevwevFrame(Gtk::manage(new Gtk::Frame())),
    forcebw(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_BWFORCE")))),

    sigmoidFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SIGFRA")))),
    sigq(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_SIGFRA")))),
    sigmoidldacie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMOIDLAMBDA"), 0.0, 1., 0.01, 0.5))),
    sigmoidthcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMOIDTH"), 0.1, 4., 0.01, 1., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    sigmoidblcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMOIDBL"), 0.5, 1.5, 0.01, 1.))),
    sigmoidqjcie(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_SIGMOIDQJ")))),
    logcie(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_LOGCIE")))),
    sigmoidjzFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SIGJZFRA")))),
    sigjz(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_SIGJZFRA")))),
    sigmoidldajzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMOIDLAMBDA"), 0., 1.0, 0.01, 0.5))),
    sigmoidthjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMOIDTH"), 0.1, 4., 0.01, 1., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    sigmoidbljzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SIGMOIDBL"), 0.5, 1.5, 0.01, 1.))),
    colorflcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCOLORFL"), -100., 100., 0.5, 0.))),
    saturlcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SATURV"), -100., 100., 0.5, 0.))),
    rstprotectcie(Gtk::manage(new Adjuster(M("TP_COLORAPP_RSTPRO"), 0., 100., 0.1, 0.))),
    chromlcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROML"), -100., 100., 0.5, 0.))),
    huecie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_HUECIE"), -100., 100., 0.1, 0.))),
    cieCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVES_CIE"))),
    toneMethodcie(Gtk::manage(new MyComboBoxText())),
    shapecie(static_cast<DiagonalCurveEditor*>(cieCurveEditorG->addCurve(CT_Diagonal, "", toneMethodcie))),
    cieCurveEditorG2(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_COLOR_CIE"))),
    toneMethodcie2(Gtk::manage(new MyComboBoxText())),
    shapecie2(static_cast<DiagonalCurveEditor*>(cieCurveEditorG2->addCurve(CT_Diagonal, "", toneMethodcie2))),

    chromjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZCHROM"), -100., 100., 0.5, 0.))),
    saturjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZSAT"), -100., 100., 0.5, 0.))),
    huejzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZHUECIE"), -100., 100., 0.1, 0.))),
    jz1CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    shapejz(static_cast<DiagonalCurveEditor*>(jz1CurveEditorG->addCurve(CT_Diagonal, "Jz(J)"))),
    shapecz(static_cast<DiagonalCurveEditor*>(jz1CurveEditorG->addCurve(CT_Diagonal, "Cz(C)"))),

    HFramejz(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_JZHFRA")))),
    JzHFramejz(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_JZHJZFRA")))),
    jz2CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    jz3CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    shapeczjz(static_cast<DiagonalCurveEditor*>(jz1CurveEditorG->addCurve(CT_Diagonal, "Cz(J)"))),
    HHshapejz(static_cast<FlatCurveEditor*>(jz3CurveEditorG->addCurve(CT_Flat, "Hz(Hz)", nullptr, false, true))),
    CHshapejz(static_cast<FlatCurveEditor*>(jz3CurveEditorG->addCurve(CT_Flat, "Cz(Hz)", nullptr, false, true))),
    LHshapejz(static_cast<FlatCurveEditor*>(jz2CurveEditorG->addCurve(CT_Flat, "Jz(Hz)", nullptr, false, true))),
    softjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZSOFTCIE"), 0., 100., 0.1, 0.))),
    thrhjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZTHRHCIE"), 40., 150., 0.5, 60.))),
    chjzcie(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_JZCH")))),
    strsoftjzcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_JZSTRSOFTCIE"), 0, 100., 0.5, 100.))),
    
/*
    ciezFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_ZCAMFRA")))),

    lightlzcam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGLIGHTL"), -100., 100., 0.5, 0.))),
    lightqzcam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGLIGHTQ"), -100., 100., 0.05, 0.))),
    contlzcam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCONTL"), -100., 100., 0.5, 0.))),
    contqzcam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCONQL"), -100., 100., 0.5, 0.))),
    contthreszcam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ZCAMTHRES"), 0., 1., 0.01, 0.))),
    colorflzcam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LOGCOLORFL"), -100., 100., 0.5, 0.))),
    saturzcam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SATURV"), -100., 100., 0.5, 0.))),
    chromzcam(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROML"), -100., 100., 0.5, 0.))),
*/    
    expLcie(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_CIETOOLEXP")))),
    cie2Frame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LOG2FRA")))),
    targetGraycie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_TARGET_GRAY"), 5.0, 80.0, 0.1, 18.0))),
    targabscie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOURCE_ABS"), 0.01, 16384.0, 0.01, 16.0))),
    detailcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DETAIL"), 0., 100., 0.1, 0.))),
    catadcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CATAD"), -100., 100., 0.5, 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-orange-small.png"))))),
    surroundcie(Gtk::manage (new MyComboBoxText ())),
    surrHBoxcie(Gtk::manage(new Gtk::Box())),
    exprecovcie(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_DENOI2_EXP")))),
    maskusablecie(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUSABLE")))),
    maskunusablecie(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_MASKUNUSABLE")))),
    recothrescie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKRECOTHRES"), 0., 2., 0.01, 1.))),
    lowthrescie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHRLOW"), 1., 80., 0.5, 12.))),
    higthrescie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKLCTHR"), 20., 99., 0.5, 85.))),
    decaycie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MASKDDECAY"), 0.5, 4., 0.1, 2.))),
    expmaskcie(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWS")))),
    showmaskcieMethod(Gtk::manage(new MyComboBoxText())),
    enacieMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
//    maskSHCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    maskcieCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, "", 1)),
    CCmaskcieshape(static_cast<FlatCurveEditor*>(maskcieCurveEditorG->addCurve(CT_Flat, "C", nullptr, false, false))),
    LLmaskcieshape(static_cast<FlatCurveEditor*>(maskcieCurveEditorG->addCurve(CT_Flat, "L", nullptr, false, false))),
    HHmaskcieshape(static_cast<FlatCurveEditor*>(maskcieCurveEditorG->addCurve(CT_Flat, "LC(h)", nullptr, false, true))),
    blendmaskcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    lapmaskcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskcie(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2cieCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmaskcieshape(static_cast<DiagonalCurveEditor*>(mask2cieCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
   
   
    {
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
    // Parameter Ciecam specific widgets
    const LocallabParams::LocallabSpot defSpot;
    reparcie->setAdjusterListener(this);
    sensicie->setAdjusterListener(this);


    pack_start(*sensicie);
    pack_start(*reparcie);
    modeHBoxcam->set_spacing (2);
    //modeHBoxcam->set_tooltip_markup (M ("TP_LOCALLAB_CAMMODE_TOOLTIP"));
    Gtk::Label* modeLabelcam = Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_CAMMODE") + ":"));
    modeHBoxcam->pack_start (*modeLabelcam, Gtk::PACK_SHRINK);
    modecam->append (M ("TP_LOCALLAB_CAMMODE_CAM16"));
    modecam->append (M ("TP_LOCALLAB_CAMMODE_JZ"));
 //   modecam->append (M ("TP_LOCALLAB_CAMMODE_ALL"));
//    modecam->append (M ("TP_LOCALLAB_CAMMODE_ZCAM"));
    modecam->set_active (0);
    modeHBoxcam->pack_start (*modecam);
    modecamconn = modecam->signal_changed().connect ( sigc::mem_fun (*this, &Locallabcie::modecamChanged) );
    pack_start(*modeHBoxcam);

    modeHBoxcie->set_spacing (2);
    modeHBoxcie->set_tooltip_markup (M ("TP_LOCALLAB_CIEMODE_TOOLTIP"));
    Gtk::Label* modeLabel = Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_CIEMODE") + ":"));
    modeHBoxcie->pack_start (*modeLabel, Gtk::PACK_SHRINK);
    modecie->append (M ("TP_LOCALLAB_CIEMODE_COM"));
    modecie->append (M ("TP_LOCALLAB_CIEMODE_TM"));
    modecie->append (M ("TP_LOCALLAB_CIEMODE_WAV"));
    modecie->append (M ("TP_LOCALLAB_CIEMODE_DR"));
//    modecie->append (M ("TP_LOCALLAB_CIEMODE_LOG"));
    modecie->set_active (0);
    modeHBoxcie->pack_start (*modecie);
    modecieconn = modecie->signal_changed().connect ( sigc::mem_fun (*this, &Locallabcie::modecieChanged) );
    pack_start(*modeHBoxcie);

    surHBoxcie->set_spacing (2);
    surHBoxcie->set_tooltip_markup (M ("TP_LOCALLAB_LOGSURSOUR_TOOLTIP"));
    Gtk::Label* surLabel = Gtk::manage (new Gtk::Label (M ("TP_COLORAPP_SURROUND") + ":"));
    surHBoxcie->pack_start (*surLabel, Gtk::PACK_SHRINK);
    sursourcie->append (M ("TP_COLORAPP_SURROUND_AVER"));
    sursourcie->append (M ("TP_COLORAPP_SURROUND_DIM"));
    sursourcie->append (M ("TP_COLORAPP_SURROUND_DARK"));
    sursourcie->set_active (0);
    surHBoxcie->pack_start (*sursourcie);
    sursourcieconn = sursourcie->signal_changed().connect ( sigc::mem_fun (*this, &Locallabcie::sursourcieChanged) );

    cieFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const cieFBox = Gtk::manage(new ToolParamBlock());
    cieFBox->pack_start(*Autograycie);
    cieFBox->pack_start(*sourceGraycie);
    cieFBox->pack_start(*sourceabscie);
    cieFBox->pack_start(*pqremapcam16);
    PQFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const PQFBox = Gtk::manage(new ToolParamBlock());
    PQFBox->pack_start(*adapjzcie);
    PQFBox->pack_start(*jz100);
    PQFBox->pack_start(*pqremap);
//    PQFBox->pack_start(*forcejz);
//    PQFBox->pack_start(*contthreszcam);
    PQFrame->add(*PQFBox);
    cieFBox->pack_start (*PQFrame);
    logjzFrame->set_label_align(0.025, 0.5);
    logjzFrame->set_label_widget(*logjz);
  //  Gtk::Separator* const separatorjz = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    ToolParamBlock* const logjzBox = Gtk::manage(new ToolParamBlock());
    //logjzBox->pack_start(*blackEvjz);
   // logjzBox->pack_start(*whiteEvjz);
   // logjzBox->pack_start(*separatorjz);
    logjzBox->pack_start(*targetjz);
    logjzFrame->add(*logjzBox);
    cieFBox->pack_start (*logjzFrame);
    bevwevFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const bevwevBox = Gtk::manage(new ToolParamBlock());
    bevwevBox->pack_start(*blackEvjz);
    bevwevBox->pack_start(*whiteEvjz);
    bevwevFrame->add(*bevwevBox);
    cieFBox->pack_start (*bevwevFrame);

    sigmoidFrame->set_label_align(0.025, 0.5);
    sigmoidFrame->set_label_widget(*sigq);
    ToolParamBlock* const sigBox = Gtk::manage(new ToolParamBlock());
    Gtk::Separator* const separatorsig = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    
    sigBox->pack_start(*sigmoidldacie);
    sigBox->pack_start(*sigmoidthcie);
    sigBox->pack_start(*sigmoidblcie);
    sigBox->pack_start(*sigmoidqjcie);
    sigBox->pack_start(*separatorsig);    
    sigBox->pack_start(*logcie);
    sigmoidFrame->add(*sigBox);
    cieFBox->pack_start(*sigmoidFrame);


    sigmoidjzFrame->set_label_align(0.025, 0.5);
    sigmoidjzFrame->set_label_widget(*sigjz);
    ToolParamBlock* const sigjzBox = Gtk::manage(new ToolParamBlock());
    sigjzBox->pack_start(*sigmoidldajzcie);
    sigjzBox->pack_start(*sigmoidthjzcie);
    sigjzBox->pack_start(*sigmoidbljzcie);
    sigjzBox->pack_start(*forcebw);
    sigmoidjzFrame->add(*sigjzBox);
    
  //  jzBox->pack_start(*sigmoidjzFrame);
    cieFBox->pack_start(*sigmoidjzFrame);
    
    cieFBox->pack_start (*surHBoxcie);
    cieFrame->add(*cieFBox);
    pack_start(*cieFrame);

    ToolParamBlock* const jzallBox = Gtk::manage(new ToolParamBlock());
    Gtk::Box *TittleVBoxjz;
    TittleVBoxjz = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    TittleVBoxjz->set_spacing(2);
    Gtk::Box* const LCTitleHBoxjz = Gtk::manage(new Gtk::Box());
    Gtk::Label* const LCLabeljz = Gtk::manage(new Gtk::Label());
    LCLabeljz->set_markup(Glib::ustring("<b>") + (M("TP_LOCALLAB_JZFRA")) + Glib::ustring("</b>"));
    LCTitleHBoxjz->pack_start(*LCLabeljz, Gtk::PACK_SHRINK);
    TittleVBoxjz->pack_start(*LCTitleHBoxjz, Gtk::PACK_SHRINK);
    expjz->setLabel(TittleVBoxjz);

    setExpandAlignProperties(expjz, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    float R, G, B;
    std::vector<GradientMilestone> six_shape;
    for (int i = 0; i < 6; i++) {
        const float x = static_cast<float>(i) * (1.f / 6.f);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        six_shape.emplace_back(x, R, G, B);
    }
    std::vector<GradientMilestone> milestone;
    milestone.push_back ( GradientMilestone (0., 0., 0., 0.) );
    milestone.push_back ( GradientMilestone (1., 1., 1., 1.) );

    jz1CurveEditorG->setCurveListener(this);
    shapejz->setResetCurve(DiagonalCurveType(defSpot.jzcurve.at(0)), defSpot.jzcurve);
    shapejz->setBottomBarBgGradient (milestone);
    shapejz->setLeftBarBgGradient (milestone);

    shapecz->setResetCurve(DiagonalCurveType(defSpot.czcurve.at(0)), defSpot.czcurve);

    std::vector<GradientMilestone> shapeczMilestones;
//    float R, G, B;
    shapecz->setBottomBarColorProvider (this, 1);
    shapecz->setLeftBarColorProvider (this, 1);
    shapecz->setRangeDefaultMilestones (0.05, 0.2, 0.58);

    for (int i = 0; i < 7; i++) {
        float x = float (i) * (1.0f / 6.f);
        Color::hsv2rgb01 (x, 0.5f, 0.5f, R, G, B);
        shapeczMilestones.push_back ( GradientMilestone (double (x), double (R), double (G), double (B)) );
    }

    shapecz->setBottomBarBgGradient (shapeczMilestones);
    shapecz->setLeftBarBgGradient (shapeczMilestones);
    shapecz->setRangeDefaultMilestones (0.05, 0.2, 0.58);

    shapeczjz->setLeftBarColorProvider (this, 1);
    shapeczjz->setRangeDefaultMilestones (0.05, 0.2, 0.58);
    shapeczjz->setResetCurve(DiagonalCurveType(defSpot.czjzcurve.at(0)), defSpot.czjzcurve);
    shapeczjz->setBottomBarBgGradient (milestone);
    shapeczjz->setLeftBarBgGradient (shapeczMilestones);
    shapeczjz->setRangeDefaultMilestones (0.05, 0.2, 0.58);


    jz1CurveEditorG->curveListComplete();
/*
    jz2CurveEditorG->setCurveListener(this);
    shapeczjz->setLeftBarColorProvider (this, 1);
    shapeczjz->setRangeDefaultMilestones (0.05, 0.2, 0.58);
    shapeczjz->setResetCurve(DiagonalCurveType(defSpot.czjzcurve.at(0)), defSpot.czjzcurve);
    shapeczjz->setBottomBarBgGradient (milestone);
    shapeczjz->setLeftBarBgGradient (shapeczMilestones);
    shapeczjz->setRangeDefaultMilestones (0.05, 0.2, 0.58);
    jz2CurveEditorG->curveListComplete();
*/
    jz2CurveEditorG->setCurveListener(this);
    LHshapejz->setIdentityValue(0.);
    LHshapejz->setResetCurve(FlatCurveType(defSpot.LHcurvejz.at(0)), defSpot.LHcurvejz);
  //  LHshapejz->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    LHshapejz->setCurveColorProvider(this, 3);
    LHshapejz->setBottomBarBgGradient(six_shape);
    jz2CurveEditorG->curveListComplete();

    jz3CurveEditorG->setCurveListener(this);

    CHshapejz->setIdentityValue(0.);
    CHshapejz->setResetCurve(FlatCurveType(defSpot.CHcurvejz.at(0)), defSpot.CHcurvejz);
  //  CHshapejz->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    CHshapejz->setCurveColorProvider(this, 3);
    CHshapejz->setBottomBarBgGradient(six_shape);

    HHshapejz->setIdentityValue(0.);
    HHshapejz->setResetCurve(FlatCurveType(defSpot.HHcurvejz.at(0)), defSpot.HHcurvejz);
  //  HHshapejz->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    HHshapejz->setCurveColorProvider(this, 3);
    HHshapejz->setBottomBarBgGradient(six_shape);


    jz3CurveEditorG->curveListComplete();
    
    jzFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const jzBox = Gtk::manage(new ToolParamBlock());
    jzBox->pack_start(*qtoj);
    czlightFrame->set_label_align(0.025, 0.5);
    czcolorFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const ciePzlightBox = Gtk::manage(new ToolParamBlock());
    ciePzlightBox->pack_start(*lightjzcie);
    ciePzlightBox->pack_start(*contjzcie);
    czlightFrame->add(*ciePzlightBox);
    jzBox->pack_start(*czlightFrame);

    ToolParamBlock* const ciePzcolorBox = Gtk::manage(new ToolParamBlock());
    ciePzcolorBox->pack_start(*chromjzcie);
    ciePzcolorBox->pack_start(*saturjzcie);
    ciePzcolorBox->pack_start(*huejzcie);
    czcolorFrame->add(*ciePzcolorBox);
    jzBox->pack_start(*czcolorFrame);
    
    jzBox->pack_start(*jz1CurveEditorG, Gtk::PACK_SHRINK, 4);
    HFramejz->set_label_align(0.025, 0.5);
    ToolParamBlock* const jzHHBox = Gtk::manage(new ToolParamBlock());
    JzHFramejz->set_label_align(0.025, 0.5);
    ToolParamBlock* const jzHBox = Gtk::manage(new ToolParamBlock());

    jzHBox->pack_start(*jz2CurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    jzHBox->pack_start(*thrhjzcie);
    JzHFramejz->add(*jzHBox);
    jzHHBox->pack_start(*JzHFramejz);
    
    jzHHBox->pack_start(*jz3CurveEditorG, Gtk::PACK_SHRINK, 4); //   jzBox->pack_start(*adapjzcie);
    jzHHBox->pack_start(*softjzcie);
    HFramejz->add(*jzHHBox);
    jzBox->pack_start(*HFramejz);
    /*
    sigmoidjzFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const sigjzBox = Gtk::manage(new ToolParamBlock());
    sigjzBox->pack_start(*sigmoidldajzcie);
    sigjzBox->pack_start(*sigmoidthjzcie);
    sigjzBox->pack_start(*sigmoidbljzcie);
    sigjzBox->pack_start(*jabcie);
    sigmoidjzFrame->add(*sigjzBox);
    
  //  jzBox->pack_start(*sigmoidjzFrame);
 */   
    jzshFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const jzshBox = Gtk::manage(new ToolParamBlock());
    jzshBox->pack_start(*hljzcie);
    jzshBox->pack_start(*hlthjzcie);
    jzshBox->pack_start(*shjzcie);
    jzshBox->pack_start(*shthjzcie);
    jzshBox->pack_start(*radjzcie);
    jzshFrame->add(*jzshBox);
    jzBox->pack_start(*jzshFrame);

    setExpandAlignProperties(expwavjz, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    contFramejz->set_label_align(0.025, 0.5);

    sigmalcjz->setAdjusterListener(this);

    LocalcurveEditorwavjz->setCurveListener(this);
    wavshapejz->setIdentityValue(0.);
    wavshapejz->setResetCurve(FlatCurveType(defSpot.locwavcurvejz.at(0)), defSpot.locwavcurvejz);
    LocalcurveEditorwavjz->curveListComplete();
    csThresholdjz->setAdjusterListener(this);

    ToolParamBlock* const coBox2jz = Gtk::manage(new ToolParamBlock());
    coBox2jz->pack_start(*csThresholdjz);
    ToolParamBlock* const coBoxjz = Gtk::manage(new ToolParamBlock());
    coBoxjz->pack_start(*sigmalcjz);
    coBoxjz->pack_start(*LocalcurveEditorwavjz, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    contFramejz->add(*coBoxjz);
    coBox2jz->pack_start(*contFramejz);

    clarilresjz->setAdjusterListener(this);
    claricresjz->setAdjusterListener(this);
    clarisoftjz->setAdjusterListener(this);
    
    clariFramejz->set_label_align(0.025, 0.5);
    ToolParamBlock* const coBox3jz = Gtk::manage(new ToolParamBlock());
    coBox3jz->pack_start(*clarilresjz);
    coBox3jz->pack_start(*claricresjz);
    coBox3jz->pack_start(*clarisoftjz);
    clariFramejz->add(*coBox3jz);
    coBox2jz->pack_start(*clariFramejz);
    expwavjz->add(*coBox2jz, false);
    
    jzBox->pack_start(*expwavjz, false, false);

    jzallBox->add(*jzBox);

    expjz->add(*jzallBox, false);
    
    jabcieConn = jabcie->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::jabcieChanged));
    AutograycieConn = Autograycie->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::AutograycieChanged));
    sigmoidqjcieconn = sigmoidqjcie->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::sigmoidqjcieChanged));
    logcieconn = logcie->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::logcieChanged));
    logjzconn = logjz->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::logjzChanged));
    sigjzconn = sigjz->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::sigjzChanged));
    sigqconn = sigq->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::sigqChanged));
    forcejzConn = forcejz->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::forcejzChanged));
    qtojConn = qtoj->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::qtojChanged));
    chjzcieconn = chjzcie->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::chjzcieChanged));
    forcebwConn = forcebw->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::forcebwChanged));

    sourceGraycie->setAdjusterListener(this);
    sourceGraycie->setLogScale(10, 18, true);

    sourceabscie->setLogScale(500, 0);

    sourceabscie->setAdjusterListener(this);
    setExpandAlignProperties(expLcie, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    saturlcie->setAdjusterListener(this);

    rstprotectcie->setAdjusterListener(this);

    chromlcie->setAdjusterListener(this);
    huecie->setAdjusterListener(this);


    cieCurveEditorG->setCurveListener(this);
    toneMethodcie->append (M ("TP_COLORAPP_TCMODE_LIGHTNESS"));
    toneMethodcie->append (M ("TP_COLORAPP_TCMODE_BRIGHTNESS"));
    toneMethodcie->set_active (0);
    toneMethodcieConn = toneMethodcie->signal_changed().connect(sigc::mem_fun(*this, &Locallabcie::toneMethodcieChanged));
    shapecie->setResetCurve(DiagonalCurveType(defSpot.ciecurve.at(0)), defSpot.ciecurve);
    shapecie->setBottomBarBgGradient (milestone);
    shapecie->setLeftBarBgGradient (milestone);
    cieCurveEditorG->curveListComplete();

    cieCurveEditorG2->setCurveListener(this);
    toneMethodcie2->append (M ("TP_COLORAPP_TCMODE_CHROMA"));
    toneMethodcie2->append (M ("TP_COLORAPP_TCMODE_SATUR"));
    toneMethodcie2->append (M ("TP_COLORAPP_TCMODE_COLORF"));
    toneMethodcie2->set_active (0);
    toneMethodcieConn2 = toneMethodcie2->signal_changed().connect(sigc::mem_fun(*this, &Locallabcie::toneMethodcie2Changed));
    shapecie2->setResetCurve(DiagonalCurveType(defSpot.ciecurve2.at(0)), defSpot.ciecurve2);
    shapecie2->setBottomBarColorProvider (this, 1);
    shapecie2->setLeftBarColorProvider (this, 1);
    shapecie2->setRangeDefaultMilestones (0.05, 0.2, 0.58);
    
    std::vector<GradientMilestone> shape2Milestones;
//    float R, G, B;

    for (int i = 0; i < 7; i++) {
        float x = float (i) * (1.0f / 6.f);
        Color::hsv2rgb01 (x, 0.5f, 0.5f, R, G, B);
        shape2Milestones.push_back ( GradientMilestone (double (x), double (R), double (G), double (B)) );
    }

    shapecie2->setBottomBarBgGradient (shape2Milestones);
    shapecie2->setLeftBarBgGradient (shape2Milestones);

    shapecie2->setRangeDefaultMilestones (0.05, 0.2, 0.58);
    
    cieCurveEditorG2->curveListComplete();


    chromjzcie->setAdjusterListener(this);
    saturjzcie->setAdjusterListener(this);
    huejzcie->setAdjusterListener(this);
    softjzcie->setAdjusterListener(this);
    thrhjzcie->setAdjusterListener(this);
    strsoftjzcie->setAdjusterListener(this);

    lightlcie->setAdjusterListener(this);
    lightjzcie->setAdjusterListener(this);

    lightqcie->setAdjusterListener(this);
    contlcie->setAdjusterListener(this);
    contjzcie->setAdjusterListener(this);
    adapjzcie->setAdjusterListener(this);
    jz100->setAdjusterListener(this);
    pqremap->setAdjusterListener(this);
    pqremapcam16->setAdjusterListener(this);
    hljzcie->setAdjusterListener(this);
    hlthjzcie->setAdjusterListener(this);
    shjzcie->setAdjusterListener(this);
    shthjzcie->setAdjusterListener(this);
    radjzcie->setAdjusterListener(this);
    contthrescie->setAdjusterListener(this);
    blackEvjz->setLogScale(2, -8);
    blackEvjz->setAdjusterListener(this);
    whiteEvjz->setLogScale(16, 0);
    whiteEvjz->setAdjusterListener(this);
    targetjz->setAdjusterListener(this);
    targetjz->setLogScale(10, 18, true);
    sigmoidldacie->setAdjusterListener(this);
    sigmoidthcie->setAdjusterListener(this);
    sigmoidblcie->setAdjusterListener(this);
    sigmoidldajzcie->setAdjusterListener(this);
    sigmoidthjzcie->setAdjusterListener(this);
    sigmoidbljzcie->setAdjusterListener(this);

    contqcie->setAdjusterListener(this);
    colorflcie->setAdjusterListener(this);
/*
    lightlzcam->setAdjusterListener(this);
    lightqzcam->setAdjusterListener(this);
    contlzcam->setAdjusterListener(this);
    contqzcam->setAdjusterListener(this);
    contthreszcam->setAdjusterListener(this);
    colorflzcam->setAdjusterListener(this);
    saturzcam->setAdjusterListener(this);
    chromzcam->setAdjusterListener(this);
*/
    targetGraycie->setAdjusterListener(this);
    targetGraycie->setLogScale(10, 18, true);
    targabscie->setLogScale(500, 0);

    targabscie->setAdjusterListener(this);

    detailcie->setAdjusterListener(this);

    catadcie->setAdjusterListener(this);

    Gtk::Box *TittleVBoxcam16;
    TittleVBoxcam16 = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    TittleVBoxcam16->set_spacing(2);
    Gtk::Box* const LCTitleHBoxcam16 = Gtk::manage(new Gtk::Box());
    Gtk::Label* const LCLabelcam16 = Gtk::manage(new Gtk::Label());
    LCLabelcam16->set_markup(Glib::ustring("<b>") + (M("TP_LOCALLAB_CAM16_FRA")) + Glib::ustring("</b>"));
    LCTitleHBoxcam16->pack_start(*LCLabelcam16, Gtk::PACK_SHRINK);
    TittleVBoxcam16->pack_start(*LCTitleHBoxcam16, Gtk::PACK_SHRINK);
    expcam16->setLabel(TittleVBoxcam16);

    setExpandAlignProperties(expcam16, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    surrHBoxcie->set_spacing (2);
    surrHBoxcie->set_tooltip_markup (M ("TP_COLORAPP_SURROUND_TOOLTIP"));
    Gtk::Label* surrLabelcie = Gtk::manage (new Gtk::Label (M ("TP_COLORAPP_SURROUND") + ":"));
    surrHBoxcie->pack_start (*surrLabelcie, Gtk::PACK_SHRINK);
    surroundcie->append (M ("TP_COLORAPP_SURROUND_AVER"));
    surroundcie->append (M ("TP_COLORAPP_SURROUND_DIM"));
    surroundcie->append (M ("TP_COLORAPP_SURROUND_DARK"));
//    surroundcie->append (M ("TP_COLORAPP_SURROUND_EXDARK"));
    surroundcie->set_active (0);
    surrHBoxcie->pack_start (*surroundcie);
    surroundcieconn = surroundcie->signal_changed().connect ( sigc::mem_fun (*this, &Locallabcie::surroundcieChanged) );

    cie1Frame->set_label_align(0.025, 0.5);
    cie1lightFrame->set_label_align(0.025, 0.5);
    cie1contFrame->set_label_align(0.025, 0.5);
    cie1colorFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const cieP1Box = Gtk::manage(new ToolParamBlock());
    ToolParamBlock* const cieP1lightBox = Gtk::manage(new ToolParamBlock());
    cieP1lightBox->pack_start(*lightlcie);
    cieP1lightBox->pack_start(*lightqcie);
    cie1lightFrame->add(*cieP1lightBox);
    cieP1Box->pack_start(*cie1lightFrame);
    ToolParamBlock* const cieP1contBox = Gtk::manage(new ToolParamBlock());
    cieP1contBox->pack_start(*detailcie);
    cieP1contBox->pack_start(*contlcie);
    cieP1contBox->pack_start(*contqcie);
    cieP1contBox->pack_start(*contthrescie);
    cie1contFrame->add(*cieP1contBox);
    cieP1Box->pack_start(*cie1contFrame);
    ToolParamBlock* const cieP1colorBox = Gtk::manage(new ToolParamBlock());
    cieP1colorBox->pack_start(*chromlcie);
    cieP1colorBox->pack_start(*saturlcie);
    cieP1colorBox->pack_start(*colorflcie);
    cieP1colorBox->pack_start(*huecie);
    cieP1colorBox->pack_start(*rstprotectcie);
    cie1colorFrame->add(*cieP1colorBox);
    cieP1Box->pack_start(*cie1colorFrame);
 //   pack_start(*blackEvjz);
 //   pack_start(*whiteEvjz);
/*
    sigmoidFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const sigBox = Gtk::manage(new ToolParamBlock());
    
    sigBox->pack_start(*sigmoidldacie);
    sigBox->pack_start(*sigmoidthcie);
    sigBox->pack_start(*sigmoidblcie);
    sigBox->pack_start(*sigmoidqjcie);
    sigmoidFrame->add(*sigBox);
    cieP1Box->pack_start(*sigmoidFrame);
    */
    ToolParamBlock* const cieP11Box = Gtk::manage(new ToolParamBlock());
    cieP11Box->pack_start(*cieCurveEditorG);
    cieP11Box->pack_start(*cieCurveEditorG2);
    expLcie->add(*cieP11Box, false);
    cieP1Box->pack_start(*expLcie, false, false);
   // cie1Frame->add(*cieP1Box);
   // expcam16->pack_start(*cie1Frame);
    expcam16->add(*cieP1Box, false);
    
    pack_start(*expcam16, false, false);
    
    pack_start(*expjz, false, false);
/*
    ciezFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const ciePzBox = Gtk::manage(new ToolParamBlock());
    ciePzBox->pack_start(*lightlzcam);
    ciePzBox->pack_start(*lightqzcam);
    ciePzBox->pack_start(*contlzcam);
    ciePzBox->pack_start(*contqzcam);
//    ciePzBox->pack_start(*contthreszcam);
    ciePzBox->pack_start(*colorflzcam);
    ciePzBox->pack_start(*saturzcam);
    ciePzBox->pack_start(*chromzcam);
    ciezFrame->add(*ciePzBox);
    pack_start(*ciezFrame);
*/
    
    cie2Frame->set_label_align(0.025, 0.5);    
    ToolParamBlock* const cieP2Box = Gtk::manage(new ToolParamBlock());
    cieP2Box->pack_start(*targetGraycie);
    cieP2Box->pack_start(*targabscie);
    cieP2Box->pack_start(*catadcie);
    cieP2Box->pack_start(*surrHBoxcie);
//    cieP2Box->pack_start(*detailcie);
//    cieP2Box->pack_start(*jabcie);
    
    cie2Frame->add(*cieP2Box);
    pack_start(*cie2Frame);

    recothrescie->setAdjusterListener(this);
    lowthrescie->setAdjusterListener(this);
    higthrescie->setAdjusterListener(this);
    decaycie->setAdjusterListener(this);
    setExpandAlignProperties(exprecovcie, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);


    setExpandAlignProperties(expmaskcie, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskcieMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcieMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskcieMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskcieMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcieMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskcieMethod->set_active(0);
    showmaskcieMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskcieMethodConn = showmaskcieMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallabcie::showmaskcieMethodChanged));
    
    enacieMaskConn = enacieMask->signal_toggled().connect(sigc::mem_fun(*this, &Locallabcie::enacieMaskChanged));

    maskcieCurveEditorG->setCurveListener(this);
    CCmaskcieshape->setIdentityValue(0.);
    CCmaskcieshape->setResetCurve(FlatCurveType(defSpot.CCmaskciecurve.at(0)), defSpot.CCmaskciecurve);
    CCmaskcieshape->setBottomBarColorProvider(this, 1);

    LLmaskcieshape->setIdentityValue(0.);
    LLmaskcieshape->setResetCurve(FlatCurveType(defSpot.LLmaskciecurve.at(0)), defSpot.LLmaskciecurve);
    LLmaskcieshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskcieshape->setIdentityValue(0.);
    HHmaskcieshape->setResetCurve(FlatCurveType(defSpot.HHmaskciecurve.at(0)), defSpot.HHmaskciecurve);
    HHmaskcieshape->setCurveColorProvider(this, 2);
    HHmaskcieshape->setBottomBarColorProvider(this, 2);

    maskcieCurveEditorG->curveListComplete();
    blendmaskcie->setAdjusterListener(this);

    radmaskcie->setAdjusterListener(this);
    lapmaskcie->setAdjusterListener(this);
    gammaskcie->setAdjusterListener(this);
    slomaskcie->setAdjusterListener(this);

    chromaskcie->setAdjusterListener(this);
    mask2cieCurveEditorG->setCurveListener(this);

    Lmaskcieshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskciecurve.at(0)), defSpot.Lmaskciecurve);
    Lmaskcieshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskcieshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2cieCurveEditorG->curveListComplete();

    ToolParamBlock* const cieBox3 = Gtk::manage(new ToolParamBlock());
    cieBox3->pack_start(*maskusablecie, Gtk::PACK_SHRINK, 0);
    cieBox3->pack_start(*maskunusablecie, Gtk::PACK_SHRINK, 0);
    cieBox3->pack_start(*recothrescie);
    cieBox3->pack_start(*lowthrescie);
    cieBox3->pack_start(*higthrescie);
    cieBox3->pack_start(*decaycie);
    exprecovcie->add(*cieBox3, false);
    pack_start(*exprecovcie, false, false);

    ToolParamBlock* const maskcieBox = Gtk::manage(new ToolParamBlock());
    maskcieBox->pack_start(*showmaskcieMethod, Gtk::PACK_SHRINK, 4);
    maskcieBox->pack_start(*enacieMask, Gtk::PACK_SHRINK, 0);
    maskcieBox->pack_start(*maskcieCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskcieBox->pack_start(*blendmaskcie, Gtk::PACK_SHRINK, 0);
    maskcieBox->pack_start(*radmaskcie, Gtk::PACK_SHRINK, 0);
    maskcieBox->pack_start(*lapmaskcie, Gtk::PACK_SHRINK, 0);
    maskcieBox->pack_start(*chromaskcie, Gtk::PACK_SHRINK, 0);
    maskcieBox->pack_start(*gammaskcie, Gtk::PACK_SHRINK, 0);
    maskcieBox->pack_start(*slomaskcie, Gtk::PACK_SHRINK, 0);
    maskcieBox->pack_start(*mask2cieCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    expmaskcie->add(*maskcieBox, false);
    pack_start(*expmaskcie, false, false);
    
    

    }
Locallabcie::~Locallabcie()
{
    delete jz1CurveEditorG;
    delete jz2CurveEditorG;
    delete jz3CurveEditorG;
    delete cieCurveEditorG;
    delete cieCurveEditorG2;
    delete maskcieCurveEditorG;
    delete mask2cieCurveEditorG;
    delete LocalcurveEditorwavjz;

}

bool Locallabcie::isMaskViewActive()
{
    return ((showmaskcieMethod->get_active_row_number() != 0));
}

void Locallabcie::resetMaskView()
{
    showmaskcieMethodConn.block(true);

    showmaskcieMethod->set_active(0);

    showmaskcieMethodConn.block(false);
}

void Locallabcie::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask, int &logMask, int &maskMask, int &cieMask)
{
   cieMask = showmaskcieMethod->get_active_row_number();
}

void Locallabcie::setDefaultExpanderVisibility()
{
    expLcie->set_expanded(false);
    expjz->set_expanded(false);
    expwavjz->set_expanded(false);
    expcam16->set_expanded(false);
    expmaskcie->set_expanded(false);
    exprecovcie->set_expanded(false);
}
void Locallabcie::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        recothrescie->set_tooltip_text(M("TP_LOCALLAB_RECOTHRES02_TOOLTIP"));
        reparcie->set_tooltip_text(M("TP_LOCALLAB_LOGREPART_TOOLTIP"));
        cieFrame->set_tooltip_text(M("TP_LOCALLAB_LOGSCENE_TOOLTIP"));
        PQFrame->set_tooltip_text(M("TP_LOCALLAB_JZPQFRA_TOOLTIP"));
        qtoj->set_tooltip_text(M("TP_LOCALLAB_JZQTOJ_TOOLTIP"));
        logcie->set_tooltip_text(M("TP_LOCALLAB_LOGCIE_TOOLTIP"));        
        modecam->set_tooltip_text(M("TP_LOCALLAB_JZMODECAM_TOOLTIP"));
        adapjzcie->set_tooltip_text(M("TP_LOCALLAB_JABADAP_TOOLTIP"));
        jz100->set_tooltip_text(M("TP_LOCALLAB_JZ100_TOOLTIP"));
        pqremap->set_tooltip_text(M("TP_LOCALLAB_JZPQREMAP_TOOLTIP"));
        pqremapcam16->set_tooltip_text(M("TP_LOCALLAB_CAM16PQREMAP_TOOLTIP"));
        Autograycie->set_tooltip_text(M("TP_LOCALLAB_LOGAUTOGRAYJZ_TOOLTIP"));
        sigmalcjz->set_tooltip_text(M("TP_LOCALLAB_WAT_SIGMALC_TOOLTIP"));
        logjzFrame->set_tooltip_text(M("TP_LOCALLAB_JZLOGWB_TOOLTIP"));
        blackEvjz->set_tooltip_text(M("TP_LOCALLAB_JZLOGWBS_TOOLTIP"));
        whiteEvjz->set_tooltip_text(M("TP_LOCALLAB_JZLOGWBS_TOOLTIP"));
        clariFramejz->set_tooltip_markup(M("TP_LOCALLAB_CLARIJZ_TOOLTIP"));
        clarilresjz->set_tooltip_text(M("TP_LOCALLAB_WAT_CLARIL_TOOLTIP"));
        claricresjz->set_tooltip_text(M("TP_LOCALLAB_WAT_CLARIC_TOOLTIP"));
        clarisoftjz->set_tooltip_markup(M("TP_LOCALLAB_CLARISOFTJZ_TOOLTIP"));
        wavshapejz->setTooltip(M("TP_LOCALLAB_WAT_WAVSHAPE_TOOLTIP"));
        LocalcurveEditorwavjz->set_tooltip_markup(M("TP_LOCALLAB_WAT_LEVELLOCCONTRAST_TOOLTIP"));
        csThresholdjz->set_tooltip_markup(M("TP_LOCALLAB_WAT_THRESHOLDWAV_TOOLTIP"));
        forcejz->set_tooltip_text(M("TP_LOCALLAB_JZFORCE_TOOLTIP"));
        sourceGraycie->set_tooltip_text(M("TP_LOCALLAB_JZLOGYBOUT_TOOLTIP"));
        sourceabscie->set_tooltip_text(M("TP_COLORAPP_ADAPSCEN_TOOLTIP"));
        cie1Frame->set_tooltip_text(M("TP_LOCALLAB_LOGIMAGE_TOOLTIP"));
        sigmoidFrame->set_tooltip_text(M("TP_LOCALLAB_SIGMOID_TOOLTIP"));
        sigmoidjzFrame->set_tooltip_text(M("TP_LOCALLAB_SIGMOID_TOOLTIP"));
        contlcie->set_tooltip_text(M("TP_LOCALLAB_LOGCONTL_TOOLTIP"));
        contqcie->set_tooltip_text(M("TP_LOCALLAB_LOGCONTQ_TOOLTIP"));
        contthrescie->set_tooltip_text(M("TP_LOCALLAB_LOGCONTTHRES_TOOLTIP"));
        colorflcie->set_tooltip_text(M("TP_LOCALLAB_LOGCOLORF_TOOLTIP"));
        lightlcie->set_tooltip_text(M("TP_LOCALLAB_LOGLIGHTL_TOOLTIP"));
        lightqcie->set_tooltip_text(M("TP_LOCALLAB_LOGLIGHTQ_TOOLTIP"));
        saturlcie->set_tooltip_text(M("TP_LOCALLAB_LOGSATURL_TOOLTIP"));
        rstprotectcie->set_tooltip_text(M("TP_LOCALLAB_RSTPROTECT_TOOLTIP"));
        chromlcie->set_tooltip_text(M("TP_COLORAPP_CHROMA_TOOLTIP"));
        targabscie->set_tooltip_text(M("TP_COLORAPP_VIEWING_ABSOLUTELUMINANCE_TOOLTIP"));
        targetGraycie->set_tooltip_text(M("TP_COLORAPP_YBOUT_TOOLTIP"));
        detailcie->set_tooltip_text(M("TP_LOCALLAB_LOGDETAIL_TOOLTIP"));
        catadcie->set_tooltip_text(M("TP_LOCALLAB_LOGCATAD_TOOLTIP"));
        cie2Frame->set_tooltip_text(M("TP_LOCALLAB_LOGVIEWING_TOOLTIP"));
        sensicie->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        CCmaskcieshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskcieshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskcieshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        expmaskcie->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        blendmaskcie->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        radmaskcie->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        mask2cieCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        Lmaskcieshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        exprecovcie->set_tooltip_markup(M("TP_LOCALLAB_MASKRESH_TOOLTIP"));


    } else {
        reparcie->set_tooltip_text("");
        recothrescie->set_tooltip_text("");
        cieFrame->set_tooltip_text("");
        PQFrame->set_tooltip_text("");
        modecam->set_tooltip_text("");
        qtoj->set_tooltip_text("");
        logcie->set_tooltip_text("");
        jabcie->set_tooltip_text("");
        adapjzcie->set_tooltip_text("");
        jz100->set_tooltip_text("");
        logjzFrame->set_tooltip_text("");
        blackEvjz->set_tooltip_text("");
        whiteEvjz->set_tooltip_text("");
        pqremap->set_tooltip_text("");
        pqremapcam16->set_tooltip_text("");
        Autograycie->set_tooltip_text("");
        forcejz->set_tooltip_text("");
        sourceGraycie->set_tooltip_text("");
        sourceabscie->set_tooltip_text("");
        cie1Frame->set_tooltip_text("");
        sigmoidFrame->set_tooltip_text("");
        sigmoidjzFrame->set_tooltip_text("");
        contlcie->set_tooltip_text("");
        contqcie->set_tooltip_text("");
        contthrescie->set_tooltip_text("");
        colorflcie->set_tooltip_text("");
        lightlcie->set_tooltip_text("");
        lightqcie->set_tooltip_text("");
        saturlcie->set_tooltip_text("");
        rstprotectcie->set_tooltip_text("");
        chromlcie->set_tooltip_text("");
        targabscie->set_tooltip_text("");
        targetGraycie->set_tooltip_text("");
        detailcie->set_tooltip_text("");
        catadcie->set_tooltip_text("");
        cie2Frame->set_tooltip_text("");
        sensicie->set_tooltip_text("");
        CCmaskcieshape->setTooltip("");
        LLmaskcieshape->setTooltip("");
        HHmaskcieshape->setTooltip("");
        expmaskcie->set_tooltip_markup("");
        blendmaskcie->set_tooltip_text("");
        radmaskcie->set_tooltip_text("");
        mask2cieCurveEditorG->set_tooltip_text("");
        Lmaskcieshape->setTooltip("");
        exprecovcie->set_tooltip_markup("");
        sigmalcjz->set_tooltip_text("");
        clarilresjz->set_tooltip_text("");
        claricresjz->set_tooltip_text("");
        clarisoftjz->set_tooltip_markup("");
        clariFramejz->set_tooltip_markup("");
        wavshapejz->setTooltip("");
        LocalcurveEditorwavjz->set_tooltip_markup("");
        csThresholdjz->set_tooltip_markup("");

    }
}
void Locallabcie::disableListener()
{
    LocallabTool::disableListener();
    AutograycieConn.block(true);
    forcejzConn.block(true);
    forcebwConn.block(true);
    qtojConn.block(true);
    jabcieConn.block(true);
    sigmoidqjcieconn.block(true);
    logcieconn.block(true);
    logjzconn.block(true);
    sigjzconn.block(true);
    sigqconn.block(true);
    chjzcieconn.block(true);
    sursourcieconn.block (true);
    surroundcieconn.block (true);
    modecieconn.block (true);
    modecamconn.block (true);
    toneMethodcieConn.block(true);
    toneMethodcieConn2.block(true);
    showmaskcieMethodConn.block(true);
    enacieMaskConn.block(true);
}

void Locallabcie::enableListener()
{
    LocallabTool::enableListener();
    AutograycieConn.block(false);
    forcejzConn.block(false);
    forcebwConn.block(false);
    qtojConn.block(false);
    jabcieConn.block(false);
    sigmoidqjcieconn.block(false);
    logcieconn.block(false);
    logjzconn.block(false);
    sigjzconn.block(false);
    sigqconn.block(false);
    chjzcieconn.block(false);
    sursourcieconn.block (false);
    surroundcieconn.block (false);
    modecieconn.block (false);
    modecamconn.block (false);
    toneMethodcieConn.block(false);
    toneMethodcieConn2.block(false);
    showmaskcieMethodConn.block(false);
    enacieMaskConn.block(false);
}

void Locallabcie::showmaskcieMethodChanged()
{
    // If mask preview is activated, deactivate other Shadow highlight mask preview

    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if(exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskMethod, "");
        }
    }
}



void Locallabcie::enacieMaskChanged()
{
    if (enacieMask->get_active()) {
        maskusablecie->show();
        maskunusablecie->hide();

    } else {
        maskusablecie->hide();
        maskunusablecie->show();
    }
    
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enacieMask->get_active()) {
                listener->panelChanged(EvLocallabEnacieMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocallabEnacieMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        exp->set_visible(spot.visicie);
        exp->setEnabled(spot.expcie);
        complexity->set_active(spot.complexcie);

        reparcie->setValue(spot.reparcie);
        sensicie->setValue(spot.sensicie);

        if (spot.modecam == "cam16") {
            modecam->set_active (0);
        } else if (spot.modecam == "jz") {
            modecam->set_active (1);
//        } else if (spot.modecam == "all") {
//            modecam->set_active (2);
//        } else if (spot.modecam == "zcam") {
//            modecam->set_active (3);
        }

        if (spot.modecie == "com") {
            modecie->set_active (0);
        } else if (spot.modecie == "tm") {
            modecie->set_active (1);
        } else if (spot.modecie == "wav") {
            modecie->set_active (2);
        } else if (spot.modecie == "dr") {
            modecie->set_active (3);
//        } else if (spot.modecie == "log") {
//            modecie->set_active (4);
        }

        if (spot.toneMethodcie == "one") {
            toneMethodcie->set_active(0);
        } else if (spot.toneMethodcie == "two") {
            toneMethodcie->set_active(1);
        }

        if (spot.toneMethodcie2 == "onec") {
            toneMethodcie2->set_active(0);
        } else if (spot.toneMethodcie2 == "twoc") {
            toneMethodcie2->set_active(1);
        } else if (spot.toneMethodcie2 == "thrc") {
            toneMethodcie2->set_active(2);
        }
        
        Autograycie->set_active(spot.Autograycie);
        forcejz->set_active(spot.forcejz);
        forcebw->set_active(spot.forcebw);
        qtoj->set_active(spot.qtoj);
        sourceGraycie->setValue(spot.sourceGraycie);
        sigmoidqjcie->set_active(spot.sigmoidqjcie);
        logcie->set_active(spot.logcie);
        logjz->set_active(spot.logjz);
        sigjz->set_active(spot.sigjz);
        sigq->set_active(spot.sigq);
       // chjzcie->set_active(spot.chjzcie);
        chjzcie->set_active(true);//force to true to avoid other mode
        sourceabscie->setValue(spot.sourceabscie);
        jabcie->set_active(spot.jabcie);
        jabcieChanged();
        modecamChanged();
        
        if(logcie->get_active()) {
            sigmoidldacie->set_sensitive(false);
            sigmoidthcie->set_sensitive(false);
            sigmoidblcie->set_sensitive(false);
            sigmoidqjcie->set_sensitive(false);
        } else {
            sigmoidldacie->set_sensitive(true);
            sigmoidthcie->set_sensitive(true);
            sigmoidblcie->set_sensitive(true);
            sigmoidqjcie->set_sensitive(true);
        }
        
        if (spot.sursourcie == "Average") {
            sursourcie->set_active (0);
        } else if (spot.sursourcie == "Dim") {
            sursourcie->set_active (1);
        } else if (spot.sursourcie == "Dark") {
            sursourcie->set_active (2);
        }

        if (spot.surroundcie == "Average") {
            surroundcie->set_active (0);
        } else if (spot.surroundcie == "Dim") {
            surroundcie->set_active (1);
        } else if (spot.surroundcie == "Dark") {
            surroundcie->set_active (2);
//        } else if (spot.surroundcie == "ExtremelyDark") {
//            surroundcie->set_active (3);
        }
        shapecie->setCurve(spot.ciecurve);
        shapecie2->setCurve(spot.ciecurve2);

        shapejz->setCurve(spot.jzcurve);
        shapecz->setCurve(spot.czcurve);
        shapeczjz->setCurve(spot.czjzcurve);
        HHshapejz->setCurve(spot.HHcurvejz);
        CHshapejz->setCurve(spot.CHcurvejz);
        LHshapejz->setCurve(spot.LHcurvejz);

        saturlcie->setValue(spot.saturlcie);
        rstprotectcie->setValue(spot.rstprotectcie);
        chromlcie->setValue(spot.chromlcie);
        huecie->setValue(spot.huecie);
        chromjzcie->setValue(spot.chromjzcie);
        saturjzcie->setValue(spot.saturjzcie);
        huejzcie->setValue(spot.huejzcie);
        softjzcie->setValue(spot.softjzcie);
        strsoftjzcie->setValue(spot.strsoftjzcie);
        thrhjzcie->setValue(spot.thrhjzcie);
        lightlcie->setValue(spot.lightlcie);
        lightjzcie->setValue(spot.lightjzcie);
        lightqcie->setValue(spot.lightqcie);
        contlcie->setValue(spot.contlcie);
        contjzcie->setValue(spot.contjzcie);
        adapjzcie->setValue(spot.adapjzcie);
        jz100->setValue(spot.jz100);
        pqremap->setValue(spot.pqremap);
        pqremapcam16->setValue(spot.pqremapcam16);
        hljzcie->setValue(spot.hljzcie);
        hlthjzcie->setValue(spot.hlthjzcie);
        shjzcie->setValue(spot.shjzcie);
        shthjzcie->setValue(spot.shthjzcie);
        radjzcie->setValue(spot.radjzcie);
        sigmalcjz->setValue(spot.sigmalcjz);
        wavshapejz->setCurve(spot.locwavcurvejz);
        csThresholdjz->setValue<int>(spot.csthresholdjz);
        clarilresjz->setValue(spot.clarilresjz);
        claricresjz->setValue(spot.claricresjz);
        clarisoftjz->setValue(spot.clarisoftjz);
        contthrescie->setValue(spot.contthrescie);
        blackEvjz->setValue(spot.blackEvjz);
        whiteEvjz->setValue(spot.whiteEvjz);
        targetjz->setValue(spot.targetjz);
        sigmoidldacie->setValue(spot.sigmoidldacie);
        sigmoidthcie->setValue(spot.sigmoidthcie);
        sigmoidblcie->setValue(spot.sigmoidblcie);
        sigmoidldajzcie->setValue(spot.sigmoidldajzcie);
        sigmoidthjzcie->setValue(spot.sigmoidthjzcie);
        sigmoidbljzcie->setValue(spot.sigmoidbljzcie);
        contqcie->setValue(spot.contqcie);
        colorflcie->setValue(spot.colorflcie);
        targabscie->setValue(spot.targabscie);
        targetGraycie->setValue(spot.targetGraycie);
        detailcie->setValue(spot.detailcie);
        catadcie->setValue(spot.catadcie);
/*
        lightlzcam->setValue(spot.lightlzcam);
        lightqzcam->setValue(spot.lightqzcam);
        contlzcam->setValue(spot.contlzcam);
        contqzcam->setValue(spot.contqzcam);
        contthreszcam->setValue(spot.contthreszcam);
        colorflzcam->setValue(spot.colorflzcam);
        saturzcam->setValue(spot.saturzcam);
        chromzcam->setValue(spot.chromzcam);
*/
        enacieMask->set_active(spot.enacieMask);
        CCmaskcieshape->setCurve(spot.CCmaskciecurve);
        LLmaskcieshape->setCurve(spot.LLmaskciecurve);
        HHmaskcieshape->setCurve(spot.HHmaskciecurve);
        blendmaskcie->setValue((double)spot.blendmaskcie);
        radmaskcie->setValue(spot.radmaskcie);
        chromaskcie->setValue(spot.chromaskcie);
        lapmaskcie->setValue(spot.lapmaskcie);
        gammaskcie->setValue(spot.gammaskcie);
        slomaskcie->setValue(spot.slomaskcie);
        Lmaskcieshape->setCurve(spot.Lmaskciecurve);
        recothrescie->setValue((double)spot.recothrescie);
        lowthrescie->setValue((double)spot.lowthrescie);
        higthrescie->setValue((double)spot.higthrescie);
        decaycie->setValue((double)spot.decaycie);

       
    }
    enableListener();
    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));
    // Update Ciecam GUI
    updatecieGUI();
}

void Locallabcie::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
     const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);
        spot.expcie = exp->getEnabled();
        spot.visicie = exp->get_visible();
        spot.complexcie = complexity->get_active_row_number();

        spot.reparcie = reparcie->getValue();
        spot.sensicie = sensicie->getIntValue();

        if (modecam->get_active_row_number() == 0) {
            spot.modecam = "cam16";
        } else if (modecam->get_active_row_number() == 1) {
            spot.modecam = "jz";
//        } else if (modecam->get_active_row_number() == 2) {
//            spot.modecam = "all";
//        } else if (modecam->get_active_row_number() == 3) {
//            spot.modecam = "zcam";
        }

        if (modecie->get_active_row_number() == 0) {
            spot.modecie = "com";
        } else if (modecie->get_active_row_number() == 1) {
            spot.modecie = "tm";
        } else if (modecie->get_active_row_number() == 2) {
            spot.modecie = "wav";
        } else if (modecie->get_active_row_number() == 3) {
            spot.modecie = "dr";
//        } else if (modecie->get_active_row_number() == 4) {
//            spot.modecie = "log";
        }

        if (toneMethodcie->get_active_row_number() == 0) {
            spot.toneMethodcie = "one";
        } else if (toneMethodcie->get_active_row_number() == 1) {
            spot.toneMethodcie = "two";
        }

        if (toneMethodcie2->get_active_row_number() == 0) {
            spot.toneMethodcie2 = "onec";
        } else if (toneMethodcie2->get_active_row_number() == 1) {
            spot.toneMethodcie2 = "twoc";
        } else if (toneMethodcie2->get_active_row_number() == 2) {
            spot.toneMethodcie2 = "thrc";
        }

        spot.Autograycie = Autograycie->get_active();
        spot.forcejz = forcejz->get_active();
        spot.forcebw = forcebw->get_active();
        spot.qtoj = qtoj->get_active();
        spot.jabcie = jabcie->get_active();
        spot.sourceGraycie = sourceGraycie->getValue();
        spot.sourceabscie = sourceabscie->getValue();
        spot.sigmoidqjcie = sigmoidqjcie->get_active();
        spot.logcie = logcie->get_active();
        spot.logjz = logjz->get_active();
        spot.sigjz = sigjz->get_active();
        spot.chjzcie = chjzcie->get_active();
        spot.sigq = sigq->get_active();

        if(sursourcie->get_active_row_number() == 0) {
            spot.sursourcie = "Average";
        } else if (sursourcie->get_active_row_number() == 1) {
            spot.sursourcie = "Dim";
        } else if (sursourcie->get_active_row_number() == 2) {
            spot.sursourcie = "Dark";
        }

        if (surroundcie->get_active_row_number() == 0) {
            spot.surroundcie = "Average";
        } else if (surroundcie->get_active_row_number() == 1) {
            spot.surroundcie = "Dim";
        } else if (surroundcie->get_active_row_number() == 2) {
            spot.surroundcie = "Dark";
//        } else if (surroundcie->get_active_row_number() == 3) {
//            spot.surroundcie = "ExtremelyDark";
        }
        spot.jzcurve = shapejz->getCurve();
        spot.czcurve = shapecz->getCurve();
        spot.czjzcurve = shapeczjz->getCurve();
        spot.HHcurvejz = HHshapejz->getCurve();
        spot.CHcurvejz = CHshapejz->getCurve();
        spot.LHcurvejz = LHshapejz->getCurve();
        spot.ciecurve = shapecie->getCurve();
        spot.ciecurve2 = shapecie2->getCurve();

        spot.saturlcie = saturlcie->getValue();
        spot.rstprotectcie = rstprotectcie->getValue();
        spot.chromlcie = chromlcie->getValue();
        spot.huejzcie = huejzcie->getValue();
        spot.softjzcie = softjzcie->getValue();
        spot.strsoftjzcie = strsoftjzcie->getValue();
        spot.thrhjzcie = thrhjzcie->getValue();
        spot.chromjzcie = chromjzcie->getValue();
        spot.saturjzcie = saturjzcie->getValue();
        spot.huecie = huecie->getValue();
        spot.lightlcie = lightlcie->getValue();
        spot.lightjzcie = lightjzcie->getValue();
        spot.lightqcie = lightqcie->getValue();
        spot.contlcie = contlcie->getValue();
        spot.contjzcie = contjzcie->getValue();
        spot.adapjzcie = adapjzcie->getValue();
        spot.jz100 = jz100->getValue();
        spot.pqremap = pqremap->getValue();
        spot.pqremapcam16 = pqremapcam16->getValue();
        spot.hljzcie = hljzcie->getValue();
        spot.hlthjzcie = hlthjzcie->getValue();
        spot.shjzcie = shjzcie->getValue();
        spot.shthjzcie = shthjzcie->getValue();
        spot.radjzcie = radjzcie->getValue();
        spot.sigmalcjz = sigmalcjz->getValue();
        spot.locwavcurvejz = wavshapejz->getCurve();
        spot.csthresholdjz = csThresholdjz->getValue<int>();
        spot.clarilresjz = clarilresjz->getValue();
        spot.claricresjz = claricresjz->getValue();
        spot.clarisoftjz = clarisoftjz->getValue();
        spot.contthrescie = contthrescie->getValue();
        spot.blackEvjz = blackEvjz->getValue();
        spot.whiteEvjz = whiteEvjz->getValue();
        spot.targetjz = targetjz->getValue();
        spot.sigmoidldacie = sigmoidldacie->getValue();
        spot.sigmoidthcie = sigmoidthcie->getValue();
        spot.sigmoidblcie = sigmoidblcie->getValue();
        spot.sigmoidldajzcie = sigmoidldajzcie->getValue();
        spot.sigmoidthjzcie = sigmoidthjzcie->getValue();
        spot.sigmoidbljzcie = sigmoidbljzcie->getValue();
        spot.contqcie = contqcie->getValue();
        spot.colorflcie = colorflcie->getValue();
        spot.targabscie = targabscie->getValue();
        spot.targetGraycie = targetGraycie->getValue();
        spot.catadcie = catadcie->getValue();
        spot.detailcie = detailcie->getValue();
/*
        spot.lightlzcam = lightlzcam->getValue();
        spot.lightqzcam = lightqzcam->getValue();
        spot.contlzcam = contlzcam->getValue();
        spot.contqzcam = contqzcam->getValue(); 
        spot.contthreszcam = contthreszcam->getValue();
        spot.colorflzcam = colorflzcam->getValue();
        spot.saturzcam = saturzcam->getValue();
        spot.chromzcam = chromzcam->getValue();
*/
        spot.enacieMask = enacieMask->get_active();
        spot.LLmaskciecurve = LLmaskcieshape->getCurve();
        spot.CCmaskciecurve = CCmaskcieshape->getCurve();
        spot.HHmaskciecurve = HHmaskcieshape->getCurve();
        spot.blendmaskcie = blendmaskcie->getIntValue();
        spot.radmaskcie = radmaskcie->getValue();
        spot.chromaskcie = chromaskcie->getValue();
        spot.lapmaskcie = lapmaskcie->getValue();
        spot.gammaskcie = gammaskcie->getValue();
        spot.slomaskcie = slomaskcie->getValue();
        spot.Lmaskciecurve = Lmaskcieshape->getCurve();
        spot.recothrescie = recothrescie->getValue();
        spot.lowthrescie = lowthrescie->getValue();
        spot.higthrescie = higthrescie->getValue();
        spot.decaycie = decaycie->getValue();

    }
}

void Locallabcie::toneMethodcieChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvLocallabtoneMethodcie,
                                   toneMethodcie->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}

void Locallabcie::toneMethodcie2Changed()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvLocallabtoneMethodcie2,
                                   toneMethodcie2->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}


void Locallabcie::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer, const double normHuerjz)
{
    idle_register.add(
    [this, normHuerjz, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        HHshapejz->updateLocallabBackground(normHuerjz);
        CHshapejz->updateLocallabBackground(normHuerjz);
        LHshapejz->updateLocallabBackground(normHuerjz);
        shapejz->updateLocallabBackground(normLumar);
        shapecz->updateLocallabBackground(normChromar);
        shapeczjz->updateLocallabBackground(normLumar);
        shapecie->updateLocallabBackground(normLumar);
        shapecie2->updateLocallabBackground(normChromar);
        CCmaskcieshape->updateLocallabBackground(normChromar);
        LLmaskcieshape->updateLocallabBackground(normLumar);
        HHmaskcieshape->updateLocallabBackground(normHuer);
        Lmaskcieshape->updateLocallabBackground(normLumar);

        return false;
    }
    );
}


void Locallabcie::updateAutocompute(const float blackev, const float whiteev, const float sourceg, const float sourceab, const float targetg, const float jz1)
{    

    if (Autograycie->get_active()) {
        idle_register.add(
        [this, blackev, whiteev, sourceg, sourceab, jz1]() -> bool {
            GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

            // Update adjuster values according to autocomputed ones
            disableListener();
            blackEvjz->setValue(blackev);
            whiteEvjz->setValue(whiteev);
            sourceGraycie->setValue(sourceg);
            sourceabscie->setValue(sourceab);
            float sour = std::min((double) sourceab, 10000.) / 10000.f;
            float pal = std::max(10. * (double) sqrt(sour), 1.5);
            adapjzcie->setValue(pal);//max = 10 and min 1.5
            jz100->setValue(jz1);
            enableListener();

            return false;
        }
        );
    }
}



void Locallabcie::AutograycieChanged()
{

    if (Autograycie->get_active()) {
        sourceGraycie->set_sensitive(false);
        sourceabscie->set_sensitive(false);
        adapjzcie->set_sensitive(false);
        jz100->set_sensitive(false);
        blackEvjz->set_sensitive(false);
        whiteEvjz->set_sensitive(false);
    } else {
        sourceGraycie->set_sensitive(true);
        sourceabscie->set_sensitive(true);
        adapjzcie->set_sensitive(true);
        jz100->set_sensitive(true);
        blackEvjz->set_sensitive(true);
        whiteEvjz->set_sensitive(true);
      //  adapjzcie->set_sensitive(false);
      //  jz100->set_sensitive(false);
    }
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (Autograycie->get_active()) {
                listener->panelChanged(EvlocallabAutograycie,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvlocallabAutograycie,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::forcejzChanged()
{

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (forcejz->get_active()) {
                listener->panelChanged(Evlocallabforcejz,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabforcejz,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::forcebwChanged()
{

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (forcebw->get_active()) {
                listener->panelChanged(Evlocallabforcebw,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabforcebw,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::qtojChanged()
{

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (qtoj->get_active()) {
                listener->panelChanged(Evlocallabqtoj,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabqtoj,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::jabcieChanged()
{ 
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (jabcie->get_active()) {
                listener->panelChanged(Evlocallabjabcie,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabjabcie,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::sigmoidqjcieChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (sigmoidqjcie->get_active()) {
                listener->panelChanged(Evlocallabsigmoidqjcie,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabsigmoidqjcie,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::logcieChanged()
{
    
    if(logcie->get_active()) {
        sigmoidldacie->set_sensitive(false);
        sigmoidthcie->set_sensitive(false);
        sigmoidblcie->set_sensitive(false);
        sigmoidqjcie->set_sensitive(false);
    } else {
        sigmoidldacie->set_sensitive(true);
        sigmoidthcie->set_sensitive(true);
        sigmoidblcie->set_sensitive(true);
        sigmoidqjcie->set_sensitive(true);
    }
    
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (logcie->get_active()) {
                listener->panelChanged(Evlocallablogcie,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallablogcie,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::logjzChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (logjz->get_active()) {
                listener->panelChanged(Evlocallablogjz,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallablogjz,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::sigjzChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (sigjz->get_active()) {
                listener->panelChanged(Evlocallabsigjz,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabsigjz,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::sigqChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (sigq->get_active()) {
                listener->panelChanged(Evlocallabsigq,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabsigq,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::chjzcieChanged()
{
    if (chjzcie->get_active()) {
        thrhjzcie->set_sensitive(true);
    } else {
        thrhjzcie->set_sensitive(false);
    }
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (chjzcie->get_active()) {
                listener->panelChanged(Evlocallabchjzcie,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(Evlocallabchjzcie,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}

void Locallabcie::modecamChanged()
{
    const int mode = complexity->get_active_row_number();
    
    if (modecam->get_active_row_number() == 1 || modecam->get_active_row_number() == 2) {
        expjz->show();
        jzFrame->show();
        adapjzcie->show();
        jz100->show();
        pqremap->show();
        jabcie->show();
        PQFrame->show();
        logjzFrame->show();
        bevwevFrame->show();
        sigmoidjzFrame->show();
        sigmoidFrame->hide();
        forcejz->hide();
        
    } else {
        expjz->hide();
        jzFrame->hide();
        adapjzcie->hide();
        jz100->hide();
        pqremap->hide();
        pqremapcam16->show();
        jabcie->hide();
        PQFrame->hide();
        logjzFrame->hide();
        if (modecam->get_active_row_number() == 0){ 
            bevwevFrame->show();
            sigmoidFrame->show();
        }
        sigmoidjzFrame->hide();
        forcejz->hide();
        catadcie->show();
    }
    surHBoxcie->show();
    cie1Frame->show();
    expcam16->show();
    cie2Frame->show();
    sourceGraycie->show();
    cieFrame->show();

    if (modecam->get_active_row_number() == 1) {
        surHBoxcie->show();
        cie1Frame->hide();
        expcam16->hide();
        targetGraycie->hide();
        targabscie->hide();
        surrHBoxcie->hide();
        forcejz->hide();
        pqremapcam16->hide();
        catadcie->hide();
        cie2Frame->hide();
        exprecovcie->hide();
        expmaskcie->hide();
        if(mode == Expert) {
            exprecovcie->show();
            expmaskcie->show();
        }
        
        }
    if (modecam->get_active_row_number() == 3) {
        if(mode == Expert) {
            cieFrame->show();
            cie1Frame->hide();
            expcam16->hide();
            cie2Frame->show();
            targetGraycie->show();
            targabscie->show();
            surrHBoxcie->show();
            PQFrame->show();
            logjzFrame->show();
            adapjzcie->hide();
            jz100->hide();
            forcejz->hide();
            pqremap->show();
            pqremapcam16->hide();
            catadcie->hide();
            cie2Frame->hide();

        } else {
            cieFrame->hide();
            cie1Frame->hide();
            expcam16->hide();
            cie2Frame->hide();
            catadcie->hide();
            cie2Frame->hide();

        }
    } 

    if(mode != Expert) {
        expjz->hide();
        jzFrame->hide();
        adapjzcie->hide();
        jz100->hide();
        pqremap->show();
        jabcie->hide();
        PQFrame->hide();
        logjzFrame->hide();
        sigmoidjzFrame->hide();
        sigmoidFrame->hide();
        bevwevFrame->hide();
        if (modecam->get_active_row_number() == 0){ 
            bevwevFrame->show();
            sigmoidFrame->show();
        }

        forcejz->hide();
        pqremapcam16->show();
        catadcie->show();
        sourceGraycie->show();

        if (modecam->get_active_row_number() == 1  || modecam->get_active_row_number() == 3) {
            pqremapcam16->hide();
            cieFrame->hide();
            cie1Frame->hide();
            expcam16->hide();
            cie2Frame->hide();
            catadcie->hide();
            cie2Frame->hide();
        }
    } else {
        cieFrame->show();
        cie2Frame->show();
        if (modecam->get_active_row_number() == 0){ 
            bevwevFrame->show();
            sigmoidjzFrame->hide();
            
        }
        if (modecam->get_active_row_number() == 1) {
            targetGraycie->hide();
            targabscie->hide();
            surrHBoxcie->hide();
            forcejz->hide();
            pqremapcam16->hide();
            PQFrame->show();
            logjzFrame->show();
            sigmoidjzFrame->show();
            sigmoidFrame->hide();
            bevwevFrame->show();
            catadcie->hide();
            cie2Frame->hide();
            if (chjzcie->get_active()) {
                thrhjzcie->set_sensitive(true);
            } else {
                thrhjzcie->set_sensitive(false);
            }
            
            
        }
        if (modecam->get_active_row_number() == 3) {
            cieFrame->show();
            cie2Frame->show();
            targetGraycie->show();
            targabscie->show();
            surrHBoxcie->show();
            PQFrame->show();
            logjzFrame->show();
            adapjzcie->hide();
            jz100->hide();
            forcejz->hide();
            pqremap->show();
            pqremapcam16->hide();
            catadcie->hide();
            cie2Frame->hide();
        }
        
    }
        if (modecam->get_active_row_number() == 0 || modecam->get_active_row_number() == 2) {
            targetGraycie->show();
            targabscie->show();
            surrHBoxcie->show();
            cie2Frame->show();
            pqremapcam16->show();
        }
    


    if (isLocActivated && exp->getEnabled()) {

        if (listener) {
            listener->panelChanged(Evlocallabmodecam,
                                   modecam->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}


void Locallabcie::modecieChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        
        const int mode = complexity->get_active_row_number();
                exprecovcie->show();
                expmaskcie->show();
        
        if (modecie->get_active_row_number() > 0) {
            sensicie->hide();
            reparcie->hide();
            exprecovcie->hide();
            expmaskcie->hide();
            
        } else {
            sensicie->show();
            reparcie->show();
            if(mode == Expert) {
                exprecovcie->show();
                expmaskcie->show();
            }
        }
        if (mode == Simple || mode == Normal) { // Keep widget hidden in Normal and Simple mode
        
            modecie->set_active (0);
            sensicie->show();
            reparcie->show();
            
        }

        if (listener) {
            listener->panelChanged(Evlocallabmodecie,
                                   modecie->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}



void Locallabcie::sursourcieChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(Evlocallabsursourcie,
                                   sursourcie->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}

void Locallabcie::surroundcieChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(Evlocallabsurroundcie,
                                   surroundcie->get_active_text() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}

void Locallabcie::updateGUIToMode(const modeType new_type)
{
    switch (new_type) {
        case Simple:
            catadcie->show();
            saturlcie->show();
            rstprotectcie->show();
            chromlcie->hide();
            huecie->hide();
            lightlcie->show();
            lightqcie->hide();
            contlcie->show();
            contthrescie->show();
            contqcie->hide();
            colorflcie->hide();
            surrHBoxcie->show();
            expLcie->hide();
            surHBoxcie->show();
            sourceabscie->show();
            targabscie->show();
            detailcie->hide();
            jabcie->hide();
            modeHBoxcie->hide();
            sensicie->show();
            reparcie->show();
            sigmoidblcie->hide();
            
            expjz->hide();
            jzFrame->hide();
            adapjzcie->hide();
            jz100->hide();
            pqremap->show();
            pqremapcam16->show();
            jabcie->hide();
            targetGraycie->show();
            targabscie->show();
            surrHBoxcie->show();
            forcejz->hide();
            sourceGraycie->show();
            cieFrame->show();
            exprecovcie->hide();
            maskusablecie->hide();
            maskunusablecie->hide();
            decaycie->hide();
            expmaskcie->hide();
            expmaskcie->hide();

            if (modecam->get_active_row_number() == 2) {
                PQFrame->hide();
                logjzFrame->hide();
                sigmoidjzFrame->hide();
                bevwevFrame->hide();
                sigmoidFrame->hide();
            }
            if (modecam->get_active_row_number() == 0){ 
                bevwevFrame->show();
                sigmoidFrame->show();
            }

            if (modecam->get_active_row_number() == 1) {
                cieFrame->hide();
                cie1Frame->hide();
                expcam16->hide();
                forcejz->hide();
                pqremapcam16->hide();
                PQFrame->hide();
                logjzFrame->hide();
                bevwevFrame->hide();
                sigmoidjzFrame->hide();
                sigmoidFrame->hide();
                catadcie->hide();
                cie2Frame->hide();
                maskusablecie->hide();
                maskunusablecie->hide();
                decaycie->hide();
                expmaskcie->hide();
            }
            if (modecam->get_active_row_number() == 3) {
                cieFrame->hide();
                cie1Frame->hide();
                expcam16->hide();
                cie2Frame->hide();
                pqremapcam16->hide();
                PQFrame->hide();
                logjzFrame->hide();
                catadcie->hide();
            }
            
            break;
        case Normal:
            // Expert mode widgets are hidden in Normal mode

            catadcie->show();
            saturlcie->show();
            rstprotectcie->show();
            chromlcie->hide();
            huecie->hide();
            lightlcie->show();
            lightqcie->hide();
            contlcie->show();
            contthrescie->show();
            contqcie->hide();
            colorflcie->hide();
            surrHBoxcie->show();
            expLcie->hide();
            surHBoxcie->show();
            sourceabscie->show();
            targabscie->show();
            detailcie->hide();
            jabcie->hide();
            modeHBoxcie->hide();
            sensicie->show();
            reparcie->show();
            sigmoidblcie->show();
            expjz->hide();
            forcejz->hide();

            jzFrame->hide();
            adapjzcie->hide();
            jz100->hide();
            pqremap->show();
            jabcie->hide();
            targetGraycie->show();
            targabscie->show();
            surrHBoxcie->show();
            pqremapcam16->show();
            sourceGraycie->show();
            cieFrame->show();
            exprecovcie->show();
            expmaskcie->show();
            decaycie->hide();
            lapmaskcie->hide();
            gammaskcie->hide();
            slomaskcie->hide();
            if (enacieMask->get_active()) {
                maskusablecie->show();
                maskunusablecie->hide();
                
            } else {
                maskusablecie->hide();
                maskunusablecie->show();
            }
            if (modecam->get_active_row_number() == 0){ 
                bevwevFrame->show();
                sigmoidFrame->show();
            }

            if (modecam->get_active_row_number() == 2) {
                PQFrame->hide();
                logjzFrame->hide();
                sigmoidjzFrame->hide();
                bevwevFrame->hide();
            }

            if (modecam->get_active_row_number() == 1) {
                cieFrame->hide();
                cie1Frame->hide();
                expcam16->hide();
                forcejz->hide();
                pqremapcam16->hide();
                PQFrame->hide();
                logjzFrame->hide();
                sigmoidjzFrame->hide();
                bevwevFrame->hide();
                sigmoidFrame->hide();
                catadcie->hide();
                cie2Frame->hide();
                exprecovcie->hide();
                expmaskcie->hide();
                maskusablecie->hide();
                maskunusablecie->hide();
                
            }
            if (modecam->get_active_row_number() == 3) {
                cieFrame->hide();
                cie1Frame->hide();
                expcam16->hide();
                cie2Frame->hide();
                pqremapcam16->hide();
                PQFrame->hide();
                catadcie->hide();
                logjzFrame->hide();
            }
            if (modecie->get_active_row_number() > 0) {
                exprecovcie->hide();
                expmaskcie->hide();
            }

            break;

        case Expert:
            // Show widgets hidden in Normal and Simple mode
            catadcie->show();
            saturlcie->show();
            rstprotectcie->show();
            chromlcie->show();
            huecie->show();
            lightlcie->show();
            lightqcie->show();
            contlcie->show();
            contthrescie->show();
            contqcie->show();
            colorflcie->show();
            surrHBoxcie->show();
            expLcie->show();
            surHBoxcie->show();
            sourceabscie->show();
            targabscie->show();
            detailcie->show();
            modeHBoxcie->show();
            sigmoidblcie->show();
            targetGraycie->show();
            targabscie->show();
            surrHBoxcie->show();
            forcejz->hide();
            pqremapcam16->show();
            sourceGraycie->show();
            cieFrame->show();
            exprecovcie->show();
            decaycie->show();
            lapmaskcie->show();
            gammaskcie->show();
            slomaskcie->show();
            expmaskcie->show();
            exprecovcie->show();

            if (enacieMask->get_active()) {
                maskusablecie->show();
                maskunusablecie->hide();
                
            } else {
                maskusablecie->hide();
                maskunusablecie->show();
            }
            if (modecam->get_active_row_number() == 0){ 
                bevwevFrame->show();
            }

            if (modecam->get_active_row_number() == 1 || modecam->get_active_row_number() == 2) {
                jabcie->show();
                expjz->show();
                jzFrame->show();
                adapjzcie->show();
                jz100->show();
                pqremap->show();
                PQFrame->show();
                logjzFrame->show();
                bevwevFrame->show();
                sigmoidjzFrame->show();
                sigmoidFrame->hide();
                forcejz->hide();
                
            }
                cieFrame->show();
                cie2Frame->show();

            if (modecam->get_active_row_number() == 0 || modecam->get_active_row_number() == 2) {
                targetGraycie->show();
                targabscie->show();
                surrHBoxcie->show();
                pqremapcam16->show();
                PQFrame->hide();
                logjzFrame->hide();
                sigmoidjzFrame->hide();
                bevwevFrame->hide();
                if (modecam->get_active_row_number() == 0){ 
                    bevwevFrame->show();
                    sigmoidFrame->show();
                }
                
            }
            if (modecam->get_active_row_number() == 2) {
                PQFrame->show();
                logjzFrame->hide();
                sigmoidjzFrame->hide();
                bevwevFrame->hide();
            }

            if (modecam->get_active_row_number() == 1) {
                surHBoxcie->show();
                targetGraycie->hide();
                targabscie->hide();
                surrHBoxcie->hide();
                pqremapcam16->hide();
                PQFrame->show();
                logjzFrame->show();
                sigmoidjzFrame->show();
                sigmoidFrame->hide();
                bevwevFrame->show();
                catadcie->hide();
                cie2Frame->hide();
                exprecovcie->show();
                expmaskcie->show();
                maskusablecie->show();
                maskunusablecie->show();
                if (chjzcie->get_active()) {
                    thrhjzcie->set_sensitive(true);
                } else {
                    thrhjzcie->set_sensitive(false);
                }
                
           }
                
            if (modecam->get_active_row_number() == 3) {
                cieFrame->show();
                cie1Frame->hide();
                expcam16->hide();
                cie2Frame->show();
                targetGraycie->show();
                targabscie->show();
                surrHBoxcie->show();
                PQFrame->show();
                logjzFrame->show();
                adapjzcie->hide();
                jz100->hide();
                forcejz->hide();
                pqremap->show();
                pqremapcam16->hide();
                catadcie->hide();
            }
            if (modecie->get_active_row_number() > 0) {
                exprecovcie->hide();
                expmaskcie->hide();
            } 

    }
}

void Locallabcie::updatecieGUI()
{
    const int mode = complexity->get_active_row_number();
    expmaskcie->show();
    exprecovcie->show();
    if (modecie->get_active_row_number() > 0) {
        sensicie->hide();
        reparcie->hide();
        exprecovcie->hide();
        expmaskcie->hide();
    } else {
        sensicie->show();
        reparcie->show();
        exprecovcie->show();
        expmaskcie->show();
    }
        surHBoxcie->show();
        cie1Frame->show();
        cie2Frame->show();
        expcam16->show();
    if (modecam->get_active_row_number() == 0){ 
       bevwevFrame->show();
    }
        
    if (modecam->get_active_row_number() == 2  && mode == Expert) {
        PQFrame->show();
        logjzFrame->show();
        sigmoidjzFrame->show();
        bevwevFrame->show();
    }
        sourceGraycie->show();
        cieFrame->show();

    if (enacieMask->get_active() && mode != Simple) {
        maskusablecie->show();
        maskunusablecie->hide();
                
    } else {
        maskusablecie->hide();
        maskunusablecie->show();
    }

    if (modecam->get_active_row_number() == 1) {
        surHBoxcie->show();
        cie1Frame->hide();
        expcam16->hide();
        targetGraycie->hide();
        targabscie->hide();
        surrHBoxcie->hide();
        pqremapcam16->hide();
        PQFrame->show();
        logjzFrame->show();
        sigmoidjzFrame->show();
        bevwevFrame->show();
        sigmoidFrame->hide();
        catadcie->hide();
        cie2Frame->hide();
        if(mode != Expert) {
            cieFrame->hide();
            cie1Frame->hide();
            expcam16->hide();
            cie2Frame->hide();
            PQFrame->hide();
            logjzFrame->hide();
            sigmoidjzFrame->hide();
            sigmoidFrame->hide();
            bevwevFrame->hide();
            if (modecam->get_active_row_number() == 0){ 
                bevwevFrame->show();
                sigmoidFrame->show();
            }
            exprecovcie->hide();
            expmaskcie->hide();
            maskusablecie->hide();
            maskunusablecie->hide();
        }

    }
    if (modecam->get_active_row_number() == 3) {
        if(mode == Expert) {
            cie1Frame->hide();
            expcam16->hide();
            cie2Frame->show();
            targetGraycie->show();
            targabscie->show();
            surrHBoxcie->show();
            cieFrame->show();
            PQFrame->show();
            logjzFrame->show();
            adapjzcie->hide();
            jz100->hide();
            forcejz->hide();
            pqremap->show();
            pqremapcam16->hide();
            PQFrame->show();
            catadcie->hide();
        } else {
            cie1Frame->hide();
            expcam16->hide();
            cie2Frame->hide();
            PQFrame->hide();
            logjzFrame->hide();
        }
    }

    if (Autograycie->get_active()) {
        sourceGraycie->set_sensitive(false);
        sourceabscie->set_sensitive(false);
        adapjzcie->set_sensitive(false);
        jz100->set_sensitive(false);
        blackEvjz->set_sensitive(false);
        whiteEvjz->set_sensitive(false);
    } else {
        sourceGraycie->set_sensitive(true);
        sourceabscie->set_sensitive(true);
        adapjzcie->set_sensitive(true);
        blackEvjz->set_sensitive(true);
        whiteEvjz->set_sensitive(true);
        jz100->set_sensitive(true);
    }

    if (mode == Simple || mode == Normal) { // Keep widget hidden in Normal and Simple mode
        modecie->set_active (0);
        sensicie->show();
        reparcie->show();
    }
    if (modecie->get_active_row_number() > 0) {
        exprecovcie->hide();
        expmaskcie->hide();
    }
    
}


void Locallabcie::convertParamToSimple()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();
    sigmoidblcie->setValue(defSpot.sigmoidblcie);
    showmaskcieMethod->set_active(0);
    enacieMask->set_active(defSpot.enacieMask);
    modecie->set_active(0);
    // Enable all listeners
    enableListener();
}

void Locallabcie::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();
    contqcie->setValue(defSpot.contqcie);
    colorflcie->setValue(defSpot.colorflcie);
    lightqcie->setValue(defSpot.lightqcie);
    chromlcie->setValue(defSpot.chromlcie);
    huecie->setValue(defSpot.huecie);
    detailcie->setValue(defSpot.detailcie);
    jabcie->set_active(defSpot.jabcie);
    LHshapejz->setCurve(defSpot.LHcurvejz);
    CHshapejz->setCurve(defSpot.CHcurvejz);
    HHshapejz->setCurve(defSpot.HHcurvejz);
    shapejz->setCurve(defSpot.jzcurve);
    shapecz->setCurve(defSpot.czcurve);
    shapeczjz->setCurve(defSpot.czjzcurve);
    shapecie->setCurve(defSpot.ciecurve);
    shapecie2->setCurve(defSpot.ciecurve2);
    lightjzcie->setValue(defSpot.lightjzcie);
    contjzcie->setValue(defSpot.contjzcie);
    sigmoidldajzcie->setValue(defSpot.sigmoidldajzcie);
    hljzcie->setValue(defSpot.hljzcie);
    shjzcie->setValue(defSpot.shjzcie);
    chromjzcie->setValue(defSpot.chromjzcie);
    saturjzcie->setValue(defSpot.saturjzcie);
    huejzcie->setValue(defSpot.huejzcie);
    softjzcie->setValue(defSpot.softjzcie);
    strsoftjzcie->setValue(defSpot.strsoftjzcie);
    thrhjzcie->setValue(defSpot.thrhjzcie);
    modecie->set_active(0);
    if (modecam->get_active_row_number() == 1) {
        showmaskcieMethod->set_active(0);
        enacieMask->set_active(defSpot.enacieMask);
        logjz->set_active(defSpot.logjz);
        sigjz->set_active(defSpot.sigjz);
    }
    lapmaskcie->setValue(defSpot.lapmaskcie);
    gammaskcie->setValue(defSpot.gammaskcie);
    slomaskcie->setValue(defSpot.slomaskcie);
    
    // Enable all listeners
    enableListener();
    
}

void Locallabcie::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

        reparcie->setDefault(defSpot.reparcie);
        sensicie->setDefault(defSpot.sensicie);
        sourceGraycie->setDefault(defSpot.sourceGraycie);
        sourceabscie->setDefault(defSpot.sourceabscie);
        saturlcie->setDefault(defSpot.saturlcie);
        rstprotectcie->setDefault(defSpot.rstprotectcie);
        chromlcie->setDefault(defSpot.chromlcie);
        huecie->setDefault(defSpot.huecie);
        chromjzcie->setDefault(defSpot.chromjzcie);
        saturjzcie->setDefault(defSpot.saturjzcie);
        huejzcie->setDefault(defSpot.huejzcie);
        softjzcie->setDefault(defSpot.softjzcie);
        strsoftjzcie->setDefault(defSpot.strsoftjzcie);
        thrhjzcie->setDefault(defSpot.thrhjzcie);
        lightlcie->setDefault(defSpot.lightlcie);
        lightjzcie->setDefault(defSpot.lightjzcie);
        lightqcie->setDefault(defSpot.lightqcie);
        contlcie->setDefault(defSpot.contlcie);
        contjzcie->setDefault(defSpot.contjzcie);
        adapjzcie->setDefault(defSpot.adapjzcie);
        jz100->setDefault(defSpot.jz100);
        pqremap->setDefault(defSpot.pqremap);
        pqremapcam16->setDefault(defSpot.pqremapcam16);
        hljzcie->setDefault(defSpot.hljzcie);
        hlthjzcie->setDefault(defSpot.hlthjzcie);
        shjzcie->setDefault(defSpot.shjzcie);
        shthjzcie->setDefault(defSpot.shthjzcie);
        radjzcie->setDefault(defSpot.radjzcie);
        sigmalcjz->setDefault(defSpot.sigmalcjz);
        csThresholdjz->setDefault<int>(defSpot.csthresholdjz);
        clarilresjz->setDefault(defSpot.clarilresjz);
        claricresjz->setDefault(defSpot.claricresjz);
        clarisoftjz->setDefault(defSpot.clarisoftjz);
        contthrescie->setDefault(defSpot.contthrescie);
        blackEvjz->setDefault(defSpot.blackEvjz);
        whiteEvjz->setDefault(defSpot.whiteEvjz);
        targetjz->setDefault(defSpot.targetjz);
        sigmoidldacie->setDefault(defSpot.sigmoidldacie);
        sigmoidthcie->setDefault(defSpot.sigmoidthcie);
        sigmoidblcie->setDefault(defSpot.sigmoidblcie);
        sigmoidldajzcie->setDefault(defSpot.sigmoidldajzcie);
        sigmoidthjzcie->setDefault(defSpot.sigmoidthjzcie);
        sigmoidbljzcie->setDefault(defSpot.sigmoidbljzcie);
        contqcie->setDefault(defSpot.contqcie);
        colorflcie->setDefault(defSpot.colorflcie);
        targabscie->setDefault(defSpot.targabscie);
        targetGraycie->setDefault(defSpot.targetGraycie);
        catadcie->setDefault(defSpot.catadcie);
        detailcie->setDefault(defSpot.detailcie);
        blendmaskcie->setDefault((double)defSpot.blendmaskcie);
        radmaskcie->setDefault(defSpot.radmaskcie);
        chromaskcie->setDefault(defSpot.chromaskcie);
        lapmaskcie->setDefault(defSpot.lapmaskcie);
        gammaskcie->setDefault(defSpot.gammaskcie);
        slomaskcie->setDefault(defSpot.slomaskcie);
        recothrescie->setDefault((double)defSpot.recothrescie);
        lowthrescie->setDefault((double)defSpot.lowthrescie);
        higthrescie->setDefault((double)defSpot.higthrescie);
        decaycie->setDefault((double)defSpot.decaycie);

    }
}

void Locallabcie::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
       const auto spName = M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")";
       if (ce == shapejz) {
            if (listener) {
                listener->panelChanged(Evlocallabshapejz, spName);
            }
        }

       if (ce == shapecz) {
            if (listener) {
                listener->panelChanged(Evlocallabshapecz, spName);
            }
        }

       if (ce == shapeczjz) {
            if (listener) {
                listener->panelChanged(Evlocallabshapeczjz, spName);
            }
        }
        
       if (ce == HHshapejz) {
            if (listener) {
                listener->panelChanged(EvlocallabHHshapejz, spName);
            }
        }
       
        if (ce == CHshapejz) {
            if (listener) {
                listener->panelChanged(EvlocallabCHshapejz, spName);
            }
        }

        if (ce == LHshapejz) {
            if (listener) {
                listener->panelChanged(EvlocallabLHshapejz, spName);
            }
        }

        if (ce == shapecie) {
            if (listener) {
                listener->panelChanged(Evlocallabshapecie, spName);
            }
        }

        if (ce == shapecie2) {
            if (listener) {
                listener->panelChanged(Evlocallabshapecie2, spName);
            }
        }
        
        if (ce == CCmaskcieshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskcieshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == LLmaskcieshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskcieshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == HHmaskcieshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskcieshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == Lmaskcieshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskcieshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (ce == wavshapejz) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurvejz,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }


    }
}


void Locallabcie::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabcsThresholdjz,
                                   csThresholdjz->getHistoryString() + " (" + escapeHtmlChars(getSpotName()) + ")");
        }
    }
}


void Locallabcie::adjusterChanged(Adjuster* a, double newval)
{
    const LocallabParams::LocallabSpot defSpot;

    if (isLocActivated && exp->getEnabled()) {
        const auto spName = " (" + escapeHtmlChars(getSpotName()) + ")";
        if (a == reparcie) {
            if (listener) {
                listener->panelChanged(Evlocallabreparcie,
                                       reparcie->getTextValue() + spName);
            }
        }

        if (a == sensicie) {
            if (listener) {
                listener->panelChanged(Evlocallabsensicie,
                                       sensicie->getTextValue() + spName);
            }
        }

        if (a == sourceGraycie) {
            if (listener) {
                listener->panelChanged(EvlocallabsourceGraycie,
                                       sourceGraycie->getTextValue() + spName);
            }
        }

        if (a == sourceabscie) {
            float sour = std::min(sourceabscie->getValue(), 10000.) / 10000.f;
            float pal = std::max(10. * (double) sqrt(sour), 1.5);
            adapjzcie->setValue(pal);//max to 10 if La > 10000 and mini to 1.5
            jz100->setValue(defSpot.jz100);
            
            if (listener) {
                listener->panelChanged(Evlocallabsourceabscie,
                                       sourceabscie->getTextValue() + spName );
            }
        }

        if (a == saturlcie) {
            if (listener) {
                listener->panelChanged(Evlocallabsaturlcie,
                                       saturlcie->getTextValue() + spName);
            }
        }

        if (a == rstprotectcie) {
            if (listener) {
                listener->panelChanged(Evlocallabrstprotectcie,
                                       rstprotectcie->getTextValue() + spName);
            }
        }

        if (a == chromlcie) {
            if (listener) {
                listener->panelChanged(Evlocallabchromlcie,
                                       chromlcie->getTextValue() + spName);
            }
        }

        if (a == chromjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabchromjzcie,
                                       chromjzcie->getTextValue() + spName);
            }
        }

        if (a == saturjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabsaturjzcie,
                                       saturjzcie->getTextValue() + spName);
            }
        }

        if (a == huecie) {
            if (listener) {
                listener->panelChanged(Evlocallabhuecie,
                                       huecie->getTextValue() + spName);
            }
        }

        if (a == huejzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabhuejzcie,
                                       huejzcie->getTextValue() + spName);
            }
        }

        if (a == softjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftjzcie,
                                       softjzcie->getTextValue() + spName);
            }
        }

        if (a == strsoftjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabstrsoftjzcie,
                                       strsoftjzcie->getTextValue() + spName);
            }
        }

        if (a == thrhjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabthrhjzcie,
                                       thrhjzcie->getTextValue() + spName);
            }
        }

        if (a == lightlcie) {
            if (listener) {
                listener->panelChanged(Evlocallablightlcie,
                                       lightlcie->getTextValue() + spName);
            }
        }

        if (a == lightjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallablightjzcie,
                                       lightjzcie->getTextValue() + spName);
            }
        }

        if (a == lightqcie) {
            if (listener) {
                listener->panelChanged(Evlocallablightqcie,
                                       lightqcie->getTextValue() + spName);
            }
        }


        if (a == contlcie) {
            if (listener) {
                listener->panelChanged(Evlocallabcontlcie,
                                       contlcie->getTextValue()+ spName);
            }
        }

        if (a == contjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabcontjzcie,
                                       contjzcie->getTextValue() + spName);
            }
        }

        if (a == adapjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabadapjzcie,
                                       adapjzcie->getTextValue() + spName);
            }
        }

        if (a == jz100) {
            if (listener) {
                listener->panelChanged(Evlocallabjz100,
                                       jz100->getTextValue() + spName);
            }
        }

        if (a == pqremap) {
            if (listener) {
                listener->panelChanged(Evlocallabpqremap,
                                       pqremap->getTextValue()+ spName );
            }
        }

        if (a == pqremapcam16) {
            if (listener) {
                listener->panelChanged(Evlocallabpqremapcam16,
                                       pqremapcam16->getTextValue()+ spName );
            }
        }

        if (a == hljzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabhljzcie,
                                       hljzcie->getTextValue() + spName);
            }
        }

        if (a == hlthjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabhlthjzcie,
                                       hlthjzcie->getTextValue() + spName);
            }
        }

        if (a == shjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabshjzcie,
                                       shjzcie->getTextValue()+ spName );
            }
        }

        if (a == shthjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabshthjzcie,
                                       shthjzcie->getTextValue() + spName);
            }
        }

        if (a == radjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabradjzcie,
                                       radjzcie->getTextValue() + spName);
            }
        }

        if (a == sigmalcjz) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmalcjz,
                                       sigmalcjz->getTextValue() + spName);
            }
        }

        if (a == clarilresjz) {
            if (listener) {
                listener->panelChanged(Evlocallabclarilresjz,
                                       clarilresjz->getTextValue() + spName);
            }
        }

        if (a == claricresjz) {
            if (listener) {
                listener->panelChanged(Evlocallabclaricresjz,
                                       claricresjz->getTextValue() + spName);
            }
        }

        if (a == clarisoftjz) {
            if (listener) {
                listener->panelChanged(Evlocallabclarisoftjz,
                                       clarisoftjz->getTextValue() + spName);
            }
        }

        if (a == contthrescie) {
            if (listener) {
                listener->panelChanged(Evlocallabcontthrescie,
                                       contthrescie->getTextValue() + spName);
            }
        }

        if (a == blackEvjz) {
            if (listener) {
                listener->panelChanged(EvlocallabblackEvjz,
                                       blackEvjz->getTextValue() + spName);
            }
        }

        if (a == whiteEvjz) {
            if (listener) {
                listener->panelChanged(EvlocallabwhiteEvjz,
                                       whiteEvjz->getTextValue() + spName);
            }
        }

        if (a == targetjz) {
            if (listener) {
                listener->panelChanged(Evlocallabtargetjz,
                                       targetjz->getTextValue() + spName);
            }
        }

        if (a == sigmoidldacie) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmoidldacie,
                                       sigmoidldacie->getTextValue() + spName);
            }
        }

        if (a == sigmoidldajzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmoidldajzcie,
                                       sigmoidldajzcie->getTextValue() + spName);
            }
        }

        if (a == sigmoidthcie) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmoidthcie,
                                       sigmoidthcie->getTextValue() + spName);
            }
        }

        if (a == sigmoidthjzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmoidthjzcie,
                                       sigmoidthjzcie->getTextValue()+ spName );
            }
        }

        if (a == sigmoidblcie) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmoidblcie,
                                       sigmoidblcie->getTextValue() + spName);
            }
        }

        if (a == sigmoidbljzcie) {
            if (listener) {
                listener->panelChanged(Evlocallabsigmoidbljzcie,
                                       sigmoidbljzcie->getTextValue() + spName);
            }
        }

        if (a == contqcie) {
            if (listener) {
                listener->panelChanged(Evlocallabcontqcie,
                                       contqcie->getTextValue() + spName);
            }
        }

        if (a == colorflcie) {
            if (listener) {
                listener->panelChanged(Evlocallabcolorflcie,
                                       colorflcie->getTextValue()+ spName );
            }
        }

/*
        if (a == lightlzcam) {
            if (listener) {
                listener->panelChanged(Evlocallablightlzcam,
                                       lightlzcam->getTextValue()+ spName );
            }
        }

        if (a == lightqzcam) {
            if (listener) {
                listener->panelChanged(Evlocallablightqzcam,
                                       lightqzcam->getTextValue()+ spName );
            }
        }

        if (a == contlzcam) {
            if (listener) {
                listener->panelChanged(Evlocallabcontlzcam,
                                       contlzcam->getTextValue()+ spName );
            }
        }

        if (a == contqzcam) {
            if (listener) {
                listener->panelChanged(Evlocallabcontqzcam,
                                       contqzcam->getTextValue()+ spName );
            }
        }

        if (a == contthreszcam) {
            if (listener) {
                listener->panelChanged(Evlocallabcontthreszcam,
                                       contthreszcam->getTextValue()+ spName );
            }
        }

        if (a == colorflzcam) {
            if (listener) {
                listener->panelChanged(Evlocallabcolorflzcam,
                                       colorflzcam->getTextValue()+ spName );
            }
        }

        if (a == saturzcam) {
            if (listener) {
                listener->panelChanged(Evlocallabsaturzcam,
                                       saturzcam->getTextValue()+ spName );
            }
        }

        if (a == chromzcam) {
            if (listener) {
                listener->panelChanged(Evlocallabchromzcam,
                                       chromzcam->getTextValue()+ spName );
            }
        }
*/
        if (a == targabscie) {
            if (listener) {
                listener->panelChanged(Evlocallabtargabscie,
                                       targabscie->getTextValue() + spName);
            }
        }

        if (a == targetGraycie) {
            if (listener) {
                listener->panelChanged(EvlocallabtargetGraycie,
                                       targetGraycie->getTextValue() + spName);
            }
        }

        if (a == catadcie) {
            if (listener) {
                listener->panelChanged(Evlocallabcatadcie,
                                       catadcie->getTextValue() + spName);
            }
        }

        if (a == detailcie) {
            if (listener) {
                listener->panelChanged(Evlocallabdetailcie,
                                       detailcie->getTextValue() + spName);
            }
        }
        
        if (a == blendmaskcie) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskcie,
                                       blendmaskcie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == radmaskcie) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskcie,
                                       radmaskcie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == chromaskcie) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskcie,
                                       chromaskcie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lapmaskcie) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskcie,
                                       lapmaskcie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == gammaskcie) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskcie,
                                       gammaskcie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == slomaskcie) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskcie,
                                       slomaskcie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == recothrescie) {
            
            if (listener) {
                listener->panelChanged(Evlocallabrecothrescie,
                                       recothrescie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == lowthrescie) {
            if (listener) {
                listener->panelChanged(Evlocallablowthrescie,
                                       lowthrescie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == higthrescie) {
            if (listener) {
                listener->panelChanged(Evlocallabhigthrescie,
                                       higthrescie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

        if (a == decaycie) {
            if (listener) {
                listener->panelChanged(Evlocallabdecaycie,
                                       decaycie->getTextValue() + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }

    }
}

void Locallabcie::enabledChanged()
{
     if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenacie,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            } else {
                listener->panelChanged(EvLocenacie,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(getSpotName()) + ")");
            }
        }
    }
}
