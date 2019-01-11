/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
 *  2017 Jacques Desmis <jdesmis@gmail.com>
 *  2018 Pierre Cabrera <pierre.cab@gmail.com>
 */

#include "locallab.h"
#include "rtimage.h"
#include <iomanip>
#include "../rtengine/rt_math.h"
#include "options.h"
#include <cmath>
#include "edit.h"
#include "guiutils.h"
#include <string>
#include <unistd.h>
#include "../rtengine/improcfun.h"

#define MINCHRO 0.
#define MAXCHRO 150
#define MAXCHROCC 100

using namespace rtengine;

extern Options options;

Locallab::Locallab():
    FoldableToolPanel(this, "locallab", M("TP_LOCALLAB_LABEL"), false, true),

    // Expander widgets
    expsettings(new ControlSpotPanel()),
    expcolor(new MyExpander(true, M("TP_LOCALLAB_COFR"))),
    expexpose(new MyExpander(true, M("TP_LOCALLAB_EXPOSE"))),
    expvibrance(new MyExpander(true, M("TP_LOCALLAB_VIBRANCE"))),
    expsoft(new MyExpander(true, M("TP_LOCALLAB_SOFT"))),
    explabregion(new MyExpander(true, M("TP_LOCALLAB_LABREGION"))),
    expblur(new MyExpander(true, M("TP_LOCALLAB_BLUFR"))),
    exptonemap(new MyExpander(true, M("TP_LOCALLAB_TM"))),
    expreti(new MyExpander(true, M("TP_LOCALLAB_RETI"))),
    expsharp(new MyExpander(true, M("TP_LOCALLAB_SHARP"))),
    expcontrast(new MyExpander(true, M("TP_LOCALLAB_LOC_CONTRAST"))),
    expcbdl(new MyExpander(true, M("TP_LOCALLAB_CBDL"))),
    expdenoi(new MyExpander(true, M("TP_LOCALLAB_DENOIS"))),

    // CurveEditorGroup widgets
    // Color & Light
    llCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_LUM"))),
    maskCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    // Exposure
    curveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVEEDITOR_TONES_LABEL"))),
    maskexpCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    // Vibrance
    curveEditorGG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_LABEL"))),
    // Retinex
    LocalcurveEditorgainT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAIN"))),

    // Adjuster widgets
    // Color & Light
    lightness(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LIGHTNESS"), -100, 100, 1, 0))),
    contrast(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0))),
    chroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMA"), -100, 150, 1, 0))),
    sensi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    // Exposure
    expcomp(Gtk::manage(new Adjuster(M("TP_EXPOSURE_EXPCOMP"), -200, 200, 5, 0))),
    hlcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 60))),
    hlcomprthresh(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 33))),
    black(Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0))),
    shcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50))),
    warm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WARM"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-orange-small.png"))))),
    sensiex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    // Vibrance
    saturated(Gtk::manage(new Adjuster(M("TP_VIBRANCE_SATURATED"), -100., 100., 1., 0.))),
    pastels(Gtk::manage(new Adjuster(M("TP_VIBRANCE_PASTELS"), -100., 100., 1., 0.))),
    sensiv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    //Soft Light
    streng(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENG"), 1, 100, 1, 1))),
    sensisf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 1, 100, 1, 19))),
    // Blur & Noise
    radius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADIUS"), 1, 100, 1, 1))),
    strength(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENGTH"), 0, 100, 1, 0))),
    sensibn(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIBN"), 0, 100, 1, 40))),
    // Tone Mapping
    stren(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STREN"), -50, 100, 1, 0))),
    gamma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAM"), 80, 150, 1, 100))),
    estop(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ESTOP"), 10, 400, 1, 140))),
    scaltm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALTM"), 1, 100, 1, 10))),
    rewei(Gtk::manage(new Adjuster(M("TP_LOCALLAB_REWEI"), 0, 9, 1, 0))),
    sensitm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    // Retinex
    str(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STR"), 0, 100, 1, 0))),
    chrrt(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHRRT"), 0, 100, 1, 0))),
    neigh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NEIGH"), 14, 150, 1, 50))),
    vart(Gtk::manage(new Adjuster(M("TP_LOCALLAB_VART"), 50, 500, 1, 200))),
    dehaz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DEHAZ"), 0, 100, 1, 0))),
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIH"), 0, 100, 1, 19))),
    // Sharpening
    sharcontrast(Gtk::manage(new Adjuster(M("TP_SHARPENING_CONTRAST"), 0, 200, 1, 20))),
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 42, 500, 1, 4))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 100))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 75))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sharblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARBLUR"), 20, 200, 1, 20))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    // Local Contrast
    lcradius(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_RADIUS"), 20, 200, 1, 80))),
    lcamount(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_AMOUNT"), 0, 100, 1, 0))),
    lcdarkness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_DARKNESS"), 0, 300, 1, 100))),
    lclightness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_LIGHTNESS"), 0, 300, 1, 100))),
    sensilc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    // Contrast by detail levels
    chromacbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMACBDL"), 0, 300, 1, 0))),
    threshold(Gtk::manage(new Adjuster(M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 100, 1, 20))),
    sensicb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSICB"), 0, 100, 1, 19))),
    // Denoise
    noiselumf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINE"), MINCHRO, MAXCHRO, 1, 0))),
    noiselumc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMCOARSE"), MINCHRO, MAXCHROCC, 1, 0))),
    noiselumdetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMDETAIL"), 0, 100, 1, 0))),
    noiselequal(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELEQUAL"), -2, 10, 1, 7, Gtk::manage(new RTImage("circle-white-small.png")), Gtk::manage(new RTImage("circle-black-small.png"))))),
    noisechrof(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROFINE"), MINCHRO, MAXCHRO, 1, 0))),
    noisechroc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROCOARSE"), MINCHRO, MAXCHROCC, 1, 0))),
    noisechrodetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHRODETAIL"), 0, 100, 1, 0))),
    adjblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ADJ"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-red-small.png"))))),
    bilateral(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BILATERAL"), 0, 100, 1, 0))),
    sensiden(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIDEN"), 0, 100, 1, 30))),

    // ButtonCheck widgets
    // Color & Light
    curvactiv(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CURV")))),
    invers(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    // Vibrance
    protectSkins(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PROTECTSKINS")))),
    avoidColorShift(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_AVOIDCOLORSHIFT")))),
    pastSatTog(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PASTSATTOG")))),
    // Blur & Noise
    activlum(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ACTIV")))),
    // Retinex
    inversret(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    // Sharpening
    inverssha(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    // Others
    avoid(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_AVOID")))),

    // ComboBox widgets
    // Color & Light
    qualitycurveMethod(Gtk::manage(new MyComboBoxText())),
    showmaskcolMethod(Gtk::manage(new MyComboBoxText())),
    //Exposure
    showmaskexpMethod(Gtk::manage(new MyComboBoxText())),
    // Blur & Noise
    blurMethod(Gtk::manage(new MyComboBoxText())),
    // Retinex
    retinexMethod(Gtk::manage(new MyComboBoxText())),

    // ThresholdAdjuster widgets
    // Vibrance
    psThreshold(Gtk::manage(new ThresholdAdjuster(M("TP_VIBRANCE_PSTHRESHOLD"), -100., 100., 0., M("TP_VIBRANCE_PSTHRESHOLD_WEIGTHING"), 0, 0., 100., 75., M("TP_VIBRANCE_PSTHRESHOLD_SATTHRESH"), 0, this, false))),

    // Other widgets
    labqualcurv(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_QUALCURV_METHOD") + ":"))),
    lumacontrastMinusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")))),
    lumaneutralButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")))),
    lumacontrastPlusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS")))),

    // Others
    defparams(nullptr),
    defpedited(nullptr),
    pe(nullptr),
    nexthuer(0.),
    nextlumar(0.),
    nextchromar(0.)

{
    ToolVBox* const panel = Gtk::manage(new ToolVBox());

    CurveListener::setMulti(true);
    float R, G, B;

    LocallabParams::LocallabSpot defSpot;

    // Settings
    expsettings->getExpander()->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsettings->getExpander()));
    expsettings->setLevel(2);

    panel->pack_start(*expsettings->getExpander(), false, false);

    // Color & Light
    expcolor->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcolor));
    enablecolorConn = expcolor->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcolor));
    expcolor->set_tooltip_text(M("TP_LOCALLAB_EXPCOLOR_TOOLTIP"));

    curvactivConn = curvactiv->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::curvactivChanged));

    lightness->setAdjusterListener(this);

    contrast->setAdjusterListener(this);

    chroma->setAdjusterListener(this);

    sensi->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensi->setAdjusterListener(this);

    qualitycurveMethod->append(M("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVENH"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVENHSU"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVENCONTRAST"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVENSOB2"));

    qualitycurveMethod->set_active(0);
    qualitycurveMethod->set_tooltip_markup(M("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));
    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::qualitycurveMethodChanged));

    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskcolMethod->append(M("TP_LOCALLAB_USEMASK"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMASK"));

    showmaskcolMethod->set_active(0);
    showmaskcolMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskcolMethodConn  = showmaskcolMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::showmaskcolMethodChanged));
    
    llCurveEditorG->setCurveListener(this);

    llshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "L(L)"));
    llshape->setResetCurve(DiagonalCurveType(defSpot.llcurve.at(0)), defSpot.llcurve);
    llshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    std::vector<GradientMilestone> mllshape;
    mllshape.push_back(GradientMilestone(0., 0., 0., 0.));
    mllshape.push_back(GradientMilestone(1., 1., 1., 1.));
    llshape->setBottomBarBgGradient(mllshape);
    llshape->setLeftBarBgGradient(mllshape);

    ccshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "C(C)"));
    ccshape->setResetCurve(DiagonalCurveType(defSpot.cccurve.at(0)), defSpot.cccurve);
    ccshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    std::vector<GradientMilestone> mccshape;
    mccshape.push_back(GradientMilestone(0., 0., 0., 0.));
    mccshape.push_back(GradientMilestone(1., 1., 1., 1.));
    ccshape->setBottomBarBgGradient(mccshape);
    ccshape->setLeftBarBgGradient(mccshape);

    llCurveEditorG->newLine();

    LHshape = static_cast<FlatCurveEditor*>(llCurveEditorG->addCurve(CT_Flat, "L(H)", nullptr, false, true));
    LHshape->setIdentityValue(0.);
    LHshape->setResetCurve(FlatCurveType(defSpot.LHcurve.at(0)), defSpot.LHcurve);
    LHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    LHshape->setCurveColorProvider(this, 1);
    std::vector<GradientMilestone> mLHshape;

    for (int i = 0; i < 7; i++) {
        float x = float (i) * (1.0f / 6.0);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        mLHshape.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    LHshape->setBottomBarBgGradient(mLHshape);

    HHshape = static_cast<FlatCurveEditor*>(llCurveEditorG->addCurve(CT_Flat, "H(H)", nullptr, false, true));
    HHshape->setIdentityValue(0.);
    HHshape->setResetCurve(FlatCurveType(defSpot.HHcurve.at(0)), defSpot.HHcurve);
    HHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    HHshape->setCurveColorProvider(this, 1);
    std::vector<GradientMilestone> mHHshape;

    for (int i = 0; i < 7; i++) {
        float x = float (i) * (1.0f / 6.0);

        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        mHHshape.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    HHshape->setBottomBarBgGradient(mHHshape);

    llCurveEditorG->curveListComplete();
    maskCurveEditorG->setCurveListener(this);

    inversConn  = invers->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversChanged));
    transLabels = Gtk::manage (new Gtk::Label ("---"));

    CCmaskshape = static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false));
    CCmaskshape->setIdentityValue(0.);
    CCmaskshape->setResetCurve(FlatCurveType(defSpot.CCmaskcurve.at(0)), defSpot.CCmaskcurve);
    CCmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    CCmaskshape->setBottomBarColorProvider(this, 7);

    LLmaskshape = static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false));
    LLmaskshape->setIdentityValue(0.);
    LLmaskshape->setResetCurve(FlatCurveType(defSpot.LLmaskcurve.at(0)), defSpot.LLmaskcurve);
    LLmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    LLmaskshape->setBottomBarBgGradient(mllshape);
    

    maskCurveEditorG->curveListComplete();
    
    ToolParamBlock* const colorBox = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const superFrame = Gtk::manage(new Gtk::Frame());
    superFrame->set_label_align(0.025, 0.5);
    superFrame->set_label_widget(*curvactiv);
    ToolParamBlock* const superBox = Gtk::manage(new ToolParamBlock());
    superBox->pack_start(*lightness);
    superBox->pack_start(*contrast);
    superBox->pack_start(*chroma);
    superFrame->add(*superBox);
    colorBox->pack_start(*superFrame);
    colorBox->pack_start(*sensi);
    Gtk::HBox* const qualcurvbox = Gtk::manage(new Gtk::HBox());
    qualcurvbox->pack_start(*labqualcurv, Gtk::PACK_SHRINK, 4);
    qualcurvbox->pack_start(*qualitycurveMethod);
    colorBox->pack_start(*qualcurvbox);
    colorBox->pack_start(*llCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    colorBox->pack_start(*invers);
    maskcolFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHOW")));
    maskcolFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const maskcolBox = Gtk::manage(new ToolParamBlock());
    maskcolBox->pack_start(*transLabels, Gtk::PACK_SHRINK, 4);
    maskcolBox->pack_start(*showmaskcolMethod, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*maskCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor

    maskcolFrame->add(*maskcolBox);
    colorBox->pack_start(*maskcolFrame);

    expcolor->add(*colorBox);
    expcolor->setLevel(2);

    panel->pack_start(*expcolor, false, false);

    // Exposure
    expexpose->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expexpose));
    enableexposeConn = expexpose->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expexpose));

    expcomp->setAdjusterListener(this);

    hlcompr->setAdjusterListener(this);

    hlcomprthresh->setAdjusterListener(this);

    black->setAdjusterListener(this);

    shcompr->setAdjusterListener(this);

    warm->set_tooltip_text(M("TP_LOCALLAB_WARM_TOOLTIP"));
    warm->setAdjusterListener(this);

    sensiex->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensiex->setAdjusterListener(this);

    curveEditorG->setCurveListener(this);

    shapeexpos = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, ""));
    shapeexpos->setResetCurve(DiagonalCurveType(defSpot.excurve.at(0)), defSpot.excurve);
    shapeexpos->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_TONES_TOOLTIP"));
    std::vector<GradientMilestone> mshapeexpos;
    mshapeexpos.push_back(GradientMilestone(0., 0., 0., 0.));
    mshapeexpos.push_back(GradientMilestone(1., 1., 1., 1.));
    shapeexpos->setBottomBarBgGradient(mshapeexpos);
    shapeexpos->setLeftBarBgGradient(mshapeexpos);

    curveEditorG->curveListComplete();
    maskexpCurveEditorG->setCurveListener(this);
    transLabels2 = Gtk::manage (new Gtk::Label ("---"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskexpMethod->append(M("TP_LOCALLAB_USEMASK"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMASK"));

    showmaskexpMethod->set_active(0);
    showmaskexpMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskexpMethodConn  = showmaskexpMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::showmaskexpMethodChanged));

    CCmaskexpshape = static_cast<FlatCurveEditor*>(maskexpCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false));
    CCmaskexpshape->setIdentityValue(0.);
    CCmaskexpshape->setResetCurve(FlatCurveType(defSpot.CCmaskexpcurve.at(0)), defSpot.CCmaskexpcurve);
    CCmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    CCmaskexpshape->setBottomBarColorProvider(this, 7);
    const ColorToningParams default_params;

    LLmaskexpshape = static_cast<FlatCurveEditor*>(maskexpCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false));
    LLmaskexpshape->setIdentityValue(0.);
    LLmaskexpshape->setResetCurve(FlatCurveType(defSpot.LLmaskexpcurve.at(0)), defSpot.LLmaskexpcurve);   
    LLmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    LLmaskexpshape->setBottomBarBgGradient(mllshape);

    maskexpCurveEditorG->curveListComplete();
    
    
    maskexpFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHOW")));
    maskexpFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const maskexpBox = Gtk::manage(new ToolParamBlock());
    maskexpBox->pack_start(*transLabels2, Gtk::PACK_SHRINK, 4);
    maskexpBox->pack_start(*showmaskexpMethod, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*maskexpCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor

    maskexpFrame->add(*maskexpBox);
    
    ToolParamBlock* const exposeBox = Gtk::manage(new ToolParamBlock());
    exposeBox->pack_start(*expcomp);
    exposeBox->pack_start(*hlcompr);
    exposeBox->pack_start(*hlcomprthresh);
    exposeBox->pack_start(*black);
    exposeBox->pack_start(*shcompr);
    exposeBox->pack_start(*warm);
    exposeBox->pack_start(*sensiex);
    exposeBox->pack_start(*curveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    exposeBox->pack_start(*maskexpFrame);
    
    expexpose->add(*exposeBox);
    expexpose->setLevel(2);

    panel->pack_start(*expexpose, false, false);

    // Vibrance
    expvibrance->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expvibrance));
    enablevibranceConn = expvibrance->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expvibrance));

    saturated->setAdjusterListener(this);

    pastels->setAdjusterListener(this);

    psThreshold->set_tooltip_markup(M("TP_VIBRANCE_PSTHRESHOLD_TOOLTIP"));
    psThreshold->setAdjusterListener(this);

    pskinsconn = protectSkins->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::protectskins_toggled));

    ashiftconn = avoidColorShift->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::avoidcolorshift_toggled));

    pastsattogconn = pastSatTog->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::pastsattog_toggled));

    sensiv->setAdjusterListener(this);

    curveEditorGG->setCurveListener(this);

    skinTonesCurve = static_cast<DiagonalCurveEditor*>(curveEditorGG->addCurve(CT_Diagonal, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES")));
    skinTonesCurve->setTooltip(M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_TOOLTIP"));
    std::vector<GradientMilestone> mskinTonesCurve;
    // -0.1 rad < Hue < 1.6 rad
    Color::hsv2rgb01(0.92f, 0.45f, 0.6f, R, G, B);
    mskinTonesCurve.push_back(GradientMilestone(0.0, double (R), double (G), double (B)));
    Color::hsv2rgb01(0.14056f, 0.45f, 0.6f, R, G, B);
    mskinTonesCurve.push_back(GradientMilestone(1.0, double (R), double (G), double (B)));
    skinTonesCurve->setBottomBarBgGradient(mskinTonesCurve);
    skinTonesCurve->setLeftBarBgGradient(mskinTonesCurve);
    skinTonesCurve->setRangeLabels(
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE1"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE2"),
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE3"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE4")
    );
    skinTonesCurve->setRangeDefaultMilestones(0.1, 0.4, 0.85);

    curveEditorGG->curveListComplete();

    ToolParamBlock* const vibranceBox = Gtk::manage(new ToolParamBlock());
    vibranceBox->pack_start(*saturated, Gtk::PACK_SHRINK, 0);
    vibranceBox->pack_start(*pastels, Gtk::PACK_SHRINK, 0);
    vibranceBox->pack_start(*psThreshold, Gtk::PACK_SHRINK, 0);
    vibranceBox->pack_start(*protectSkins, Gtk::PACK_SHRINK, 0);
    vibranceBox->pack_start(*avoidColorShift, Gtk::PACK_SHRINK, 0);
    vibranceBox->pack_start(*pastSatTog, Gtk::PACK_SHRINK, 0);
    vibranceBox->pack_start(*sensiv, Gtk::PACK_SHRINK, 0);
    vibranceBox->pack_start(*curveEditorGG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    expvibrance->add(*vibranceBox);
    expvibrance->setLevel(2);

    panel->pack_start(*expvibrance, false, false);

    // Soft Light
    expsoft->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsoft));
    enablesoftConn = expsoft->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expsoft));

    streng->setAdjusterListener(this);

    sensisf->setAdjusterListener(this);

    ToolParamBlock* const softBox = Gtk::manage(new ToolParamBlock());
    softBox->pack_start(*streng);
    softBox->pack_start(*sensisf);
    expsoft->add(*softBox);
    expsoft->setLevel(2);

    panel->pack_start(*expsoft, false, false);

    // Lab Region
    explabregion->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), explabregion));
    enablelabregionConn = explabregion->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), explabregion));

    ToolParamBlock* const labRegionBox = Gtk::manage(new ToolParamBlock());
    explabregion->add(*labRegionBox);
    explabregion->setLevel(2);
/*    
    labRegionSlope = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_SLOPE"), 0.1, 4.0, 0.001, 1));
    labRegionSlope->setLogScale(4, 0.1);
    labRegionSlope->setAdjusterListener(this);
    labRegionBox->pack_start(*labRegionSlope);
*/
    // panel->pack_start(*explabregion, false, false);
//    labRegionSlope->delay = options.adjusterMaxDelay;

    // Blur & Noise
    expblur->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expblur));
    enableblurConn = expblur->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expblur));

    radius->setAdjusterListener(this);

    strength->setAdjusterListener(this);

    sensibn->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensibn->setAdjusterListener(this);

    blurMethod->append(M("TP_LOCALLAB_BLNORM"));
    blurMethod->append(M("TP_LOCALLAB_BLINV"));
    blurMethod->append(M("TP_LOCALLAB_BLSYM"));
    blurMethod->set_active(0);
    blurMethod->set_tooltip_markup(M("TP_LOCALLAB_BLMETHOD_TOOLTIP"));
    blurMethodConn = blurMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::blurMethodChanged));

    activlumConn  = activlum->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::activlumChanged));

    ToolParamBlock* const blurrBox = Gtk::manage(new ToolParamBlock());
    blurrBox->pack_start(*radius);
    blurrBox->pack_start(*strength);
    blurrBox->pack_start(*sensibn);
    blurrBox->pack_start(*blurMethod);
    blurrBox->pack_start(*activlum);
    expblur->add(*blurrBox);
    expblur->setLevel(2);

    panel->pack_start(*expblur, false, false);

    // Tone Mapping
    exptonemap->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), exptonemap));
    enabletonemapConn = exptonemap->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), exptonemap));

    stren->setAdjusterListener(this);

    gamma->setAdjusterListener(this);

    estop->setAdjusterListener(this);

    scaltm->setAdjusterListener(this);

    rewei->setAdjusterListener(this);

    sensitm->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensitm->setAdjusterListener(this);

    ToolParamBlock* const tmBox = Gtk::manage(new ToolParamBlock());
    tmBox->pack_start(*stren);
    tmBox->pack_start(*gamma);
    tmBox->pack_start(*estop);
    tmBox->pack_start(*scaltm);
    tmBox->pack_start(*rewei);
    tmBox->pack_start(*sensitm);
    exptonemap->add(*tmBox);
    exptonemap->setLevel(2);

    // panel->pack_start(*exptonemap, false, false);

    // Retinex
    expreti->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expreti));
    enableretiConn = expreti->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expreti));

    retinexMethod->append(M("TP_RETINEX_LOW"));
    retinexMethod->append(M("TP_RETINEX_UNIFORM"));
    retinexMethod->append(M("TP_RETINEX_HIGH"));
    retinexMethod->set_active(0);
    retinexMethod->set_tooltip_markup(M("TP_LOCRETI_METHOD_TOOLTIP"));
    retinexMethodConn = retinexMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::retinexMethodChanged));

    str->setAdjusterListener(this);

    neigh->setAdjusterListener(this);

    vart->setAdjusterListener(this);

    dehaz->setAdjusterListener(this);

    chrrt->setAdjusterListener(this);

    sensih->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensih->setAdjusterListener(this);

    LocalcurveEditorgainT->setCurveListener(this);

    cTgainshape = static_cast<FlatCurveEditor*>(LocalcurveEditorgainT->addCurve(CT_Flat, "", nullptr, false, false));
    cTgainshape->setIdentityValue(0.);
    cTgainshape->setResetCurve(FlatCurveType(defSpot.localTgaincurve.at(0)), defSpot.localTgaincurve);
    cTgainshape->setTooltip(M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));

    LocalcurveEditorgainT->curveListComplete();

    inversretConn  = inversret->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversretChanged));

    ToolParamBlock* const retiBox = Gtk::manage(new ToolParamBlock());
    retiBox->pack_start(*retinexMethod);
    retiBox->pack_start(*str);
    retiBox->pack_start(*chrrt);
    retiBox->pack_start(*neigh);
    retiBox->pack_start(*vart);
    retiBox->pack_start(*dehaz);
    retiBox->pack_start(*sensih);
    retiBox->pack_start(*LocalcurveEditorgainT, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    retiBox->pack_start(*inversret);
    expreti->add(*retiBox);
    expreti->setLevel(2);

    panel->pack_start(*expreti, false, false);

    // Sharpening
    expsharp->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsharp));
    enablesharpConn = expsharp->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expsharp));

    sharcontrast->setAdjusterListener(this);

    sharradius->setAdjusterListener(this);

    sharamount->setAdjusterListener(this);

    shardamping->setAdjusterListener(this);

    shariter->setAdjusterListener(this);

    sharblur->setAdjusterListener(this);

    sensisha->set_tooltip_text(M("TP_LOCALLAB_SENSIS_TOOLTIP"));
    sensisha->setAdjusterListener(this);

    inversshaConn  = inverssha->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversshaChanged));

    ToolParamBlock* const sharpBox = Gtk::manage(new ToolParamBlock());
    sharpBox->pack_start(*sharcontrast);
    sharpBox->pack_start(*sharradius);
    sharpBox->pack_start(*sharamount);
    sharpBox->pack_start(*shardamping);
    sharpBox->pack_start(*shariter);
    sharpBox->pack_start(*sharblur);
    sharpBox->pack_start(*sensisha);
    sharpBox->pack_start(*inverssha);
    expsharp->add(*sharpBox);
    expsharp->setLevel(2);

    panel->pack_start(*expsharp, false, false);

    // Local Contrast
    expcontrast->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcontrast));
    enablecontrastConn = expcontrast->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcontrast));

    lcradius->setAdjusterListener(this);

    lcamount->setAdjusterListener(this);

    lcdarkness->setAdjusterListener(this);

    lclightness->setAdjusterListener(this);

    sensilc->setAdjusterListener(this);

    ToolParamBlock* const contrastBox = Gtk::manage(new ToolParamBlock());
    contrastBox->pack_start(*lcradius);
    contrastBox->pack_start(*lcamount);
    contrastBox->pack_start(*lcdarkness);
    contrastBox->pack_start(*lclightness);
    contrastBox->pack_start(*sensilc);
    expcontrast->add(*contrastBox);
    expcontrast->setLevel(2);

    panel->pack_start(*expcontrast, false, false);

    // Contrast by detail levels
    expcbdl->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcbdl));
    enablecbdlConn = expcbdl->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcbdl));
    expcbdl->set_tooltip_text(M("TP_LOCALLAB_EXPCBDL_TOOLTIP"));

    for (int i = 0; i < 5; i++) {
        Glib::ustring ss;
        ss = Glib::ustring::format(i);

        if (i == 0) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMAFINEST"));
        } else if (i == 4) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        }

        multiplier[i] = Gtk::manage(new Adjuster(ss, 0, 400, 1, 100));
        multiplier[i]->setAdjusterListener(this);
    }

    chromacbdl->set_tooltip_text(M("TP_LOCALLAB_CHROMACB_TOOLTIP"));
    chromacbdl->setAdjusterListener(this);

    threshold->setAdjusterListener(this);

    sensicb->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensicb->setAdjusterListener(this);

    ToolParamBlock* const cbdlBox = Gtk::manage(new ToolParamBlock());
    Gtk::HBox* buttonBox = Gtk::manage(new Gtk::HBox(true, 10));
    buttonBox->pack_start(*lumacontrastMinusButton);
    lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumacontrastMinusPressed));
    buttonBox->pack_start(*lumaneutralButton);
    lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumaneutralPressed));
    buttonBox->pack_start(*lumacontrastPlusButton);
    lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumacontrastPlusPressed));
    cbdlBox->pack_start(*buttonBox);

    for (int i = 0; i < 5; i++) {
        cbdlBox->pack_start(*multiplier[i]);
    }

    Gtk::HSeparator *separator = Gtk::manage(new  Gtk::HSeparator());
    cbdlBox->pack_start(*separator, Gtk::PACK_SHRINK, 2);
    cbdlBox->pack_start(*chromacbdl);
    cbdlBox->pack_start(*threshold);
    cbdlBox->pack_start(*sensicb);
    expcbdl->add(*cbdlBox);
    expcbdl->setLevel(2);

    panel->pack_start(*expcbdl, false, false);

    // Denoise
    expdenoi->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expdenoi));
    enabledenoiConn = expdenoi->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expdenoi));

    noiselumf->setAdjusterListener(this);

    noiselumc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    noiselumc->setAdjusterListener(this);

    noiselumdetail->setAdjusterListener(this);

    noiselequal->setAdjusterListener(this);

    noisechrof->setAdjusterListener(this);

    noisechroc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    noisechroc->setAdjusterListener(this);

    noisechrodetail->setAdjusterListener(this);

    adjblur->setAdjusterListener(this);

    bilateral->setAdjusterListener(this);

    sensiden->setAdjusterListener(this);

    avoidConn  = avoid->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::avoidChanged));

    ToolParamBlock* const denoisBox = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const wavFrame = Gtk::manage(new Gtk::Frame());
    ToolParamBlock* const wavBox = Gtk::manage(new ToolParamBlock());
    wavBox->pack_start(*noiselumf);
    wavBox->pack_start(*noiselumc);
    wavBox->pack_start(*noiselumdetail);
    wavBox->pack_start(*noiselequal);
    wavBox->pack_start(*noisechrof);
    wavBox->pack_start(*noisechroc);
    // wavBox->pack_start(*noisechrodetail); // Uncomment this line to use the noisechrodetail adjuster
    wavBox->pack_start(*adjblur);
    wavFrame->add(*wavBox);
    denoisBox->pack_start(*wavFrame);
    denoisBox->pack_start(*bilateral);
    denoisBox->pack_start(*sensiden);
    expdenoi->add(*denoisBox);
    expdenoi->setLevel(2);

    panel->pack_start(*expdenoi, false, false);

    panel->pack_start(*avoid);

    pack_start(*panel);

    setParamEditable(false);

    show_all();
}

Locallab::~Locallab()
{
    delete llCurveEditorG;
    delete curveEditorG;
    delete curveEditorGG;
    delete LocalcurveEditorgainT;
}
void Locallab::foldAllButMe(GdkEventButton* event, MyExpander *expander)
{
    if (event->button == 3) {
        expsettings->setExpanded(expsettings->getExpander() == expander);
        expcolor->set_expanded(expcolor == expander);
        expexpose->set_expanded(expexpose == expander);
        expvibrance->set_expanded(expvibrance == expander);
        expsoft->set_expanded(expsoft == expander);
        explabregion->set_expanded(explabregion == expander);
        expblur->set_expanded(expblur == expander);
        exptonemap->set_expanded(exptonemap == expander);
        expreti->set_expanded(expreti == expander);
        expsharp->set_expanded(expsharp == expander);
        expcontrast->set_expanded(expcontrast == expander);
        expcbdl->set_expanded(expcbdl == expander);
        expdenoi->set_expanded(expdenoi == expander);

    }
}

void Locallab::enableToggled(MyExpander *expander)
{
    // printf("enableToggled\n");

    rtengine::ProcEvent event;
    sigc::connection* expConn;

    if (expander == expcolor) {
        event = EvLocenacolor;
        expConn = &enablecolorConn;
    } else if (expander == expexpose) {
        event = EvLocenaexpose;
        expConn = &enableexposeConn;
    } else if (expander == expvibrance) {
        event = EvLocenavibrance;
        expConn = &enablevibranceConn;
    } else if (expander == expsoft) {
        event = EvLocenasoft;
        expConn = &enablesoftConn;
    } else if (expander == explabregion) {
        event = EvLocenalabregion;
        expConn = &enablelabregionConn;
    } else if (expander == expblur) {
        event = EvLocenablur;
        expConn = &enableblurConn;
    } else if (expander == exptonemap) {
        event = EvLocenatonemap;
        expConn = &enabletonemapConn;
    } else if (expander == expreti) {
        event = EvLocenareti;
        expConn = &enableretiConn;
    } else if (expander == expsharp) {
        event = EvLocenasharp;
        expConn = &enablesharpConn;
    } else if (expander == expcontrast) {
        event = EvLocenacontrast;
        expConn = &enablecontrastConn;
    } else if (expander == expcbdl) {
        event = EvLocenacbdl;
        expConn = &enablecbdlConn;
    } else if (expander == expdenoi) {
        event = EvLocenadenoi;
        expConn = &enabledenoiConn;
    } else {
        return;
    }

    if (multiImage) {
        if (expander->get_inconsistent()) {
            expander->set_inconsistent(false);
            expConn->block(true);
            expander->setEnabled(false);
            expConn->block(false);
        }
    }

    if (getEnabled()) {
        if (listener) {
            if (expander->getEnabled()) {
                listener->panelChanged(event, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(event, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::writeOptions(std::vector<int> &tpOpen)
{
    tpOpen.push_back(expsettings->getExpanded());
    tpOpen.push_back(expcolor->get_expanded());
    tpOpen.push_back(expexpose->get_expanded());
    tpOpen.push_back(expvibrance->get_expanded());
    tpOpen.push_back(expsoft->get_expanded());
    tpOpen.push_back(explabregion->get_expanded());
    tpOpen.push_back(expblur->get_expanded());
    tpOpen.push_back(exptonemap->get_expanded());
    tpOpen.push_back(expreti->get_expanded());
    tpOpen.push_back(expsharp->get_expanded());
    tpOpen.push_back(expcontrast->get_expanded());
    tpOpen.push_back(expcbdl->get_expanded());
    tpOpen.push_back(expdenoi->get_expanded());

}

void Locallab::refChanged (double huer, double lumar, double chromar)
{
    nexthuer = huer;
    nextlumar = lumar / 100.f;
    nextchromar = chromar / 137.4f;
    //printf("nh=%f nl=%f nc=%f\n", nexthuer, nextlumar, nextchromar);
    const auto func = [] (gpointer data) -> gboolean {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected
        static_cast<Locallab*> (data)->refComputed_();

        return FALSE;
    };

    idle_register.add (func, this);
}

bool Locallab::refComputed_ ()
{

    disableListener ();
    enableListener ();
    updateLabel ();
    return false;
}

void Locallab::updateLabel ()
{
    if (!batchMode) {
        float nX, nY;
        
        nY = nextlumar;
        nX = nextchromar;
        {
            transLabels->set_text (
                Glib::ustring::compose (M ("TP_LOCALLAB_REFLABEL"),
                                        Glib::ustring::format (std::fixed, std::setprecision (3), nX),
                                        Glib::ustring::format (std::fixed, std::setprecision (3), nY))
            );
            transLabels2->set_text (
                Glib::ustring::compose (M ("TP_LOCALLAB_REFLABEL"),
                                        Glib::ustring::format (std::fixed, std::setprecision (3), nX),
                                        Glib::ustring::format (std::fixed, std::setprecision (3), nY))
            );
        }
    }
}

void Locallab::updateToolState(std::vector<int> &tpOpen)
{
    if (tpOpen.size() >= 13) {
        expsettings->setExpanded(tpOpen.at(0));
        expcolor->set_expanded(tpOpen.at(1));
        expexpose->set_expanded(tpOpen.at(2));
        expvibrance->set_expanded(tpOpen.at(3));
        expsoft->set_expanded(tpOpen.at(4));
        explabregion->set_expanded(tpOpen.at(5));
        expblur->set_expanded(tpOpen.at(6));
        exptonemap->set_expanded(tpOpen.at(7));
        expreti->set_expanded(tpOpen.at(8));
        expsharp->set_expanded(tpOpen.at(9));
        expcontrast->set_expanded(tpOpen.at(10));
        expcbdl->set_expanded(tpOpen.at(11));
        expdenoi->set_expanded(tpOpen.at(12));
    }
}

void Locallab::lumaneutralPressed()
{
    // printf("lumaneutralPressed\n");

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue(100);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void Locallab::lumacontrastPlusPressed()
{
    // printf("lumacontrastPlusPressed\n");

    for (int i = 0; i < 5; i++) {
        float inc = (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void Locallab::lumacontrastMinusPressed()
{
    // printf("lumacontrastMinusPressed\n");

    for (int i = 0; i < 5; i++) {
        float inc = - (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void Locallab::read(const ProcParams* pp, const ParamsEdited* pedited)
{
    // printf("Locallab read\n");

    // Disable all listeners
    disableListener();

    setEnabled(pp->locallab.enabled);

    if (pedited) {
        set_inconsistent(multiImage && !pedited->locallab.enabled);
    }

    // Delete all existent spots
    std::vector<int>* const list = expsettings->getSpotIdList();

    for (size_t i = 0; i < list->size(); i++) {
        expsettings->deleteControlSpot(list->at(i));
    }

    // Add existent spots based on pp
    ControlSpotPanel::SpotRow* const r = new ControlSpotPanel::SpotRow();

    for (int i = 0; i < pp->locallab.nbspot && i < (int)pp->locallab.spots.size(); i++) {
        r->id = pp->locallab.spots.at(i).id;
        r->name = pp->locallab.spots.at(i).name;
        r->isvisible = pp->locallab.spots.at(i).isvisible;

        if (pp->locallab.spots.at(i).shape == "ELI") {
            r->shape = 0;
        } else {
            r->shape = 1;
        }

        if (pp->locallab.spots.at(i).spotMethod == "norm") {
            r->spotMethod = 0;
        } else {
            r->spotMethod = 1;
        }

        r->sensiexclu = pp->locallab.spots.at(i).sensiexclu;
        r->struc = pp->locallab.spots.at(i).struc;

        if (pp->locallab.spots.at(i).shapeMethod == "IND") {
            r->shapeMethod = 0;
        } else if (pp->locallab.spots.at(i).shapeMethod == "SYM") {
            r->shapeMethod = 1;
        } else if (pp->locallab.spots.at(i).shapeMethod == "INDSL") {
            r->shapeMethod = 2;
        } else {
            r->shapeMethod = 3;
        }

        r->locX = pp->locallab.spots.at(i).locX;
        r->locXL = pp->locallab.spots.at(i).locXL;
        r->locY = pp->locallab.spots.at(i).locY;
        r->locYT = pp->locallab.spots.at(i).locYT;
        r->centerX = pp->locallab.spots.at(i).centerX;
        r->centerY = pp->locallab.spots.at(i).centerY;
        r->circrad = pp->locallab.spots.at(i).circrad;

        if (pp->locallab.spots.at(i).qualityMethod == "std") {
            r->qualityMethod = 0;
        } else if (pp->locallab.spots.at(i).qualityMethod == "enh") {
            r->qualityMethod = 1;
        } else {
            r->qualityMethod = 2;
        }

        r->transit = pp->locallab.spots.at(i).transit;
        r->thresh = pp->locallab.spots.at(i).thresh;
        r->iter = pp->locallab.spots.at(i).iter;

        expsettings->addControlSpot(r);
    }

    // Select active spot
    if (pp->locallab.nbspot > 0) {
        expsettings->setSelectedSpot(pp->locallab.spots.at(pp->locallab.selspot).id);
    }

    // Update Locallab tools GUI
    updateLocallabGUI(pp, pedited, pp->locallab.selspot);
    updateSpecificGUIState();

    if (pp->locallab.nbspot > 0) {
        setParamEditable(true);

        // Locallab params are not editable if nbspot, selspot or id are not coherent (batch mode)
        if (pedited) {
            if (!pedited->locallab.nbspot || !pedited->locallab.selspot || !pedited->locallab.id) {
                setParamEditable(false);
            }
        }
    } else {
        setParamEditable(false);
    }

    // Enable all listeners
    enableListener();

    // Update default values according to selected spot
    if (pp->locallab.nbspot > 0 && pp->locallab.selspot < (int)pp->locallab.spots.size()) {
        setDefaults(defparams, defpedited, pp->locallab.spots.at(pp->locallab.selspot).id);
    }
}

void Locallab::write(ProcParams* pp, ParamsEdited* pedited)
{
    // printf("Locallab write\n");

    pp->locallab.enabled = getEnabled();

    const int spotPanelEvent = expsettings->getEventType();
    int spotId;
    ControlSpotPanel::SpotRow* r;
    LocallabParams::LocallabSpot* newSpot;

    switch (spotPanelEvent) {
        case (1): // 1 = Spot creation event
            // Spot creation (default initialization)
            newSpot = new LocallabParams::LocallabSpot();
            spotId = expsettings->getNewId();
            r = new ControlSpotPanel::SpotRow();
            r->id = newSpot->id = spotId;
            r->name = newSpot->name = M("TP_LOCALLAB_SPOTNAME") + std::to_string(spotId);
            r->isvisible = newSpot->isvisible;

            if (newSpot->shape == "ELI") {
                r->shape = 0;
            } else {
                r->shape = 1;
            }

            if (newSpot->spotMethod == "norm") {
                r->spotMethod = 0;
            } else {
                r->spotMethod = 1;
            }

            r->sensiexclu = newSpot->sensiexclu;
            r->struc = newSpot->struc;

            if (newSpot->shapeMethod == "IND") {
                r->shapeMethod = 0;
            } else if (newSpot->shapeMethod == "SYM") {
                r->shapeMethod = 1;
            } else if (newSpot->shapeMethod == "INDSL") {
                r->shapeMethod = 2;
            } else {
                r->shapeMethod = 3;
            }

            r->locX = newSpot->locX;
            r->locXL = newSpot->locXL;
            r->locY = newSpot->locY;
            r->locYT = newSpot->locYT;
            r->centerX = newSpot->centerX;
            r->centerY = newSpot->centerY;
            r->circrad = newSpot->circrad;

            if (newSpot->qualityMethod == "std") {
                r->qualityMethod = 0;
            } else if (newSpot->qualityMethod == "enh") {
                r->qualityMethod = 1;
            } else {
                r->qualityMethod = 2;
            }

            r->transit = newSpot->transit;
            r->thresh = newSpot->thresh;
            r->iter = newSpot->iter;
            expsettings->addControlSpot(r);

            // ProcParams update
            pp->locallab.nbspot++;
            pp->locallab.selspot = pp->locallab.nbspot - 1;
            pp->locallab.spots.push_back(*newSpot);

            // New created spot selection
            expsettings->setSelectedSpot(spotId);

            // Update Locallab tools GUI with new created spot
            disableListener();

            if (pe) {
                pe->locallab.spots.push_back(new LocallabParamsEdited::LocallabSpotEdited(true));
            }

            updateLocallabGUI(pp, pe, pp->locallab.selspot);

            enableListener();

            if (pp->locallab.nbspot == 1) {
                setParamEditable(true);
            }

            // Update default values according to selected spot
            setDefaults(defparams, defpedited, spotId);

            // ParamsEdited update
            if (pedited) {
                pedited->locallab.nbspot = true;
                pedited->locallab.selspot = true;
                pedited->locallab.id = true;
                pedited->locallab.spots.push_back(new LocallabParamsEdited::LocallabSpotEdited(true));
            }

            break;

        case (2): // 2 = Spot deletion event
            // Get deleted spot index in ProcParams and update it
            spotId = expsettings->getSelectedSpot();

            for (int i = 0; i < pp->locallab.nbspot && i < (int)pp->locallab.spots.size(); i++) {
                if (pp->locallab.spots.at(i).id == spotId) {
                    // ProcParams update
                    pp->locallab.nbspot--;
                    pp->locallab.selspot = 0;
                    pp->locallab.spots.erase(pp->locallab.spots.begin() + i);
                    expsettings->deleteControlSpot(spotId);

                    // Select one remaining spot
                    if (pp->locallab.nbspot > 0) {
                        expsettings->setSelectedSpot(pp->locallab.spots.at(pp->locallab.selspot).id);
                    }

                    // Update Locallab tools GUI with new created spot
                    disableListener();

                    if (pe) {
                        if (i < (int)pe->locallab.spots.size()) {
                            pe->locallab.spots.erase(pe->locallab.spots.begin() + i);
                        }
                    }

                    updateLocallabGUI(pp, pe, pp->locallab.selspot);

                    enableListener();

                    if (pp->locallab.nbspot == 0) {
                        setParamEditable(false);
                    }

                    // Update default values according to selected spot
                    if (pp->locallab.nbspot > 0) {
                        setDefaults(defparams, defpedited, pp->locallab.spots.at(pp->locallab.selspot).id);
                    }

                    // ParamsEdited update
                    if (pedited) {
                        pedited->locallab.nbspot = true;
                        pedited->locallab.selspot = true;
                        pedited->locallab.id = true;

                        if (i < (int)pedited->locallab.spots.size()) {
                            pedited->locallab.spots.erase(pedited->locallab.spots.begin() + i);
                        }
                    }

                    break;
                }
            }

            break;

        case (3):  // 3 = Spot selection event
            spotId = expsettings->getSelectedSpot();

            for (int i = 0; i < pp->locallab.nbspot && i < (int)pp->locallab.spots.size(); i++) {
                if (pp->locallab.spots.at(i).id == spotId) {
                    pp->locallab.selspot = i;
                    break;
                }
            }

            // Update control spots and Locallab tools GUI with selected spot
            expsettings->setSelectedSpot(spotId);
            disableListener();
            updateLocallabGUI(pp, pe, pp->locallab.selspot);
            enableListener();

            // Update default values according to selected spot
            setDefaults(defparams, defpedited, spotId);

            // ParamsEdited update
            if (pedited) {
                pedited->locallab.selspot = true;
            }

            break;

        case (4): // 4 = Spot duplication event
            newSpot = nullptr;
            spotId = expsettings->getSelectedSpot();

            for (int i = 0; i < pp->locallab.nbspot && i < (int)pp->locallab.spots.size(); i++) {
                if (pp->locallab.spots.at(i).id == spotId) {
                    newSpot = new LocallabParams::LocallabSpot(pp->locallab.spots.at(i));
                    break;
                }
            }

            if (!newSpot) {
                break;
            }

            // Spot creation (initialization at currently selected spot)
            spotId = expsettings->getNewId();
            r = new ControlSpotPanel::SpotRow();
            r->id = newSpot->id = spotId;
            r->name = newSpot->name = newSpot->name + " - " + M("TP_LOCALLAB_DUPLSPOTNAME");
            r->isvisible = newSpot->isvisible;

            if (newSpot->shape == "ELI") {
                r->shape = 0;
            } else {
                r->shape = 1;
            }

            if (newSpot->spotMethod == "norm") {
                r->spotMethod = 0;
            } else {
                r->spotMethod = 1;
            }

            r->sensiexclu = newSpot->sensiexclu;
            r->struc = newSpot->struc;

            if (newSpot->shapeMethod == "IND") {
                r->shapeMethod = 0;
            } else if (newSpot->shapeMethod == "SYM") {
                r->shapeMethod = 1;
            } else if (newSpot->shapeMethod == "INDSL") {
                r->shapeMethod = 2;
            } else {
                r->shapeMethod = 3;
            }

            r->locX = newSpot->locX;
            r->locXL = newSpot->locXL;
            r->locY = newSpot->locY;
            r->locYT = newSpot->locYT;
            r->centerX = newSpot->centerX;
            r->centerY = newSpot->centerY;
            r->circrad = newSpot->circrad;

            if (newSpot->qualityMethod == "std") {
                r->qualityMethod = 0;
            } else if (newSpot->qualityMethod == "enh") {
                r->qualityMethod = 1;
            } else {
                r->qualityMethod = 2;
            }

            r->transit = newSpot->transit;
            r->thresh = newSpot->thresh;
            r->iter = newSpot->iter;
            expsettings->addControlSpot(r);

            // ProcParams update
            pp->locallab.nbspot++;
            pp->locallab.selspot = pp->locallab.nbspot - 1;
            pp->locallab.spots.push_back(*newSpot);

            // New created spot selection
            expsettings->setSelectedSpot(spotId);

            // Update Locallab tools GUI with new created spot
            disableListener();

            if (pe) {
                pe->locallab.spots.push_back(new LocallabParamsEdited::LocallabSpotEdited(true));
            }

            updateLocallabGUI(pp, pe, pp->locallab.selspot);

            enableListener();

            // Update default values according to selected spot
            setDefaults(defparams, defpedited, spotId);

            // ParamsEdited update
            if (pedited) {
                pedited->locallab.nbspot = true;
                pedited->locallab.selspot = true;
                pedited->locallab.id = true;
                pedited->locallab.spots.push_back(new LocallabParamsEdited::LocallabSpotEdited(true));
            }

            break;

        default: // Spot or locallab GUI updated
            if (pp->locallab.nbspot > 0) {
                r = expsettings->getSpot(expsettings->getSelectedSpot());

                // ProcParams update
                if (pp->locallab.selspot < (int)pp->locallab.spots.size()) {
                    // Control spot settings
                    pp->locallab.spots.at(pp->locallab.selspot).name = r->name;
                    pp->locallab.spots.at(pp->locallab.selspot).isvisible = r->isvisible;

                    if (r->shape == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).shape = "ELI";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).shape = "RECT";
                    }

                    if (r->spotMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).spotMethod = "norm";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).spotMethod = "exc";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).sensiexclu = r->sensiexclu;
                    pp->locallab.spots.at(pp->locallab.selspot).struc = r->struc;

                    if (r->shapeMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "IND";
                    } else if (r->shapeMethod == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "SYM";
                    } else if (r->shapeMethod == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "INDSL";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "SYMSL";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).locX = r->locX;
                    pp->locallab.spots.at(pp->locallab.selspot).locXL = r->locXL;
                    pp->locallab.spots.at(pp->locallab.selspot).locY = r->locY;
                    pp->locallab.spots.at(pp->locallab.selspot).locYT = r->locYT;
                    pp->locallab.spots.at(pp->locallab.selspot).centerX = r->centerX;
                    pp->locallab.spots.at(pp->locallab.selspot).centerY = r->centerY;
                    pp->locallab.spots.at(pp->locallab.selspot).circrad = r->circrad;

                    if (r->qualityMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "std";
                    } else if (r->qualityMethod == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "enh";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "enhden";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).transit = r->transit;
                    pp->locallab.spots.at(pp->locallab.selspot).thresh = r->thresh;
                    pp->locallab.spots.at(pp->locallab.selspot).iter = r->iter;
                    // Color & Light
                    pp->locallab.spots.at(pp->locallab.selspot).expcolor = expcolor->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).curvactiv = curvactiv->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).lightness = lightness->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).contrast = contrast->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).chroma = chroma->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensi = sensi->getIntValue();

                    if (qualitycurveMethod->get_active_row_number() == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = "none";
                    } else if (qualitycurveMethod->get_active_row_number() == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = "std";
                    } else if (qualitycurveMethod->get_active_row_number() == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = "enh";
                    }

                    if (showmaskcolMethod->get_active_row_number() == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod = "none";
                    } else if (showmaskcolMethod->get_active_row_number() == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod = "color";
                    } else if (showmaskcolMethod->get_active_row_number() == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod = "colormask";
                    } else if (showmaskcolMethod->get_active_row_number() == 3) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod = "mask";
                    } else if (showmaskcolMethod->get_active_row_number() == 4) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod = "showmask";
                    }
                 
                    pp->locallab.spots.at(pp->locallab.selspot).llcurve = llshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).cccurve = ccshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).LHcurve = LHshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).HHcurve = HHshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).CCmaskcurve = CCmaskshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).LLmaskcurve = LLmaskshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).invers = invers->get_active();
                    // Exposure
                    pp->locallab.spots.at(pp->locallab.selspot).expexpose = expexpose->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).expcomp = expcomp->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).hlcompr = hlcompr->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).hlcomprthresh = hlcomprthresh->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).black = black->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shcompr = shcompr->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).warm = warm->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensiex = sensiex->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).excurve = shapeexpos->getCurve();
                    if (showmaskexpMethod->get_active_row_number() == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod = "none";
                    } else if (showmaskexpMethod->get_active_row_number() == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod = "expo";
                    } else if (showmaskexpMethod->get_active_row_number() == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod = "expomask";
                    } else if (showmaskexpMethod->get_active_row_number() == 3) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod = "mask";
                    } else if (showmaskexpMethod->get_active_row_number() == 4) {
                        pp->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod = "showmask";
                    }
                    pp->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve = LLmaskexpshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve = CCmaskexpshape->getCurve();

                    // Vibrance
                    pp->locallab.spots.at(pp->locallab.selspot).expvibrance = expvibrance->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).saturated = saturated->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).pastels = pastels->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).psthreshold = psThreshold->getValue<int>();
                    pp->locallab.spots.at(pp->locallab.selspot).protectskins = protectSkins->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).avoidcolorshift = avoidColorShift->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).pastsattog = pastSatTog->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).sensiv = sensiv->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).skintonescurve = skinTonesCurve->getCurve();
                    // Soft Light
                    pp->locallab.spots.at(pp->locallab.selspot).expsoft = expsoft->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).streng = streng->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensisf = sensisf->getIntValue();
                    // Lab Region
                    pp->locallab.spots.at(pp->locallab.selspot).explabregion = explabregion->getEnabled();
                    // Blur & Noise
                    pp->locallab.spots.at(pp->locallab.selspot).expblur = expblur->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).radius = radius->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).strength = strength->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensibn = sensibn->getIntValue();

                    if (blurMethod->get_active_row_number() == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).blurMethod = "norm";
                    } else if (blurMethod->get_active_row_number() == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).blurMethod = "inv";
                    } else if (blurMethod->get_active_row_number() == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).blurMethod = "sym";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).activlum = activlum->get_active();
                    // Tone Mapping
                    pp->locallab.spots.at(pp->locallab.selspot).exptonemap = exptonemap->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).stren = stren->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).gamma = gamma->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).estop = estop->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).scaltm = scaltm->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).rewei = rewei->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensitm = sensitm->getIntValue();
                    // Retinex
                    pp->locallab.spots.at(pp->locallab.selspot).expreti = expreti->getEnabled();

                    if (retinexMethod->get_active_row_number() == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).retinexMethod = "low";
                    } else if (retinexMethod->get_active_row_number() == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).retinexMethod = "uni";
                    } else if (retinexMethod->get_active_row_number() == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).retinexMethod = "high";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).str = str->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).chrrt = chrrt->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).neigh = neigh->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).vart = vart->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).dehaz = dehaz->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensih = sensih->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).localTgaincurve = cTgainshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).inversret = inversret->get_active();
                    // Sharpening
                    pp->locallab.spots.at(pp->locallab.selspot).expsharp = expsharp->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).sharcontrast = sharcontrast->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sharradius = sharradius->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sharamount = sharamount->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shardamping = shardamping->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shariter = shariter->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sharblur = sharblur->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensisha = sensisha->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).inverssha = inverssha->get_active();
                    // Local Contrast
                    pp->locallab.spots.at(pp->locallab.selspot).expcontrast = expcontrast->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).lcradius = lcradius->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).lcamount = lcamount->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).lcdarkness = lcdarkness->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).lclightness = lclightness->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensilc = sensilc->getIntValue();
                    // Contrast by detail levels
                    pp->locallab.spots.at(pp->locallab.selspot).expcbdl = expcbdl->getEnabled();

                    for (int i = 0; i < 5; i++) {
                        pp->locallab.spots.at(pp->locallab.selspot).mult[i] = multiplier[i]->getIntValue();
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).chromacbdl = chromacbdl->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).threshold = threshold->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensicb = sensicb->getIntValue();
                    // Denoise
                    pp->locallab.spots.at(pp->locallab.selspot).expdenoi = expdenoi->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselumf = noiselumf->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselumc = noiselumc->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselumdetail = noiselumdetail->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselequal = noiselequal->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noisechrof = noisechrof->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noisechroc = noisechroc->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noisechrodetail = noisechrodetail->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).adjblur = adjblur->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).bilateral = bilateral->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensiden = sensiden->getIntValue();
                    // Others
                    pp->locallab.spots.at(pp->locallab.selspot).avoid = avoid->get_active();
                }

                ControlSpotPanel::SpotEdited* const se = expsettings->getEditedStates();

                if (pe) {
                    if (pp->locallab.selspot < (int)pe->locallab.spots.size()) {
                        pe->locallab.spots.at(pp->locallab.selspot).name = pe->locallab.spots.at(pp->locallab.selspot).name || se->name;
                        pe->locallab.spots.at(pp->locallab.selspot).isvisible = pe->locallab.spots.at(pp->locallab.selspot).isvisible || se->isvisible;
                        pe->locallab.spots.at(pp->locallab.selspot).shape = pe->locallab.spots.at(pp->locallab.selspot).shape || se->shape;
                        pe->locallab.spots.at(pp->locallab.selspot).spotMethod = pe->locallab.spots.at(pp->locallab.selspot).spotMethod || se->spotMethod;
                        pe->locallab.spots.at(pp->locallab.selspot).sensiexclu = pe->locallab.spots.at(pp->locallab.selspot).sensiexclu || se->sensiexclu;
                        pe->locallab.spots.at(pp->locallab.selspot).struc = pe->locallab.spots.at(pp->locallab.selspot).struc || se->struc;
                        pe->locallab.spots.at(pp->locallab.selspot).shapeMethod = pe->locallab.spots.at(pp->locallab.selspot).shapeMethod || se->shapeMethod;
                        pe->locallab.spots.at(pp->locallab.selspot).locX = pe->locallab.spots.at(pp->locallab.selspot).locX || se->locX;
                        pe->locallab.spots.at(pp->locallab.selspot).locXL = pe->locallab.spots.at(pp->locallab.selspot).locXL || se->locXL;
                        pe->locallab.spots.at(pp->locallab.selspot).locY = pe->locallab.spots.at(pp->locallab.selspot).locY || se->locY;
                        pe->locallab.spots.at(pp->locallab.selspot).locYT = pe->locallab.spots.at(pp->locallab.selspot).locYT || se->locYT;
                        pe->locallab.spots.at(pp->locallab.selspot).centerX = pe->locallab.spots.at(pp->locallab.selspot).centerX || se->centerX;
                        pe->locallab.spots.at(pp->locallab.selspot).centerY = pe->locallab.spots.at(pp->locallab.selspot).centerY || se->centerY;
                        pe->locallab.spots.at(pp->locallab.selspot).circrad = pe->locallab.spots.at(pp->locallab.selspot).circrad || se->circrad;
                        pe->locallab.spots.at(pp->locallab.selspot).qualityMethod = pe->locallab.spots.at(pp->locallab.selspot).qualityMethod || se->qualityMethod;
                        pe->locallab.spots.at(pp->locallab.selspot).transit = pe->locallab.spots.at(pp->locallab.selspot).transit || se->transit;
                        pe->locallab.spots.at(pp->locallab.selspot).thresh = pe->locallab.spots.at(pp->locallab.selspot).thresh || se->thresh;
                        pe->locallab.spots.at(pp->locallab.selspot).iter = pe->locallab.spots.at(pp->locallab.selspot).iter || se->iter;
                        // Color & Light
                        pe->locallab.spots.at(pp->locallab.selspot).expcolor = pe->locallab.spots.at(pp->locallab.selspot).expcolor || !expcolor->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).curvactiv = pe->locallab.spots.at(pp->locallab.selspot).curvactiv || !curvactiv->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).lightness = pe->locallab.spots.at(pp->locallab.selspot).lightness || lightness->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).contrast = pe->locallab.spots.at(pp->locallab.selspot).contrast || contrast->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).chroma = pe->locallab.spots.at(pp->locallab.selspot).chroma || chroma->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensi = pe->locallab.spots.at(pp->locallab.selspot).sensi || sensi->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = pe->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod || qualitycurveMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pe->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod = pe->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod || showmaskcolMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pe->locallab.spots.at(pp->locallab.selspot).llcurve = pe->locallab.spots.at(pp->locallab.selspot).llcurve || !llshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).cccurve = pe->locallab.spots.at(pp->locallab.selspot).cccurve || !ccshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).LHcurve = pe->locallab.spots.at(pp->locallab.selspot).LHcurve || !LHshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).HHcurve = pe->locallab.spots.at(pp->locallab.selspot).HHcurve || !HHshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).CCmaskcurve = pe->locallab.spots.at(pp->locallab.selspot).CCmaskcurve || !CCmaskshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).LLmaskcurve = pe->locallab.spots.at(pp->locallab.selspot).LLmaskcurve || !LLmaskshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).invers = pe->locallab.spots.at(pp->locallab.selspot).invers || !invers->get_inconsistent();
                        // Exposure
                        pe->locallab.spots.at(pp->locallab.selspot).expexpose = pe->locallab.spots.at(pp->locallab.selspot).expexpose || !expexpose->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).expcomp = pe->locallab.spots.at(pp->locallab.selspot).expcomp || expcomp->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).hlcompr = pe->locallab.spots.at(pp->locallab.selspot).hlcompr || hlcompr->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).hlcomprthresh = pe->locallab.spots.at(pp->locallab.selspot).hlcomprthresh || hlcomprthresh->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).black = pe->locallab.spots.at(pp->locallab.selspot).black || black->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).shcompr = pe->locallab.spots.at(pp->locallab.selspot).shcompr || shcompr->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).warm = pe->locallab.spots.at(pp->locallab.selspot).warm || warm->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensiex = pe->locallab.spots.at(pp->locallab.selspot).sensiex || sensiex->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).excurve = pe->locallab.spots.at(pp->locallab.selspot).excurve || !shapeexpos->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod = pe->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod || showmaskexpMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pe->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve = pe->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve || !CCmaskexpshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve = pe->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve || !LLmaskexpshape->isUnChanged();
                        // Vibrance
                        pe->locallab.spots.at(pp->locallab.selspot).expvibrance = pe->locallab.spots.at(pp->locallab.selspot).expvibrance || !expvibrance->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).saturated = pe->locallab.spots.at(pp->locallab.selspot).saturated || saturated->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).pastels = pe->locallab.spots.at(pp->locallab.selspot).pastels || pastels->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).psthreshold = pe->locallab.spots.at(pp->locallab.selspot).psthreshold || psThreshold->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).protectskins = pe->locallab.spots.at(pp->locallab.selspot).protectskins || !protectSkins->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).avoidcolorshift = pe->locallab.spots.at(pp->locallab.selspot).avoidcolorshift || !avoidColorShift->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).pastsattog = pe->locallab.spots.at(pp->locallab.selspot).pastsattog || !pastSatTog->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).sensiv = pe->locallab.spots.at(pp->locallab.selspot).sensiv || sensiv->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).skintonescurve = pe->locallab.spots.at(pp->locallab.selspot).skintonescurve || !skinTonesCurve->isUnChanged();
                        // Soft Light
                        pe->locallab.spots.at(pp->locallab.selspot).expsoft = pe->locallab.spots.at(pp->locallab.selspot).expsoft || !expsoft->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).streng = pe->locallab.spots.at(pp->locallab.selspot).streng || streng->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensisf = pe->locallab.spots.at(pp->locallab.selspot).sensisf || sensisf->getEditedState();
                        // Lab Region
                        pe->locallab.spots.at(pp->locallab.selspot).explabregion = pe->locallab.spots.at(pp->locallab.selspot).explabregion || !explabregion->get_inconsistent();
                        // Blur & Noise
                        pe->locallab.spots.at(pp->locallab.selspot).expblur = pe->locallab.spots.at(pp->locallab.selspot).expblur || !expblur->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).radius = pe->locallab.spots.at(pp->locallab.selspot).radius || radius->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).strength = pe->locallab.spots.at(pp->locallab.selspot).strength || strength->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensibn = pe->locallab.spots.at(pp->locallab.selspot).sensibn || sensibn->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).blurMethod = pe->locallab.spots.at(pp->locallab.selspot).blurMethod || blurMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pe->locallab.spots.at(pp->locallab.selspot).activlum = pe->locallab.spots.at(pp->locallab.selspot).activlum || !activlum->get_inconsistent();
                        // Tone Mapping
                        pe->locallab.spots.at(pp->locallab.selspot).exptonemap = pe->locallab.spots.at(pp->locallab.selspot).activlum || !exptonemap->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).stren = pe->locallab.spots.at(pp->locallab.selspot).stren || stren->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).gamma = pe->locallab.spots.at(pp->locallab.selspot).gamma || gamma->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).estop = pe->locallab.spots.at(pp->locallab.selspot).estop || estop->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).scaltm = pe->locallab.spots.at(pp->locallab.selspot).scaltm || scaltm->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).rewei = pe->locallab.spots.at(pp->locallab.selspot).rewei || rewei->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensitm = pe->locallab.spots.at(pp->locallab.selspot).sensitm || sensitm->getEditedState();
                        // Retinex
                        pe->locallab.spots.at(pp->locallab.selspot).expreti = pe->locallab.spots.at(pp->locallab.selspot).expreti || !expreti->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).retinexMethod = pe->locallab.spots.at(pp->locallab.selspot).retinexMethod || retinexMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pe->locallab.spots.at(pp->locallab.selspot).str = pe->locallab.spots.at(pp->locallab.selspot).str || str->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).chrrt = pe->locallab.spots.at(pp->locallab.selspot).chrrt || chrrt->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).neigh = pe->locallab.spots.at(pp->locallab.selspot).neigh || neigh->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).vart = pe->locallab.spots.at(pp->locallab.selspot).vart || vart->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).dehaz = pe->locallab.spots.at(pp->locallab.selspot).dehaz || dehaz->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensih = pe->locallab.spots.at(pp->locallab.selspot).sensih || sensih->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).localTgaincurve = pe->locallab.spots.at(pp->locallab.selspot).localTgaincurve || !cTgainshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).inversret = pe->locallab.spots.at(pp->locallab.selspot).inversret || !inversret->get_inconsistent();
                        // Sharpening
                        pe->locallab.spots.at(pp->locallab.selspot).expsharp = pe->locallab.spots.at(pp->locallab.selspot).expsharp || !expsharp->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).sharcontrast = pe->locallab.spots.at(pp->locallab.selspot).sharcontrast || sharcontrast->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sharradius = pe->locallab.spots.at(pp->locallab.selspot).sharradius || sharradius->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sharamount = pe->locallab.spots.at(pp->locallab.selspot).sharamount || sharamount->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).shardamping = pe->locallab.spots.at(pp->locallab.selspot).shardamping || shardamping->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).shariter = pe->locallab.spots.at(pp->locallab.selspot).shariter || shariter->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sharblur = pe->locallab.spots.at(pp->locallab.selspot).sharblur || sharblur->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensisha = pe->locallab.spots.at(pp->locallab.selspot).sensisha || sensisha->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).inverssha = pe->locallab.spots.at(pp->locallab.selspot).inverssha || !inverssha->get_inconsistent();
                        // Local Contrast
                        pe->locallab.spots.at(pp->locallab.selspot).expcontrast = pe->locallab.spots.at(pp->locallab.selspot).expcontrast || !expcontrast->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).lcradius = pe->locallab.spots.at(pp->locallab.selspot).lcradius || lcradius->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).lcamount = pe->locallab.spots.at(pp->locallab.selspot).lcamount || lcamount->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).lcdarkness = pe->locallab.spots.at(pp->locallab.selspot).lcdarkness || lcdarkness->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).lclightness = pe->locallab.spots.at(pp->locallab.selspot).lclightness || lclightness->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensilc = pe->locallab.spots.at(pp->locallab.selspot).sensilc || sensilc->getEditedState();
                        // Contrast by detail levels
                        pe->locallab.spots.at(pp->locallab.selspot).expcbdl = pe->locallab.spots.at(pp->locallab.selspot).expcbdl || !expcbdl->get_inconsistent();

                        for (int i = 0; i < 5; i++) {
                            pe->locallab.spots.at(pp->locallab.selspot).mult[i] = pe->locallab.spots.at(pp->locallab.selspot).mult[i] || multiplier[i]->getEditedState();
                        }

                        pe->locallab.spots.at(pp->locallab.selspot).chromacbdl = pe->locallab.spots.at(pp->locallab.selspot).chromacbdl || chromacbdl->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).threshold = pe->locallab.spots.at(pp->locallab.selspot).threshold || threshold->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensicb = pe->locallab.spots.at(pp->locallab.selspot).sensicb || sensicb->getEditedState();
                        // Denoise
                        pe->locallab.spots.at(pp->locallab.selspot).expdenoi = pe->locallab.spots.at(pp->locallab.selspot).expdenoi || !expdenoi->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselumf = pe->locallab.spots.at(pp->locallab.selspot).noiselumf || noiselumf->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselumc = pe->locallab.spots.at(pp->locallab.selspot).noiselumc || noiselumc->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselumdetail = pe->locallab.spots.at(pp->locallab.selspot).noiselumdetail || noiselumdetail->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselequal = pe->locallab.spots.at(pp->locallab.selspot).noiselequal || noiselequal->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noisechrof = pe->locallab.spots.at(pp->locallab.selspot).noisechrof || noisechrof->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noisechroc = pe->locallab.spots.at(pp->locallab.selspot).noisechroc || noisechroc->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noisechrodetail = pe->locallab.spots.at(pp->locallab.selspot).noisechrodetail || noisechrodetail->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).adjblur = pe->locallab.spots.at(pp->locallab.selspot).adjblur || adjblur->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).bilateral = pe->locallab.spots.at(pp->locallab.selspot).bilateral || bilateral->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensiden = pe->locallab.spots.at(pp->locallab.selspot).sensiden || sensiden->getEditedState();
                        // Others
                        pe->locallab.spots.at(pp->locallab.selspot).avoid = pe->locallab.spots.at(pp->locallab.selspot).avoid || !avoid->get_inconsistent();
                    }
                }

                // ParamsEdited update
                if (pedited) {
                    pedited->locallab.enabled = pedited->locallab.enabled || !get_inconsistent();

                    if (pp->locallab.selspot < (int)pedited->locallab.spots.size()) {
                        // Control spot settings
                        pedited->locallab.spots.at(pp->locallab.selspot).name = pedited->locallab.spots.at(pp->locallab.selspot).name || se->name;
                        pedited->locallab.spots.at(pp->locallab.selspot).isvisible = pedited->locallab.spots.at(pp->locallab.selspot).isvisible || se->isvisible;
                        pedited->locallab.spots.at(pp->locallab.selspot).shape = pedited->locallab.spots.at(pp->locallab.selspot).shape || se->shape;
                        pedited->locallab.spots.at(pp->locallab.selspot).spotMethod = pedited->locallab.spots.at(pp->locallab.selspot).spotMethod || se->spotMethod;
                        pedited->locallab.spots.at(pp->locallab.selspot).sensiexclu = pedited->locallab.spots.at(pp->locallab.selspot).sensiexclu || se->sensiexclu;
                        pedited->locallab.spots.at(pp->locallab.selspot).struc = pedited->locallab.spots.at(pp->locallab.selspot).struc || se->struc;
                        pedited->locallab.spots.at(pp->locallab.selspot).shapeMethod = pedited->locallab.spots.at(pp->locallab.selspot).shapeMethod || se->shapeMethod;
                        pedited->locallab.spots.at(pp->locallab.selspot).locX = pedited->locallab.spots.at(pp->locallab.selspot).locX || se->locX;
                        pedited->locallab.spots.at(pp->locallab.selspot).locXL = pedited->locallab.spots.at(pp->locallab.selspot).locXL || se->locXL;
                        pedited->locallab.spots.at(pp->locallab.selspot).locY = pedited->locallab.spots.at(pp->locallab.selspot).locY || se->locY;
                        pedited->locallab.spots.at(pp->locallab.selspot).locYT = pedited->locallab.spots.at(pp->locallab.selspot).locYT || se->locYT;
                        pedited->locallab.spots.at(pp->locallab.selspot).centerX = pedited->locallab.spots.at(pp->locallab.selspot).centerX || se->centerX;
                        pedited->locallab.spots.at(pp->locallab.selspot).centerY = pedited->locallab.spots.at(pp->locallab.selspot).centerY || se->centerY;
                        pedited->locallab.spots.at(pp->locallab.selspot).circrad = pedited->locallab.spots.at(pp->locallab.selspot).circrad || se->circrad;
                        pedited->locallab.spots.at(pp->locallab.selspot).qualityMethod = pedited->locallab.spots.at(pp->locallab.selspot).qualityMethod || se->qualityMethod;
                        pedited->locallab.spots.at(pp->locallab.selspot).transit = pedited->locallab.spots.at(pp->locallab.selspot).transit || se->transit;
                        pedited->locallab.spots.at(pp->locallab.selspot).thresh = pedited->locallab.spots.at(pp->locallab.selspot).thresh || se->thresh;
                        pedited->locallab.spots.at(pp->locallab.selspot).iter = pedited->locallab.spots.at(pp->locallab.selspot).iter || se->iter;
                        // Color & Light
                        pedited->locallab.spots.at(pp->locallab.selspot).expcolor = pedited->locallab.spots.at(pp->locallab.selspot).expcolor || !expcolor->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).curvactiv = pedited->locallab.spots.at(pp->locallab.selspot).curvactiv || !curvactiv->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).lightness = pedited->locallab.spots.at(pp->locallab.selspot).lightness || lightness->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).contrast = pedited->locallab.spots.at(pp->locallab.selspot).contrast || contrast->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).chroma = pedited->locallab.spots.at(pp->locallab.selspot).chroma || chroma->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensi = pedited->locallab.spots.at(pp->locallab.selspot).sensi || sensi->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = pedited->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod || qualitycurveMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pedited->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod = pedited->locallab.spots.at(pp->locallab.selspot).showmaskcolMethod || showmaskcolMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pedited->locallab.spots.at(pp->locallab.selspot).llcurve = pedited->locallab.spots.at(pp->locallab.selspot).llcurve || !llshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).cccurve = pedited->locallab.spots.at(pp->locallab.selspot).cccurve || !ccshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).LHcurve = pedited->locallab.spots.at(pp->locallab.selspot).LHcurve || !LHshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).HHcurve = pedited->locallab.spots.at(pp->locallab.selspot).HHcurve || !HHshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).CCmaskcurve = pedited->locallab.spots.at(pp->locallab.selspot).CCmaskcurve || !CCmaskshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).LLmaskcurve = pedited->locallab.spots.at(pp->locallab.selspot).LLmaskcurve || !LLmaskshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).invers = pedited->locallab.spots.at(pp->locallab.selspot).invers || !invers->get_inconsistent();
                        // Exposure
                        pedited->locallab.spots.at(pp->locallab.selspot).expexpose = pedited->locallab.spots.at(pp->locallab.selspot).expexpose || !expexpose->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).expcomp = pedited->locallab.spots.at(pp->locallab.selspot).expcomp || expcomp->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).hlcompr = pedited->locallab.spots.at(pp->locallab.selspot).hlcompr || hlcompr->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).hlcomprthresh = pedited->locallab.spots.at(pp->locallab.selspot).hlcomprthresh || hlcomprthresh->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).black = pedited->locallab.spots.at(pp->locallab.selspot).black || black->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).shcompr = pedited->locallab.spots.at(pp->locallab.selspot).shcompr || shcompr->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).warm = pedited->locallab.spots.at(pp->locallab.selspot).warm || warm->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensiex = pedited->locallab.spots.at(pp->locallab.selspot).sensiex || sensiex->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).excurve = pedited->locallab.spots.at(pp->locallab.selspot).excurve || !shapeexpos->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod = pedited->locallab.spots.at(pp->locallab.selspot).showmaskexpMethod || showmaskexpMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pedited->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve = pedited->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve || !CCmaskexpshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve = pedited->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve || !LLmaskexpshape->isUnChanged();
                        // Vibrance
                        pedited->locallab.spots.at(pp->locallab.selspot).expvibrance = pedited->locallab.spots.at(pp->locallab.selspot).expvibrance || !expvibrance->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).saturated = pedited->locallab.spots.at(pp->locallab.selspot).saturated || saturated->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).pastels = pedited->locallab.spots.at(pp->locallab.selspot).pastels || pastels->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).psthreshold = pedited->locallab.spots.at(pp->locallab.selspot).psthreshold || psThreshold->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).protectskins = pedited->locallab.spots.at(pp->locallab.selspot).protectskins || !protectSkins->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).avoidcolorshift = pedited->locallab.spots.at(pp->locallab.selspot).avoidcolorshift || !avoidColorShift->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).pastsattog = pedited->locallab.spots.at(pp->locallab.selspot).pastsattog || !pastSatTog->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensiv = pedited->locallab.spots.at(pp->locallab.selspot).sensiv || sensiv->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).skintonescurve = pedited->locallab.spots.at(pp->locallab.selspot).skintonescurve || !skinTonesCurve->isUnChanged();
                        // Soft Light
                        pedited->locallab.spots.at(pp->locallab.selspot).expsoft = pedited->locallab.spots.at(pp->locallab.selspot).expsoft || !expsoft->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).streng = pedited->locallab.spots.at(pp->locallab.selspot).streng || streng->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensisf = pedited->locallab.spots.at(pp->locallab.selspot).sensisf || sensisf->getEditedState();
                        // Lab Region
                        pedited->locallab.spots.at(pp->locallab.selspot).explabregion = pedited->locallab.spots.at(pp->locallab.selspot).explabregion || !explabregion->get_inconsistent();
                        // Blur & Noise
                        pedited->locallab.spots.at(pp->locallab.selspot).expblur = pedited->locallab.spots.at(pp->locallab.selspot).expblur || !expblur->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).radius = pedited->locallab.spots.at(pp->locallab.selspot).radius || radius->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).strength = pedited->locallab.spots.at(pp->locallab.selspot).strength || strength->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensibn = pedited->locallab.spots.at(pp->locallab.selspot).sensibn || sensibn->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).blurMethod = pedited->locallab.spots.at(pp->locallab.selspot).blurMethod || blurMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pedited->locallab.spots.at(pp->locallab.selspot).activlum = pedited->locallab.spots.at(pp->locallab.selspot).activlum || !activlum->get_inconsistent();
                        // Tone Mapping
                        pedited->locallab.spots.at(pp->locallab.selspot).exptonemap = pedited->locallab.spots.at(pp->locallab.selspot).exptonemap || !exptonemap->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).stren = pedited->locallab.spots.at(pp->locallab.selspot).stren || stren->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).gamma = pedited->locallab.spots.at(pp->locallab.selspot).gamma || gamma->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).estop = pedited->locallab.spots.at(pp->locallab.selspot).estop || estop->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).scaltm = pedited->locallab.spots.at(pp->locallab.selspot).scaltm || scaltm->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).rewei = pedited->locallab.spots.at(pp->locallab.selspot).rewei || rewei->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensitm = pedited->locallab.spots.at(pp->locallab.selspot).sensitm || sensitm->getEditedState();
                        // Retinex
                        pedited->locallab.spots.at(pp->locallab.selspot).expreti = pedited->locallab.spots.at(pp->locallab.selspot).expreti || !expreti->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).retinexMethod = pedited->locallab.spots.at(pp->locallab.selspot).retinexMethod || retinexMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pedited->locallab.spots.at(pp->locallab.selspot).str = pedited->locallab.spots.at(pp->locallab.selspot).str || str->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).chrrt = pedited->locallab.spots.at(pp->locallab.selspot).chrrt || chrrt->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).neigh = pedited->locallab.spots.at(pp->locallab.selspot).neigh || neigh->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).vart = pedited->locallab.spots.at(pp->locallab.selspot).vart || vart->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).dehaz = pedited->locallab.spots.at(pp->locallab.selspot).dehaz || dehaz->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensih = pedited->locallab.spots.at(pp->locallab.selspot).sensih || sensih->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).localTgaincurve = pedited->locallab.spots.at(pp->locallab.selspot).localTgaincurve || !cTgainshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).inversret = pedited->locallab.spots.at(pp->locallab.selspot).inversret || !inversret->get_inconsistent();
                        // Sharpening
                        pedited->locallab.spots.at(pp->locallab.selspot).expsharp = pedited->locallab.spots.at(pp->locallab.selspot).expsharp || !expsharp->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).sharcontrast = pedited->locallab.spots.at(pp->locallab.selspot).sharcontrast || sharcontrast->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sharradius = pedited->locallab.spots.at(pp->locallab.selspot).sharradius || sharradius->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sharamount = pedited->locallab.spots.at(pp->locallab.selspot).sharamount || sharamount->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).shardamping = pedited->locallab.spots.at(pp->locallab.selspot).shardamping || shardamping->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).shariter = pedited->locallab.spots.at(pp->locallab.selspot).shariter || shariter->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sharblur = pedited->locallab.spots.at(pp->locallab.selspot).sharblur || sharblur->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensisha = pedited->locallab.spots.at(pp->locallab.selspot).sensisha || sensisha->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).inverssha = pedited->locallab.spots.at(pp->locallab.selspot).inverssha || !inverssha->get_inconsistent();
                        // Local Contrast
                        pedited->locallab.spots.at(pp->locallab.selspot).expcontrast = pedited->locallab.spots.at(pp->locallab.selspot).expcontrast || !expcontrast->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).lcradius = pedited->locallab.spots.at(pp->locallab.selspot).lcradius || lcradius->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).lcamount = pedited->locallab.spots.at(pp->locallab.selspot).lcamount || lcamount->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).lcdarkness = pedited->locallab.spots.at(pp->locallab.selspot).lcdarkness || lcdarkness->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).lclightness = pedited->locallab.spots.at(pp->locallab.selspot).lclightness || lclightness->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensilc = pedited->locallab.spots.at(pp->locallab.selspot).sensilc || sensilc->getEditedState();
                        // Contrast by detail levels
                        pedited->locallab.spots.at(pp->locallab.selspot).expcbdl = pedited->locallab.spots.at(pp->locallab.selspot).expcbdl || !expcbdl->get_inconsistent();

                        for (int i = 0; i < 5; i++) {
                            pedited->locallab.spots.at(pp->locallab.selspot).mult[i] = pedited->locallab.spots.at(pp->locallab.selspot).mult[i] || multiplier[i]->getEditedState();
                        }

                        pedited->locallab.spots.at(pp->locallab.selspot).chromacbdl = pedited->locallab.spots.at(pp->locallab.selspot).chromacbdl || chromacbdl->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).threshold = pedited->locallab.spots.at(pp->locallab.selspot).threshold || threshold->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensicb = pedited->locallab.spots.at(pp->locallab.selspot).sensicb || sensicb->getEditedState();
                        // Denoise
                        pedited->locallab.spots.at(pp->locallab.selspot).expdenoi = pedited->locallab.spots.at(pp->locallab.selspot).expdenoi || !expdenoi->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselumf = pedited->locallab.spots.at(pp->locallab.selspot).noiselumf || noiselumf->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselumc = pedited->locallab.spots.at(pp->locallab.selspot).noiselumc || noiselumc->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselumdetail = pedited->locallab.spots.at(pp->locallab.selspot).noiselumdetail || noiselumdetail->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselequal = pedited->locallab.spots.at(pp->locallab.selspot).noiselequal || noiselequal->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noisechrof = pedited->locallab.spots.at(pp->locallab.selspot).noisechrof || noisechrof->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noisechroc = pedited->locallab.spots.at(pp->locallab.selspot).noisechroc || noisechroc->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noisechrodetail = pedited->locallab.spots.at(pp->locallab.selspot).noisechrodetail || noisechrodetail->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).adjblur = pedited->locallab.spots.at(pp->locallab.selspot).adjblur || adjblur->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).bilateral = pedited->locallab.spots.at(pp->locallab.selspot).bilateral || bilateral->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensiden = pedited->locallab.spots.at(pp->locallab.selspot).sensiden || sensiden->getEditedState();
                        // Others
                        pedited->locallab.spots.at(pp->locallab.selspot).avoid = pedited->locallab.spots.at(pp->locallab.selspot).avoid || !avoid->get_inconsistent();
                    }
                }
            }
    }

    // Update Locallab tools GUI
    disableListener();
    updateSpecificGUIState();
    enableListener();
}

void Locallab::protectskins_toggled()
{
    // printf("protectskins_toggled\n");

    if (multiImage) {
        if (protectSkins->get_inconsistent()) {
            protectSkins->set_inconsistent(false);
            pskinsconn.block(true);
            protectSkins->set_active(false);
            pskinsconn.block(false);
        }
    }

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            if (protectSkins->get_active()) {
                listener->panelChanged(EvlocallabProtectSkins, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvlocallabProtectSkins, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::avoidcolorshift_toggled()
{
    // printf("avoidcolorshift_toggled\n");

    if (multiImage) {
        if (avoidColorShift->get_inconsistent()) {
            avoidColorShift->set_inconsistent(false);
            ashiftconn.block(true);
            avoidColorShift->set_active(false);
            ashiftconn.block(false);
        }
    }

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            if (avoidColorShift->get_active()) {
                listener->panelChanged(EvlocallabAvoidColorShift, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvlocallabAvoidColorShift, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::pastsattog_toggled()
{
    // printf("pastsattog_toggled\n");

    if (multiImage) {
        if (pastSatTog->get_inconsistent()) {
            pastSatTog->set_inconsistent(false);
            pastsattogconn.block(true);
            pastSatTog->set_active(false);
            pastsattogconn.block(false);
        }
    }

    // Update Vibrance GUI according to pastsattog button state (to be compliant with updateSpecificGUIState function)
    if (multiImage && pastSatTog->get_inconsistent()) {
        psThreshold->set_sensitive(true);
        saturated->set_sensitive(true);
    } else if (pastSatTog->get_active()) {
        // Link both slider, so we set saturated and psThresholds unsensitive
        psThreshold->set_sensitive(false);
        saturated->set_sensitive(false);
        saturated->setValue(pastels->getValue()); // Pastels and Saturated are linked
    } else {
        // Separate sliders, so we set saturated and psThresholds sensitive again
        psThreshold->set_sensitive(true);
        saturated->set_sensitive(true);
    }

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            if (pastSatTog->get_active()) {
                listener->panelChanged(EvlocallabPastSatTog, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvlocallabPastSatTog, M("GENERAL_DISABLED"));
            }
        }
    }
}


void Locallab::curveChanged(CurveEditor* ce)
{
    // Color & Light
    if (getEnabled() && expcolor->getEnabled()) {
        if (ce == llshape) {
            if (listener) {
                listener->panelChanged(Evlocallabllshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == ccshape) {
            if (listener) {
                listener->panelChanged(Evlocallabccshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == LHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLHshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == HHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == CCmaskshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == LLmaskshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

    }

    // Exposure
    if (getEnabled() && expexpose->getEnabled()) {
        if (ce == shapeexpos) {
            if (listener) {
                listener->panelChanged(Evlocallabshapeexpos, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == CCmaskexpshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskexpshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == LLmaskexpshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskexpshape, M("HISTORY_CUSTOMCURVE"));
            }
        }
        
    }

    // Vibrance
    if (getEnabled() && expvibrance->getEnabled()) {
        if (ce == skinTonesCurve) {
            if (listener) {
                listener->panelChanged(EvlocallabSkinTonesCurve, M("HISTORY_CUSTOMCURVE"));
            }
        }
    }

    // Retinex
    if (getEnabled() && expreti->getEnabled()) {
        if (ce == cTgainshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCTgainCurve, M("HISTORY_CUSTOMCURVE"));
            }
        }
    }
}

void Locallab::retinexMethodChanged()
{
    // printf("retinexMethodChanged\n");

    if (getEnabled() && expreti->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabretinexMethod, retinexMethod->get_active_text());
        }
    }
}

void Locallab::blurMethodChanged()
{
    // printf("blurMethodChanged\n");

    // Update Blur & Noise GUI according to blurMethod combobox (to be compliant with updateSpecificGUIState function)
    if (multiImage && blurMethod->get_active_text() == M("GENERAL_UNCHANGED")) {
        sensibn->show();
    } else if (blurMethod->get_active_row_number() == 0 || blurMethod->get_active_row_number() == 2) {
        sensibn->show();
    } else {
        sensibn->hide();
    }

    if (blurMethod->get_active_row_number() == 2) {
        strength->hide();
    } else {
        strength->show();
    }

    if (getEnabled() && expblur->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabblurMethod, blurMethod->get_active_text());
        }
    }
}



void Locallab::qualitycurveMethodChanged()
{
    // printf("qualitycurveMethodChanged\n");

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text());
        }
    } 
}

void Locallab::showmaskcolMethodChanged()
{
    if((showmaskcolMethod->get_active_row_number() == 1 || showmaskcolMethod->get_active_row_number() == 2 || showmaskcolMethod->get_active_row_number() == 4) && expcolor->getEnabled()) {
        showmaskexpMethod->set_active(0);
        expexpose->setEnabled(false); 
    }
    
    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskcolMethod , showmaskcolMethod->get_active_text());
        }
    }
}

void Locallab::showmaskexpMethodChanged()
{
    
    if((showmaskexpMethod->get_active_row_number() == 1 || showmaskexpMethod->get_active_row_number() == 2 || showmaskexpMethod->get_active_row_number() == 4) && expexpose->getEnabled()) {
        showmaskcolMethod->set_active(0);
        expcolor->setEnabled(false); 
    }
    
    if (getEnabled() && expexpose->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshowmaskexpMethod , showmaskexpMethod->get_active_text());
        }
    }
}

void Locallab::inversChanged()
{
    // printf("inversChanged\n");

    if (multiImage) {
        if (invers->get_inconsistent()) {
            invers->set_inconsistent(false);
            inversConn.block(true);
            invers->set_active(false);
            inversConn.block(false);
        }
    }

    // Update Color & Light GUI according to invers button state (to be compliant with updateSpecificGUIState function)
    if (multiImage && invers->get_inconsistent()) {
        sensi->show();
        llCurveEditorG->show();
        curvactiv->show();
        qualitycurveMethod->show();
        labqualcurv->show();
        maskcolFrame->show();
    } else if (invers->get_active()) {
        sensi->show();
        llCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        labqualcurv->hide();
        maskcolFrame->hide();
    } else {
        sensi->show();
        llCurveEditorG->show();
        curvactiv->show();
        qualitycurveMethod->show();
        labqualcurv->show();
        maskcolFrame->show();
    }

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            if (invers->get_active()) {
                listener->panelChanged(Evlocallabinvers, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinvers, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::curvactivChanged()
{
    // printf("curvactivChanged\n");

    if (multiImage) {
        if (curvactiv->get_inconsistent()) {
            curvactiv->set_inconsistent(false);
            curvactivConn.block(true);
            curvactiv->set_active(false);
            curvactivConn.block(false);
        }
    }

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            if (curvactiv->get_active()) {
                listener->panelChanged(Evlocallabcurvactiv, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabcurvactiv, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::activlumChanged()
{
    // printf("activlumChanged\n");

    if (multiImage) {
        if (activlum->get_inconsistent()) {
            activlum->set_inconsistent(false);
            activlumConn.block(true);
            activlum->set_active(false);
            activlumConn.block(false);
        }
    }

    if (getEnabled() && expblur->getEnabled()) {
        if (listener) {
            if (activlum->get_active()) {
                listener->panelChanged(Evlocallabactivlum, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabactivlum, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::inversshaChanged()
{
    // printf("inversshaChanged\n");

    if (multiImage) {
        if (inverssha->get_inconsistent()) {
            inverssha->set_inconsistent(false);
            inversshaConn.block(true);
            inverssha->set_active(false);
            inversshaConn.block(false);
        }
    }

    // Update Sharpening GUI according to inverssha button state (to be compliant with updateSpecificGUIState function)
    if (multiImage && inverssha->get_inconsistent()) {
        sensisha->show();
    } else if (inverssha->get_active()) {
        sensisha->hide();
    } else {
        sensisha->show();
    }

    if (getEnabled() && expsharp->getEnabled()) {
        if (listener) {
            if (inverssha->get_active()) {
                listener->panelChanged(Evlocallabinverssha, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinverssha, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::inversretChanged()
{
    // printf("inversretChanged\n");

    if (multiImage) {
        if (inversret->get_inconsistent()) {
            inversret->set_inconsistent(false);
            inversretConn.block(true);
            inversret->set_active(false);
            inversretConn.block(false);
        }
    }

    // Update Retinex GUI according to inversret button state (to be compliant with updateSpecificGUIState function)
    if (multiImage && inversret->get_inconsistent()) {
        sensih->show();
    } else if (inversret->get_active()) {
        sensih->hide();
        dehaz->hide();
    } else {
        sensih->show();
        dehaz->show();
    }

    if (getEnabled() && expreti->getEnabled()) {
        if (listener) {
            if (inversret->get_active()) {
                listener->panelChanged(Evlocallabinversret, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinversret, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::setParamEditable(bool cond)
{
    // printf("setParamEditable: %d\n", cond);

    // Update params editable state for controlspotpanel
    expsettings->setParamEditable(cond);

    // Color & Light
    expcolor->set_sensitive(cond);
    // Exposure
    expexpose->set_sensitive(cond);
    // Vibrance
    expvibrance->set_sensitive(cond);
    // Soft Light
    expsoft->set_sensitive(cond);
    // Lab Region
    explabregion->set_sensitive(cond);
    // Blur & Noise
    expblur->set_sensitive(cond);
    // Tone Mapping
    exptonemap->set_sensitive(cond);
    // Retinex
    expreti->set_sensitive(cond);
    // Sharpening
    expsharp->set_sensitive(cond);
    // Local Contrast
    expcontrast->set_sensitive(cond);
    // Contrast by detail levels
    expcbdl->set_sensitive(cond);
    // Denoise
    expdenoi->set_sensitive(cond);
    // Others
    avoid->set_sensitive(cond);
}

void Locallab::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    defparams = defParams;
    defpedited = pedited;

    if (pedited) {
        pe = new ParamsEdited(*pedited);
    } else {
        pe = nullptr;
    }
}

void Locallab::setDefaults(const ProcParams * defParams, const ParamsEdited * pedited, int id)
{
    // printf("setDefaults\n");

    // Set default values and edited states for controlspotpanel
    expsettings->setDefaults(defParams, pedited, id);

    // Find vector index of given spot id (index = -1 if not found)
    int index = -1;

    for (int i = 0; i < (int)defParams->locallab.spots.size(); i++) {
        if (defParams->locallab.spots.at(i).id == id) {
            index = i;

            break;
        }
    }

    // Set default values for adjusters
    const LocallabParams::LocallabSpot* defSpot = new LocallabParams::LocallabSpot();

    if (index != -1 && index < (int)defParams->locallab.spots.size()) {
        defSpot = &defParams->locallab.spots.at(index);
    }

    // Color & Light
    lightness->setDefault((double)defSpot->lightness);
    contrast->setDefault((double)defSpot->contrast);
    chroma->setDefault((double)defSpot->chroma);
    sensi->setDefault((double)defSpot->sensi);
    // Exposure
    expcomp->setDefault((double)defSpot->expcomp);
    hlcompr->setDefault((double)defSpot->hlcompr);
    hlcomprthresh->setDefault((double)defSpot->hlcomprthresh);
    black->setDefault((double)defSpot->black);
    shcompr->setDefault((double)defSpot->shcompr);
    warm->setDefault((double)defSpot->warm);
    sensiex->setDefault((double)defSpot->sensiex);
    // Vibrance
    saturated->setDefault((double)defSpot->saturated);
    pastels->setDefault((double)defSpot->pastels);
    psThreshold->setDefault<int>(defSpot->psthreshold);
    sensiv->setDefault((double)defSpot->sensiv);
    // Soft Light
    streng->setDefault((double)defSpot->streng);
    sensisf->setDefault((double)defSpot->sensisf);
    // Blur & Noise
    radius->setDefault((double)defSpot->radius);
    strength->setDefault((double)defSpot->strength);
    sensibn->setDefault((double)defSpot->sensibn);
    // Tone Mapping
    stren->setDefault((double)defSpot->stren);
    gamma->setDefault((double)defSpot->gamma);
    estop->setDefault((double)defSpot->estop);
    scaltm->setDefault((double)defSpot->scaltm);
    rewei->setDefault((double)defSpot->rewei);
    sensitm->setDefault((double)defSpot->sensitm);
    // Retinex
    str->setDefault((double)defSpot->str);
    chrrt->setDefault((double)defSpot->chrrt);
    neigh->setDefault((double)defSpot->neigh);
    vart->setDefault((double)defSpot->vart);
    dehaz->setDefault((double)defSpot->dehaz);
    sensih->setDefault((double)defSpot->sensih);
    // Sharpening
    sharcontrast->setDefault((double)defSpot->sharcontrast);
    sharradius->setDefault((double)defSpot->sharradius);
    sharamount->setDefault((double)defSpot->sharamount);
    shardamping->setDefault((double)defSpot->shardamping);
    shariter->setDefault((double)defSpot->shariter);
    sharblur->setDefault((double)defSpot->sharblur);
    sensisha->setDefault((double)defSpot->sensisha);
    // Local Contrast
    lcradius->setDefault((double)defSpot->lcradius);
    lcamount->setDefault((double)defSpot->lcamount);
    lcdarkness->setDefault((double)defSpot->lcdarkness);
    lclightness->setDefault((double)defSpot->lclightness);
    sensilc->setDefault((double)defSpot->sensilc);
    // Contrast by detail levels
    for (int i = 0; i < 5; i++) {
        multiplier[i]->setDefault(defSpot->mult[i]);
    }
    chromacbdl->setDefault((double)defSpot->chromacbdl);
    threshold->setDefault(defSpot->threshold);
    sensicb->setDefault((double)defSpot->sensicb);
    // Denoise
    noiselumf->setDefault((double)defSpot->noiselumf);
    noiselumc->setDefault((double)defSpot->noiselumc);
    noiselumdetail->setDefault((double)defSpot->noiselumdetail);
    noiselequal->setDefault((double)defSpot->noiselequal);
    noisechrof->setDefault((double)defSpot->noisechrof);
    noisechroc->setDefault((double)defSpot->noisechroc);
    noisechrodetail->setDefault((double)defSpot->noisechrodetail);
    adjblur->setDefault((double)defSpot->adjblur);
    bilateral->setDefault((double)defSpot->bilateral);
    sensiden->setDefault((double)defSpot->sensiden);

    // Set default edited states for adjusters and threshold adjusters
    if (!pedited) {
        // Color & Light
        lightness->setDefaultEditedState(Irrelevant);
        contrast->setDefaultEditedState(Irrelevant);
        chroma->setDefaultEditedState(Irrelevant);
        sensi->setDefaultEditedState(Irrelevant);
        // Exposure
        expcomp->setDefaultEditedState(Irrelevant);
        hlcompr->setDefaultEditedState(Irrelevant);
        hlcomprthresh->setDefaultEditedState(Irrelevant);
        black->setDefaultEditedState(Irrelevant);
        shcompr->setDefaultEditedState(Irrelevant);
        warm->setDefaultEditedState(Irrelevant);
        sensiex->setDefaultEditedState(Irrelevant);
        // Vibrance
        saturated->setDefaultEditedState(Irrelevant);
        pastels->setDefaultEditedState(Irrelevant);
        psThreshold->setDefaultEditedState(Irrelevant);
        sensiv->setDefaultEditedState(Irrelevant);
        // Soft Light
        streng->setDefaultEditedState(Irrelevant);
        sensisf->setDefaultEditedState(Irrelevant);
        // Blur & Noise
        radius->setDefaultEditedState(Irrelevant);
        strength->setDefaultEditedState(Irrelevant);
        sensibn->setDefaultEditedState(Irrelevant);
        // Tone Mapping
        stren->setDefaultEditedState(Irrelevant);
        gamma->setDefaultEditedState(Irrelevant);
        estop->setDefaultEditedState(Irrelevant);
        scaltm->setDefaultEditedState(Irrelevant);
        rewei->setDefaultEditedState(Irrelevant);
        sensitm->setDefaultEditedState(Irrelevant);
        // Retinex
        str->setDefaultEditedState(Irrelevant);
        chrrt->setDefaultEditedState(Irrelevant);
        neigh->setDefaultEditedState(Irrelevant);
        vart->setDefaultEditedState(Irrelevant);
        dehaz->setDefaultEditedState(Irrelevant);
        sensih->setDefaultEditedState(Irrelevant);
        // Sharpening
        sharcontrast->setDefaultEditedState(Irrelevant);
        sharradius->setDefaultEditedState(Irrelevant);
        sharamount->setDefaultEditedState(Irrelevant);
        shardamping->setDefaultEditedState(Irrelevant);
        shariter->setDefaultEditedState(Irrelevant);
        sharblur->setDefaultEditedState(Irrelevant);
        sensisha->setDefaultEditedState(Irrelevant);
        // Local Contrast
        lcradius->setDefaultEditedState(Irrelevant);
        lcamount->setDefaultEditedState(Irrelevant);
        lcdarkness->setDefaultEditedState(Irrelevant);
        lclightness->setDefaultEditedState(Irrelevant);
        sensilc->setDefaultEditedState(Irrelevant);
        // Contrast by detail levels
        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState(Irrelevant);
        }
        chromacbdl->setDefaultEditedState(Irrelevant);
        threshold->setDefaultEditedState(Irrelevant);
        sensicb->setDefaultEditedState(Irrelevant);
        // Denoise
        noiselumf->setDefaultEditedState(Irrelevant);
        noiselumc->setDefaultEditedState(Irrelevant);
        noiselumdetail->setDefaultEditedState(Irrelevant);
        noiselequal->setDefaultEditedState(Irrelevant);
        noisechrof->setDefaultEditedState(Irrelevant);
        noisechroc->setDefaultEditedState(Irrelevant);
        noisechrodetail->setDefaultEditedState(Irrelevant);
        adjblur->setDefaultEditedState(Irrelevant);
        bilateral->setDefaultEditedState(Irrelevant);
        sensiden->setDefaultEditedState(Irrelevant);
    } else {
        const LocallabParamsEdited::LocallabSpotEdited* defSpotState = new LocallabParamsEdited::LocallabSpotEdited(true);

        if (index != 1 && index < (int)pedited->locallab.spots.size()) {
            defSpotState = &pedited->locallab.spots.at(index);
        }

        // Color & Light
        lightness->setDefaultEditedState(defSpotState->lightness ? Edited : UnEdited);
        contrast->setDefaultEditedState(defSpotState->contrast ? Edited : UnEdited);
        chroma->setDefaultEditedState(defSpotState->chroma ? Edited : UnEdited);
        sensi->setDefaultEditedState(defSpotState->sensi ? Edited : UnEdited);
        // Exposure
        expcomp->setDefaultEditedState(defSpotState->expcomp ? Edited : UnEdited);
        hlcompr->setDefaultEditedState(defSpotState->hlcompr ? Edited : UnEdited);
        hlcomprthresh->setDefaultEditedState(defSpotState->hlcomprthresh ? Edited : UnEdited);
        black->setDefaultEditedState(defSpotState->black ? Edited : UnEdited);
        shcompr->setDefaultEditedState(defSpotState->shcompr ? Edited : UnEdited);
        warm->setDefaultEditedState(defSpotState->warm ? Edited : UnEdited);
        sensiex->setDefaultEditedState(defSpotState->sensiex ? Edited : UnEdited);
        // Vibrance
        saturated->setDefaultEditedState(defSpotState->saturated ? Edited : UnEdited);
        pastels->setDefaultEditedState(defSpotState->pastels ? Edited : UnEdited);
        psThreshold->setDefaultEditedState(defSpotState->psthreshold ? Edited : UnEdited);
        sensiv->setDefaultEditedState(defSpotState->sensiv ? Edited : UnEdited);
        // Soft Light
        streng->setDefaultEditedState(defSpotState->streng ? Edited : UnEdited);
        sensisf->setDefaultEditedState(defSpotState->sensisf ? Edited : UnEdited);
        // Blur & Noise
        radius->setDefaultEditedState(defSpotState->radius ? Edited : UnEdited);
        strength->setDefaultEditedState(defSpotState->strength ? Edited : UnEdited);
        sensibn->setDefaultEditedState(defSpotState->sensibn ? Edited : UnEdited);
        // Tone Mapping
        stren->setDefaultEditedState(defSpotState->stren ? Edited : UnEdited);
        gamma->setDefaultEditedState(defSpotState->gamma ? Edited : UnEdited);
        estop->setDefaultEditedState(defSpotState->estop ? Edited : UnEdited);
        scaltm->setDefaultEditedState(defSpotState->scaltm ? Edited : UnEdited);
        rewei->setDefaultEditedState(defSpotState->rewei ? Edited : UnEdited);
        sensitm->setDefaultEditedState(defSpotState->sensitm ? Edited : UnEdited);
        // Retinex
        str->setDefaultEditedState(defSpotState->str ? Edited : UnEdited);
        chrrt->setDefaultEditedState(defSpotState->chrrt ? Edited : UnEdited);
        neigh->setDefaultEditedState(defSpotState->neigh ? Edited : UnEdited);
        vart->setDefaultEditedState(defSpotState->vart ? Edited : UnEdited);
        dehaz->setDefaultEditedState(defSpotState->dehaz ? Edited : UnEdited);
        sensih->setDefaultEditedState(defSpotState->sensih ? Edited : UnEdited);
        // Sharpening
        sharcontrast->setDefaultEditedState(defSpotState->sharcontrast ? Edited : UnEdited);
        sharradius->setDefaultEditedState(defSpotState->sharradius ? Edited : UnEdited);
        sharamount->setDefaultEditedState(defSpotState->sharamount ? Edited : UnEdited);
        shardamping->setDefaultEditedState(defSpotState->shardamping ? Edited : UnEdited);
        shariter->setDefaultEditedState(defSpotState->shariter ? Edited : UnEdited);
        sharblur->setDefaultEditedState(defSpotState->sharblur ? Edited : UnEdited);
        sensisha->setDefaultEditedState(defSpotState->sensisha ? Edited : UnEdited);
        // Local Contrast
        lcradius->setDefaultEditedState(defSpotState->lcradius ? Edited : UnEdited);
        lcamount->setDefaultEditedState(defSpotState->lcamount ? Edited : UnEdited);
        lcdarkness->setDefaultEditedState(defSpotState->lcdarkness ? Edited : UnEdited);
        lclightness->setDefaultEditedState(defSpotState->lclightness ? Edited : UnEdited);
        sensilc->setDefaultEditedState(defSpotState->sensilc ? Edited : UnEdited);
        // Contrast by detail levels
        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState(defSpotState->mult[i] ? Edited : UnEdited);
        }
        chromacbdl->setDefaultEditedState(defSpotState->chromacbdl ? Edited : UnEdited);
        threshold->setDefaultEditedState(defSpotState->threshold ? Edited : UnEdited);
        sensicb->setDefaultEditedState(defSpotState->sensicb ? Edited : UnEdited);
        // Denoise
        noiselumf->setDefaultEditedState(defSpotState->noiselumf ? Edited : UnEdited);
        noiselumc->setDefaultEditedState(defSpotState->noiselumc ? Edited : UnEdited);
        noiselumdetail->setDefaultEditedState(defSpotState->noiselumdetail ? Edited : UnEdited);
        noiselequal->setDefaultEditedState(defSpotState->noiselequal ? Edited : UnEdited);
        noisechrof->setDefaultEditedState(defSpotState->noisechrof ? Edited : UnEdited);
        noisechroc->setDefaultEditedState(defSpotState->noisechroc ? Edited : UnEdited);
        noisechrodetail->setDefaultEditedState(defSpotState->noisechrodetail ? Edited : UnEdited);
        adjblur->setDefaultEditedState(defSpotState->adjblur ? Edited : UnEdited);
        bilateral->setDefaultEditedState(defSpotState->bilateral ? Edited : UnEdited);
        sensiden->setDefaultEditedState(defSpotState->sensiden ? Edited : UnEdited);
    }
}

void Locallab::adjusterAutoToggled(Adjuster* a, bool newval)
{
    // Not used
}

void Locallab::adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop)
{
    // Not used
}
void Locallab::adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight)
{
    // Not used
}

void Locallab::adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight)
{
    // Not used
}

void Locallab::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
    // Not used
}

void Locallab::adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop)
{
    // printf("adjusterChangedTS\n");

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabPastSatThreshold, psThreshold->getHistoryString());
        }
    }
}

void Locallab::adjusterChanged(Adjuster * a, double newval)
{
    // printf("adjusterChanged\n");

    // Color & Light
    if (getEnabled() && expcolor->getEnabled()) {
        if (a == lightness) {
            if (listener) {
                listener->panelChanged(Evlocallablightness, lightness->getTextValue());
            }
        }

        if (a == contrast) {
            if (listener) {
                listener->panelChanged(Evlocallabcontrast, contrast->getTextValue());
            }
        }

        if (a == chroma) {
            if (listener) {
                listener->panelChanged(Evlocallabchroma, chroma->getTextValue());
            }
        }

        if (a == sensi) {
            if (listener) {
                listener->panelChanged(Evlocallabsensi, sensi->getTextValue());
            }
        }
    }

    // Exposure
    if (a == black) {
        // Update Exposure GUI according to black adjuster state (to be compliant with updateSpecificGUIState function)
        if (multiImage && black->getEditedState() != Edited) {
            shcompr->set_sensitive(true);
        } else {
            shcompr->set_sensitive(!((int)black->getValue() == 0)); // At black = 0, shcompr value has no effect
        }
    }

    if (getEnabled() && expexpose->getEnabled()) {
        if (a == expcomp) {
            if (listener) {
                listener->panelChanged(Evlocallabexpcomp, expcomp->getTextValue());
            }
        }

        if (a == hlcompr) {
            if (listener) {
                listener->panelChanged(Evlocallabhlcompr, hlcompr->getTextValue());
            }
        }

        if (a == hlcomprthresh) {
            if (listener) {
                listener->panelChanged(Evlocallabhlcomprthresh, hlcomprthresh->getTextValue());
            }
        }

        if (a == black) {
            if (listener) {
                listener->panelChanged(Evlocallabblack, black->getTextValue());
            }
        }

        if (a == shcompr) {
            if (listener) {
                listener->panelChanged(Evlocallabshcompr, shcompr->getTextValue());
            }
        }

        if (a == warm) {
            if (listener) {
                listener->panelChanged(Evlocallabwarm, warm->getTextValue());
            }
        }

        if (a == sensiex) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiex, sensiex->getTextValue());
            }
        }
    }

    // Vibrance
    if (a == pastels && pastSatTog->get_active() && !(multiImage && pastSatTog->get_inconsistent())) {
        saturated->setValue(newval);
    }

    if (getEnabled() && expvibrance->getEnabled()) {
        if (a == saturated && !pastSatTog->get_active()) {
            if (listener) {
                listener->panelChanged(EvlocallabSaturated, saturated->getTextValue());
            }
        }

        if (a == pastels) {
            if (listener) {
                listener->panelChanged(EvlocallabPastels, pastels->getTextValue());
            }
        }

        if (a == sensiv) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiv, sensiv->getTextValue());
            }
        }
    }

    // Soft Light
    if (getEnabled() && expsoft->getEnabled()) {
        if (a == streng) {
            if (listener) {
                listener->panelChanged(Evlocallabstreng, streng->getTextValue());
            }
        }

        if (a == sensisf) {
            if (listener) {
                listener->panelChanged(Evlocallabsensisf, sensisf->getTextValue());
            }
        }
    }

    //Lab region
    /*
    if (getEnabled() && explabregion->getEnabled()) {
        if (a == labRegionSlope) {
            if (listener) {
                listener->panelChanged(EvlocallablabRegionSlope, labRegionSlope->getTextValue());
            }
        }

    }
    */

    // Blur & Noise
    if (getEnabled() && expblur->getEnabled()) {
        if (a == radius) {
            if (listener) {
                listener->panelChanged(Evlocallabradius, radius->getTextValue());
            }
        }

        if (a == strength) {
            if (listener) {
                listener->panelChanged(Evlocallabstrength, strength->getTextValue());
            }
        }

        if (a == sensibn) {
            if (listener) {
                listener->panelChanged(Evlocallabsensibn, sensibn->getTextValue());
            }
        }
    }

    // Tone Mapping
    if (getEnabled() && exptonemap->getEnabled()) {
        if (a == stren) {
            if (listener) {
                listener->panelChanged(Evlocallabstren, stren->getTextValue());
            }
        }

        if (a == gamma) {
            if (listener) {
                listener->panelChanged(Evlocallabgamma, gamma->getTextValue());
            }
        }

        if (a == estop) {
            if (listener) {
                listener->panelChanged(Evlocallabestop, estop->getTextValue());
            }
        }

        if (a == scaltm) {
            if (listener) {
                listener->panelChanged(Evlocallabscaltm, scaltm->getTextValue());
            }
        }

        if (a == rewei) {
            if (listener) {
                listener->panelChanged(Evlocallabrewei, rewei->getTextValue());
            }
        }

        if (a == sensitm) {
            if (listener) {
                listener->panelChanged(Evlocallabsensitm, sensitm->getTextValue());
            }
        }
    }

    // Retinex
    if (getEnabled() && expreti->getEnabled()) {
        if (a == str) {
            if (listener) {
                listener->panelChanged(Evlocallabstr, str->getTextValue());
            }
        }

        if (a == chrrt) {
            if (listener) {
                listener->panelChanged(Evlocallabchrrt, chrrt->getTextValue());
            }
        }

        if (a == neigh) {
            if (listener) {
                listener->panelChanged(Evlocallabneigh, neigh->getTextValue());
            }
        }

        if (a == vart) {
            if (listener) {
                listener->panelChanged(Evlocallabvart, vart->getTextValue());
            }
        }

        if (a == dehaz) {
            if (listener) {
                listener->panelChanged(Evlocallabdehaz, dehaz->getTextValue());
            }
        }

        if (a == sensih) {
            if (listener) {
                listener->panelChanged(Evlocallabsensih, sensih->getTextValue());
            }
        }
    }

    // Sharpening
    if (getEnabled() && expsharp->getEnabled()) {
        if (a == sharcontrast) {
            if (listener) {
                listener->panelChanged(Evlocallabsharcontrast, sharcontrast->getTextValue());
            }
        }

        if (a == sharradius) {
            if (listener) {
                listener->panelChanged(Evlocallabsharradius, sharradius->getTextValue());
            }
        }

        if (a == sharamount) {
            if (listener) {
                listener->panelChanged(Evlocallabsharamount, sharamount->getTextValue());
            }
        }

        if (a == shardamping) {
            if (listener) {
                listener->panelChanged(Evlocallabshardamping, shardamping->getTextValue());
            }
        }

        if (a == shariter) {
            if (listener) {
                listener->panelChanged(Evlocallabshariter, shariter->getTextValue());
            }
        }

        if (a == sharblur) {
            if (listener) {
                listener->panelChanged(Evlocallabsharblur, sharblur->getTextValue());
            }
        }

        if (a == sensisha) {
            if (listener) {
                listener->panelChanged(Evlocallabsensis, sensisha->getTextValue());
            }
        }
    }

    // Local Contrast
    if (getEnabled() && expcontrast->getEnabled()) {
        if (a == lcradius) {
            if (listener) {
                listener->panelChanged(Evlocallablcradius, lcradius->getTextValue());
            }
        }

        if (a == lcamount) {
            if (listener) {
                listener->panelChanged(Evlocallablcamount, lcamount->getTextValue());
            }
        }

        if (a == lcdarkness) {
            if (listener) {
                listener->panelChanged(Evlocallablcdarkness, lcdarkness->getTextValue());
            }
        }

        if (a == lclightness) {
            if (listener) {
                listener->panelChanged(Evlocallablclightness, lclightness->getTextValue());
            }
        }

        if (a == sensilc) {
            if (listener) {
                listener->panelChanged(Evlocallabsensilc, sensilc->getTextValue());
            }
        }
    }

    // Contrast by detail levels
    if (getEnabled() && expcbdl->getEnabled()) {
        if (a == multiplier[0] || a == multiplier[1] || a == multiplier[2] || a == multiplier[3] || a == multiplier[4]) {
            if (listener) {
                listener->panelChanged(EvlocallabEqualizer,
                                       Glib::ustring::compose("%1, %2, %3, %4, %5",
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[0]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[1]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[2]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[3]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[4]->getValue())));
            }
        }

        if (a == chromacbdl) {
            if (listener) {
                listener->panelChanged(Evlocallabchromacbdl, chromacbdl->getTextValue());
            }
        }

        if (a == threshold) {
            if (listener) {
                listener->panelChanged(EvlocallabThresho, threshold->getTextValue());
            }
        }

        if (a == sensicb) {
            if (listener) {
                listener->panelChanged(Evlocallabsensicb, sensicb->getTextValue());
            }
        }
    }

    // Denoise
    if (getEnabled() && expdenoi->getEnabled()) {
        if (a == noiselumf) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf, noiselumf->getTextValue());
            }
        }

        if (a == noiselumc) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumc, noiselumc->getTextValue());
            }
        }

        if (a == noiselumdetail) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumdetail, noiselumdetail->getTextValue());
            }
        }

        if (a == noiselequal) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselequal, noiselequal->getTextValue());
            }
        }

        if (a == noisechrof) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechrof, noisechrof->getTextValue());
            }
        }

        if (a == noisechroc) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechroc, noisechroc->getTextValue());
            }
        }

        if (a == noisechrodetail) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechrodetail, noisechrodetail->getTextValue());
            }
        }

        if (a == adjblur) {
            if (listener) {
                listener->panelChanged(Evlocallabadjblur, adjblur->getTextValue());
            }
        }

        if (a == bilateral) {
            if (listener) {
                listener->panelChanged(Evlocallabbilateral, bilateral->getTextValue());
            }
        }

        if (a == sensiden) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiden, sensiden->getTextValue());
            }
        }
    }
}

void Locallab::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvlocallabEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvlocallabEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvlocallabEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::avoidChanged()
{
    // printf("avoidChanged\n");

    if (multiImage) {
        if (avoid->get_inconsistent()) {
            avoid->set_inconsistent(false);
            avoidConn.block(true);
            avoid->set_active(false);
            avoidConn.block(false);
        }
    }

    if (getEnabled()) {
        if (listener) {
            if (avoid->get_active()) {
                listener->panelChanged(Evlocallabavoid, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabavoid, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::trimValues(rtengine::procparams::ProcParams * pp)
{
    // TODO
}

void Locallab::setBatchMode(bool batchMode)
{
    // printf("BatchMode : %d\n", batchMode);

    ToolPanel::setBatchMode(batchMode);

    // Set batch mode for controlspotpanel
    expsettings->setBatchMode(batchMode);

    // Set batch mode for adjusters and threshold adjusters
    // Color & Light
    lightness->showEditedCB();
    contrast->showEditedCB();
    chroma->showEditedCB();
    sensi->showEditedCB();
    // Exposure
    expcomp->showEditedCB();
    hlcompr->showEditedCB();
    hlcomprthresh->showEditedCB();
    black->showEditedCB();
    shcompr->showEditedCB();
    warm->showEditedCB();
    sensiex->showEditedCB();
    // Vibrance
    saturated->showEditedCB();
    pastels->showEditedCB();
    psThreshold->showEditedCB();
    sensiv->showEditedCB();
    // Soft Light
    streng->showEditedCB();
    sensisf->showEditedCB();
    // Blur & Noise
    radius->showEditedCB();
    streng->showEditedCB();
    strength->showEditedCB();
    sensibn->showEditedCB();
    // Tone Mapping
    stren->showEditedCB();
    gamma->showEditedCB();
    estop->showEditedCB();
    scaltm->showEditedCB();
    rewei->showEditedCB();
    sensitm->showEditedCB();
    // Retinex
    str->showEditedCB();
    chrrt->showEditedCB();
    neigh->showEditedCB();
    vart->showEditedCB();
    dehaz->showEditedCB();
    sensih->showEditedCB();
    // Sharpening
    sharradius->showEditedCB();
    sharamount->showEditedCB();
    shardamping->showEditedCB();
    shariter->showEditedCB();
    sensisha->showEditedCB();
    // Local Contrast
    lcradius->showEditedCB();
    lcamount->showEditedCB();
    lcdarkness->showEditedCB();
    lclightness->showEditedCB();
    sensilc->showEditedCB();
    // Contrast by detail levels
    for (int i = 0; i < 5; i++) {
        multiplier[i]->showEditedCB();
    }

    chromacbdl->showEditedCB();
    threshold->showEditedCB();
    sensicb->showEditedCB();
    // Denoise
    noiselumf->showEditedCB();
    noiselumc->showEditedCB();
    noiselumdetail->showEditedCB();
    noiselequal->showEditedCB();
    noiselumf->showEditedCB();
    noisechroc->showEditedCB();
    noisechrodetail->showEditedCB();
    adjblur->showEditedCB();
    bilateral->showEditedCB();
    sensiden->showEditedCB();

    // Set batch mode for comboBoxText
    // Color & Light
    qualitycurveMethod->append(M("GENERAL_UNCHANGED"));
    showmaskcolMethod->append(M("GENERAL_UNCHANGED"));
    //Exposure
    showmaskexpMethod->append(M("GENERAL_UNCHANGED"));
    // Blur & Noise
    blurMethod->append(M("GENERAL_UNCHANGED"));
    // Retinex
    retinexMethod->append(M("GENERAL_UNCHANGED"));
}

std::vector<double> Locallab::getCurvePoints(ThresholdSelector* tAdjuster) const
{
    std::vector<double> points;
    double threshold, transitionWeighting;
    tAdjuster->getPositions<double> (transitionWeighting, threshold); // ( range -100;+100,   range 0;+100 )
    transitionWeighting /= 100.; // range -1., +1.
    threshold /= 100.;      // range  0., +1.

    // Initial point
    points.push_back(0.);
    points.push_back(0.);

    double p2 = 3.0 * threshold / 4.0;             // same one than in ipvibrance.cc
    double s0 = threshold + (1.0 - threshold) / 4.0; // same one than in ipvibrance.cc

    // point at the beginning of the first linear transition
    points.push_back(p2);
    points.push_back(0.);

    // Y value of the chroma mean point, calculated to get a straight line between p2 and s0
    double chromaMean = (threshold / 4.0) / (s0 - p2);

    // move chromaMean up or down depending on transitionWeighting
    if (transitionWeighting > 0.0) {
        // positive values -> give more weight to Saturated
        chromaMean = (1.0 - chromaMean) * transitionWeighting + chromaMean;
    } else if (transitionWeighting < 0.0) {
        // negative values -> give more weight to Pastels
        chromaMean =      chromaMean  * transitionWeighting + chromaMean;
    }

    // point at the location of the Top cursor, at the end of the first linear transition and the beginning of the second one
    points.push_back(threshold);
    points.push_back(chromaMean);

    if (threshold < 1.0) {

        // point at the end of the second linear transition
        points.push_back(s0);
        points.push_back(1.0);

        // end point
        points.push_back(1.0);
        points.push_back(1.0);
    }

    return points;
}

void Locallab::setEditProvider(EditDataProvider * provider)
{
    expsettings->setEditProvider(provider);
}

void Locallab::subscribe()
{
    expsettings->subscribe();
}

void Locallab::unsubscribe()
{
    expsettings->unsubscribe();
}

void Locallab::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R, G, B;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve

        Color::hsv2rgb01(float (valX), float (valY), 0.5f, R, G, B);
    } else if (callerId == 2) {  // cc - bottom bar

        float value = (1.f - 0.7f) * float (valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float (valY), float (valX), value, R, G, B);
    } else if (callerId == 3) {  // lc - bottom bar

        float value = (1.f - 0.7f) * float (valX) + 0.7f;
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float (valY), float (valX), value, R, G, B);
    } else if (callerId == 4) {  // LH - bottom bar
        Color::hsv2rgb01(float (valX), 0.5f, float (valY), R, G, B);
    } else if (callerId == 5) {  // HH - bottom bar
        float h = float ((valY - 0.5) * 0.3 + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01(h, 0.5f, 0.5f, R, G, B);
    } else if (callerId == 6) {
        // TODO
        float x = valX - 1.f/6.f;
        if (x < 0.f) {
            x += 1.f;
        }
        x = log2lin(x, 3.f);
        // float x = valX;
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);        
    } else if (callerId == 7) {
        Color::hsv2rgb01(float(valY), float(valX), 0.5f, R, G, B);        
    }

    caller->ccRed = double (R);
    caller->ccGreen = double (G);
    caller->ccBlue = double (B);
}

void Locallab::setListener(ToolPanelListener* tpl)
{
    this->listener = tpl;
    expsettings->setListener(tpl);
}

void Locallab::enableListener()
{
    // printf("enableListener\n");

    FoldableToolPanel::enableListener();
    // Color & Light
    enablecolorConn.block(false);
    curvactivConn.block(false);
    qualitycurveMethodConn.block(false);
    inversConn.block(false);
    showmaskcolMethodConn.block(false);
    // Exposure
    enableexposeConn.block(false);
    showmaskexpMethodConn.block(false);
    // Vibrance
    enablevibranceConn.block(false);
    pskinsconn.block(false);
    ashiftconn.block(false);
    pastsattogconn.block(false);
    // Soft Light
    enablesoftConn.block(false);
    // Lab Region
    enablelabregionConn.block(false);
    // Blur & Noise
    enableblurConn.block(false);
    blurMethodConn.block(false);
    activlumConn.block(false);
    // Tone Mapping
    enabletonemapConn.block(false);
    // Retinex
    enableretiConn.block(false);
    retinexMethodConn.block(false);
    inversretConn.block(false);
    // Sharpening
    enablesharpConn.block(false);
    inversshaConn.block(false);
    // Local Contrast
    enablecontrastConn.block(false);
    // Contrast by detail levels
    enablecbdlConn.block(false);
    // Denoise
    enabledenoiConn.block(false);
    avoidConn.block(false);
}

void Locallab::disableListener()
{
    // printf("disableListener\n");

    FoldableToolPanel::disableListener();
    // Color & Light
    enablecolorConn.block(true);
    curvactivConn.block(true);
    qualitycurveMethodConn.block(true);
    showmaskcolMethodConn.block(true);
    inversConn.block(true);
    // Exposure
    enableexposeConn.block(true);
    showmaskexpMethodConn.block(true);
    // Vibrance
    enablevibranceConn.block(true);
    pskinsconn.block(true);
    ashiftconn.block(true);
    pastsattogconn.block(true);
    // Soft Light
    enablesoftConn.block(true);
    // Lab Region
    enablelabregionConn.block(true);
    // Blur & Noise
    enableblurConn.block(true);
    blurMethodConn.block(true);
    activlumConn.block(true);
    // Tone Mapping
    enabletonemapConn.block(true);
    // Retinex
    enableretiConn.block(true);
    retinexMethodConn.block(true);
    inversretConn.block(true);
    // Sharpening
    enablesharpConn.block(true);
    inversshaConn.block(true);
    // Local Contrast
    enablecontrastConn.block(true);
    // Contrast by detail levels
    enablecbdlConn.block(true);
    // Denoise
    enabledenoiConn.block(true);
    avoidConn.block(true);
}

void Locallab::updateLocallabGUI(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited, int index)
{
    // printf("updateLocallabGUI\n");

    // Update GUI values
    if (index < pp->locallab.nbspot && index < (int)pp->locallab.spots.size()) {
        // Color & Light
        expcolor->setEnabled(pp->locallab.spots.at(index).expcolor);
        curvactiv->set_active(pp->locallab.spots.at(index).curvactiv);
        lightness->setValue(pp->locallab.spots.at(index).lightness);
        contrast->setValue(pp->locallab.spots.at(index).contrast);
        chroma->setValue(pp->locallab.spots.at(index).chroma);
        sensi->setValue(pp->locallab.spots.at(index).sensi);

        if (pp->locallab.spots.at(index).qualitycurveMethod == "none") {
            qualitycurveMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).qualitycurveMethod == "std") {
            qualitycurveMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).qualitycurveMethod == "enh") {
            qualitycurveMethod->set_active(2);
        } else if (pp->locallab.spots.at(index).qualitycurveMethod == "enhsup") {
            qualitycurveMethod->set_active(3);
        } else if (pp->locallab.spots.at(index).qualitycurveMethod == "contr") {
            qualitycurveMethod->set_active(4);
        } else if (pp->locallab.spots.at(index).qualitycurveMethod == "sob2") {
            qualitycurveMethod->set_active(5);
        }

        if (pp->locallab.spots.at(index).showmaskcolMethod == "none") {
            showmaskcolMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).showmaskcolMethod == "color") {
            showmaskcolMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).showmaskcolMethod == "colormask") {
            showmaskcolMethod->set_active(2);
        } else if (pp->locallab.spots.at(index).showmaskcolMethod == "mask") {
            showmaskcolMethod->set_active(3);
        } else if (pp->locallab.spots.at(index).showmaskcolMethod == "showmask") {
            showmaskcolMethod->set_active(4);
        }
        llshape->setCurve(pp->locallab.spots.at(index).llcurve);
        ccshape->setCurve(pp->locallab.spots.at(index).cccurve);
        LHshape->setCurve(pp->locallab.spots.at(index).LHcurve);
        HHshape->setCurve(pp->locallab.spots.at(index).HHcurve);
        CCmaskshape->setCurve(pp->locallab.spots.at(index).CCmaskcurve);
        LLmaskshape->setCurve(pp->locallab.spots.at(index).LLmaskcurve);
        invers->set_active(pp->locallab.spots.at(index).invers);

        // Exposure
        expexpose->setEnabled(pp->locallab.spots.at(index).expexpose);
        expcomp->setValue(pp->locallab.spots.at(index).expcomp);
        hlcompr->setValue(pp->locallab.spots.at(index).hlcompr);
        hlcomprthresh->setValue(pp->locallab.spots.at(index).hlcomprthresh);
        black->setValue(pp->locallab.spots.at(index).black);
        shcompr->setValue(pp->locallab.spots.at(index).shcompr);
        warm->setValue(pp->locallab.spots.at(index).warm);
        sensiex->setValue(pp->locallab.spots.at(index).sensiex);
        shapeexpos->setCurve(pp->locallab.spots.at(index).excurve);
        if (pp->locallab.spots.at(index).showmaskexpMethod == "none") {
            showmaskexpMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).showmaskexpMethod == "expo") {
            showmaskexpMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).showmaskexpMethod == "expormask") {
            showmaskexpMethod->set_active(2);
        } else if (pp->locallab.spots.at(index).showmaskexpMethod == "mask") {
            showmaskexpMethod->set_active(3);
        } else if (pp->locallab.spots.at(index).showmaskexpMethod == "showmask") {
            showmaskexpMethod->set_active(4);
        }
        CCmaskexpshape->setCurve(pp->locallab.spots.at(index).CCmaskexpcurve);
        LLmaskexpshape->setCurve(pp->locallab.spots.at(index).LLmaskexpcurve);

        // Vibrance
        expvibrance->setEnabled(pp->locallab.spots.at(index).expvibrance);
        saturated->setValue(pp->locallab.spots.at(index).saturated);
        pastels->setValue(pp->locallab.spots.at(index).pastels);
        psThreshold->setValue<int>(pp->locallab.spots.at(index).psthreshold);
        protectSkins->set_active(pp->locallab.spots.at(index).protectskins);
        avoidColorShift->set_active(pp->locallab.spots.at(index).avoidcolorshift);
        pastSatTog->set_active(pp->locallab.spots.at(index).pastsattog);
        sensiv->setValue(pp->locallab.spots.at(index).sensiv);
        skinTonesCurve->setCurve(pp->locallab.spots.at(index).skintonescurve);

        // Soft Light
        expsoft->setEnabled(pp->locallab.spots.at(index).expsoft);
        streng->setValue(pp->locallab.spots.at(index).streng);
        sensisf->setValue(pp->locallab.spots.at(index).sensisf);

        // Lab Region
        explabregion->setEnabled(pp->locallab.spots.at(index).explabregion);

        // Blur & Noise
        expblur->setEnabled(pp->locallab.spots.at(index).expblur);
        radius->setValue(pp->locallab.spots.at(index).radius);
        strength->setValue(pp->locallab.spots.at(index).strength);
        sensibn->setValue(pp->locallab.spots.at(index).sensibn);

        if (pp->locallab.spots.at(index).blurMethod == "norm") {
            blurMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).blurMethod == "inv") {
            blurMethod->set_active(1);
        } else {
            blurMethod->set_active(2);
        }

        activlum->set_active(pp->locallab.spots.at(index).activlum);

        // Tone Mapping
        exptonemap->setEnabled(pp->locallab.spots.at(index).exptonemap);
        stren->setValue(pp->locallab.spots.at(index).stren);
        gamma->setValue(pp->locallab.spots.at(index).gamma);
        estop->setValue(pp->locallab.spots.at(index).estop);
        scaltm->setValue(pp->locallab.spots.at(index).scaltm);
        rewei->setValue(pp->locallab.spots.at(index).rewei);
        sensitm->setValue(pp->locallab.spots.at(index).sensitm);

        // Retinex
        expreti->setEnabled(pp->locallab.spots.at(index).expreti);

        if (pp->locallab.spots.at(index).retinexMethod == "low") {
            retinexMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).retinexMethod == "uni") {
            retinexMethod->set_active(1);
        } else {
            retinexMethod->set_active(2);
        }

        str->setValue(pp->locallab.spots.at(index).str);
        chrrt->setValue(pp->locallab.spots.at(index).chrrt);
        neigh->setValue(pp->locallab.spots.at(index).neigh);
        vart->setValue(pp->locallab.spots.at(index).vart);
        dehaz->setValue(pp->locallab.spots.at(index).dehaz);
        sensih->setValue(pp->locallab.spots.at(index).sensih);
        cTgainshape->setCurve(pp->locallab.spots.at(index).localTgaincurve);
        inversret->set_active(pp->locallab.spots.at(index).inversret);

        // Sharpening
        expsharp->setEnabled(pp->locallab.spots.at(index).expsharp);
        sharcontrast->setValue(pp->locallab.spots.at(index).sharcontrast);
        sharradius->setValue(pp->locallab.spots.at(index).sharradius);
        sharamount->setValue(pp->locallab.spots.at(index).sharamount);
        shardamping->setValue(pp->locallab.spots.at(index).shardamping);
        shariter->setValue(pp->locallab.spots.at(index).shariter);
        sharblur->setValue(pp->locallab.spots.at(index).sharblur);
        sensisha->setValue(pp->locallab.spots.at(index).sensisha);
        inverssha->set_active(pp->locallab.spots.at(index).inverssha);

        // Local Contrast
        expcontrast->setEnabled(pp->locallab.spots.at(index).expcontrast);
        lcradius->setValue(pp->locallab.spots.at(index).lcradius);
        lcamount->setValue(pp->locallab.spots.at(index).lcamount);
        lcdarkness->setValue(pp->locallab.spots.at(index).lcdarkness);
        lclightness->setValue(pp->locallab.spots.at(index).lclightness);
        sensilc->setValue(pp->locallab.spots.at(index).sensilc);

        // Contrast by detail levels
        expcbdl->setEnabled(pp->locallab.spots.at(index).expcbdl);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setValue(pp->locallab.spots.at(index).mult[i]);
        }

        chromacbdl->setValue(pp->locallab.spots.at(index).chromacbdl);
        threshold->setValue(pp->locallab.spots.at(index).threshold);
        sensicb->setValue(pp->locallab.spots.at(index).sensicb);

        // Denoise
        expdenoi->setEnabled(pp->locallab.spots.at(index).expdenoi);
        noiselumf->setValue(pp->locallab.spots.at(index).noiselumf);
        noiselumc->setValue(pp->locallab.spots.at(index).noiselumc);
        noiselumdetail->setValue(pp->locallab.spots.at(index).noiselumdetail);
        noiselequal->setValue(pp->locallab.spots.at(index).noiselequal);
        noisechrof->setValue(pp->locallab.spots.at(index).noisechrof);
        noisechroc->setValue(pp->locallab.spots.at(index).noisechroc);
        noisechrodetail->setValue(pp->locallab.spots.at(index).noisechrodetail);
        adjblur->setValue(pp->locallab.spots.at(index).adjblur);
        bilateral->setValue(pp->locallab.spots.at(index).bilateral);
        sensiden->setValue(pp->locallab.spots.at(index).sensiden);

        // Others
        avoid->set_active(pp->locallab.spots.at(index).avoid);

        if (pedited) {
            if (index < (int)pedited->locallab.spots.size()) {
                const LocallabParamsEdited::LocallabSpotEdited* spotState = &pedited->locallab.spots.at(index);

                // Control spot settings
                ControlSpotPanel::SpotEdited* const se = new ControlSpotPanel::SpotEdited();

                if (pedited->locallab.nbspot && pedited->locallab.id) {
                    se->nbspot = true;
                } else {
                    se->nbspot = false;
                }

                se->selspot = pedited->locallab.selspot;
                se->name = spotState->name;
                se->isvisible = spotState->isvisible;
                se->shape = spotState->shape;
                se->spotMethod = spotState->spotMethod;
                se->sensiexclu = spotState->sensiexclu;
                se->struc = spotState->struc;
                se->shapeMethod = spotState->shapeMethod;
                se->locX = spotState->locX;
                se->locXL = spotState->locXL;
                se->locY = spotState->locY;
                se->locYT = spotState->locYT;
                se->centerX = spotState->centerX;
                se->centerY = spotState->centerY;
                se->circrad = spotState->circrad;
                se->qualityMethod = spotState->qualityMethod;
                se->transit = spotState->transit;
                se->thresh = spotState->thresh;
                se->iter = spotState->iter;
                expsettings->setEditedStates(se);

                // Color & Light
                expcolor->set_inconsistent(!spotState->expcolor);
                curvactiv->set_inconsistent(multiImage && !spotState->curvactiv);
                lightness->setEditedState(spotState->lightness ? Edited : UnEdited);
                contrast->setEditedState(spotState->contrast ? Edited : UnEdited);
                chroma->setEditedState(spotState->chroma ? Edited : UnEdited);
                sensi->setEditedState(spotState->sensi ? Edited : UnEdited);

                if (!spotState->qualitycurveMethod) {
                    qualitycurveMethod->set_active_text(M("GENERAL_UNCHANGED"));
                }
                if (!spotState->showmaskcolMethod) {
                    showmaskcolMethod->set_active_text(M("GENERAL_UNCHANGED"));
                }

                llshape->setUnChanged(!spotState->llcurve);
                ccshape->setUnChanged(!spotState->cccurve);
                LHshape->setUnChanged(!spotState->LHcurve);
                HHshape->setUnChanged(!spotState->HHcurve);
                CCmaskshape->setUnChanged(!spotState->CCmaskcurve);
                LLmaskshape->setUnChanged(!spotState->LLmaskcurve);
                invers->set_inconsistent(multiImage && !spotState->invers);

                // Exposure
                expexpose->set_inconsistent(!spotState->expexpose);
                expcomp->setEditedState(spotState->expcomp ? Edited : UnEdited);
                hlcompr->setEditedState(spotState->hlcompr ? Edited : UnEdited);
                hlcomprthresh->setEditedState(spotState->hlcomprthresh ? Edited : UnEdited);
                black->setEditedState(spotState->black ? Edited : UnEdited);
                warm->setEditedState(spotState->warm ? Edited : UnEdited);
                shcompr->setEditedState(spotState->shcompr ? Edited : UnEdited);
                sensiex->setEditedState(spotState->sensiex ? Edited : UnEdited);
                shapeexpos->setUnChanged(!spotState->excurve);
                if (!spotState->showmaskexpMethod) {
                    showmaskexpMethod->set_active_text(M("GENERAL_UNCHANGED"));
                }
                CCmaskexpshape->setUnChanged(!spotState->CCmaskexpcurve);
                LLmaskexpshape->setUnChanged(!spotState->LLmaskexpcurve);

                // Vibrance
                expvibrance->set_inconsistent(!spotState->expvibrance);
                saturated->setEditedState(spotState->saturated ? Edited : UnEdited);
                pastels->setEditedState(spotState->pastels ? Edited : UnEdited);
                psThreshold->setEditedState(spotState->psthreshold ? Edited : UnEdited);
                protectSkins->set_inconsistent(multiImage && !spotState->protectskins);
                avoidColorShift->set_inconsistent(multiImage && !spotState->avoidcolorshift);
                pastSatTog->set_inconsistent(multiImage && !spotState->pastsattog);
                sensiv->setEditedState(spotState->sensiv ? Edited : UnEdited);
                skinTonesCurve->setUnChanged(!spotState->skintonescurve);

                // Soft Light
                expsoft->set_inconsistent(!spotState->expsoft);
                streng->setEditedState(spotState->streng ? Edited : UnEdited);
                sensisf->setEditedState(spotState->sensisf ? Edited : UnEdited);

                // Lab Region
                explabregion->set_inconsistent(!spotState->explabregion);

                // Blur & Noise
                expblur->set_inconsistent(!spotState->expblur);
                radius->setEditedState(spotState->radius ? Edited : UnEdited);
                strength->setEditedState(spotState->strength ? Edited : UnEdited);
                sensibn->setEditedState(spotState->sensibn ? Edited : UnEdited);

                if (!spotState->blurMethod) {
                    blurMethod->set_active_text(M("GENERAL_UNCHANGED"));
                }

                activlum->set_inconsistent(multiImage && !spotState->activlum);

                // Tone Mapping
                exptonemap->set_inconsistent(!spotState->exptonemap);
                stren->setEditedState(spotState->stren ? Edited : UnEdited);
                gamma->setEditedState(spotState->gamma ? Edited : UnEdited);
                estop->setEditedState(spotState->estop ? Edited : UnEdited);
                scaltm->setEditedState(spotState->scaltm ? Edited : UnEdited);
                rewei->setEditedState(spotState->rewei ? Edited : UnEdited);
                sensitm->setEditedState(spotState->sensitm ? Edited : UnEdited);

                // Retinex
                expreti->set_inconsistent(!spotState->expreti);

                if (!spotState->retinexMethod) {
                    retinexMethod->set_active_text(M("GENERAL_UNCHANGED"));
                }

                str->setEditedState(spotState->str ? Edited : UnEdited);
                chrrt->setEditedState(spotState->chrrt ? Edited : UnEdited);
                neigh->setEditedState(spotState->neigh ? Edited : UnEdited);
                vart->setEditedState(spotState->vart ? Edited : UnEdited);
                dehaz->setEditedState(spotState->dehaz ? Edited : UnEdited);
                sensih->setEditedState(spotState->sensih ? Edited : UnEdited);
                cTgainshape->setUnChanged(!spotState->localTgaincurve);
                inversret->set_inconsistent(multiImage && !spotState->inversret);

                // Sharpening
                expsharp->set_inconsistent(!spotState->expsharp);
                sharcontrast->setEditedState(spotState->sharcontrast ? Edited : UnEdited);
                sharradius->setEditedState(spotState->sharradius ? Edited : UnEdited);
                sharamount->setEditedState(spotState->sharamount ? Edited : UnEdited);
                shardamping->setEditedState(spotState->shardamping ? Edited : UnEdited);
                shariter->setEditedState(spotState->shariter ? Edited : UnEdited);
                sharblur->setEditedState(spotState->sharblur ? Edited : UnEdited);
                sensisha->setEditedState(spotState->sensisha ? Edited : UnEdited);
                inverssha->set_inconsistent(multiImage && !spotState->inverssha);

                // Local Contrast
                expcontrast->set_inconsistent(!spotState->expcontrast);
                lcradius->setEditedState(spotState->lcradius ? Edited : UnEdited);
                lcamount->setEditedState(spotState->lcamount ? Edited : UnEdited);
                lcdarkness->setEditedState(spotState->lcdarkness ? Edited : UnEdited);
                lclightness->setEditedState(spotState->lclightness ? Edited : UnEdited);
                sensilc->setEditedState(spotState->sensilc ? Edited : UnEdited);

                // Contrast by detail levels
                expcbdl->set_inconsistent(!spotState->expcbdl);

                for (int i = 0; i < 5; i++) {
                    multiplier[i]->setEditedState(spotState->mult[i] ? Edited : UnEdited);
                }

                chromacbdl->setEditedState(spotState->chromacbdl ? Edited : UnEdited);
                threshold->setEditedState(spotState->threshold ? Edited : UnEdited);
                sensicb->setEditedState(spotState->sensicb ? Edited : UnEdited);

                // Denoise
                expdenoi->set_inconsistent(!spotState->expdenoi);
                noiselumf->setEditedState(spotState->noiselumf ? Edited : UnEdited);
                noiselumc->setEditedState(spotState->noiselumc ? Edited : UnEdited);
                noiselumdetail->setEditedState(spotState->noiselumdetail ? Edited : UnEdited);
                noiselequal->setEditedState(spotState->noiselequal ? Edited : UnEdited);
                noisechrof->setEditedState(spotState->noisechrof ? Edited : UnEdited);
                noisechroc->setEditedState(spotState->noisechroc ? Edited : UnEdited);
                noisechrodetail->setEditedState(spotState->noisechrodetail ? Edited : UnEdited);
                adjblur->setEditedState(spotState->adjblur ? Edited : UnEdited);
                bilateral->setEditedState(spotState->bilateral ? Edited : UnEdited);
                sensiden->setEditedState(spotState->sensiden ? Edited : UnEdited);

                // Others
                avoid->set_inconsistent(multiImage && !spotState->avoid);
            }
        }
    }
}

void Locallab::updateSpecificGUIState()
{
    // Update Color & Light GUI according to invers button state (to be compliant with inversChanged function)
    if (multiImage && invers->get_inconsistent()) {
        sensi->show();
        llCurveEditorG->show();
        curvactiv->show();
        qualitycurveMethod->show();
        labqualcurv->show();
        maskcolFrame->show();
    } else if (invers->get_active()) {
        sensi->show();
        llCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        maskcolFrame->hide();
        labqualcurv->hide();
    } else {
        sensi->show();
        llCurveEditorG->show();
        curvactiv->show();
        qualitycurveMethod->show();
        labqualcurv->show();
        maskcolFrame->show();
    }

    // Update Exposure GUI according to black adjuster state (to be compliant with adjusterChanged function)
    if (multiImage && black->getEditedState() != Edited) {
        shcompr->set_sensitive(true);
    } else {
        shcompr->set_sensitive(!((int)black->getValue() == 0)); // At black = 0, shcompr value has no effect
    }

    // Update Vibrance GUI according to pastsattog button state (to be compliant with pastsattog_toggled function)
    if (multiImage && pastSatTog->get_inconsistent()) {
        psThreshold->set_sensitive(true);
        saturated->set_sensitive(true);
    } else if (pastSatTog->get_active()) {
        // Link both slider, so we set saturated and psThresholds unsensitive
        psThreshold->set_sensitive(false);
        saturated->set_sensitive(false);
        saturated->setValue(pastels->getValue()); // Pastels and Saturated are linked
    } else {
        // Separate sliders, so we set saturated and psThresholds sensitive again
        psThreshold->set_sensitive(true);
        saturated->set_sensitive(true);
    }

    // Update Blur & Noise GUI according to blurMethod combobox (to be compliant with blurMethodChanged function)
    if (multiImage && blurMethod->get_active_text() == M("GENERAL_UNCHANGED")) {
        sensibn->show();
    } else if (blurMethod->get_active_row_number() == 0 || blurMethod->get_active_row_number() == 2) {
        sensibn->show();
    } else {
        sensibn->hide();
    }

    if (blurMethod->get_active_row_number() == 2) {
        strength->hide();
    } else {
        strength->show();
    }
    
    //update showmethod
    if (multiImage && showmaskcolMethod->get_active_text() == M("GENERAL_UNCHANGED")) {
        showmaskexpMethod->set_active(0);
    } else if((showmaskcolMethod->get_active_row_number() == 1 || showmaskcolMethod->get_active_row_number() == 2 || showmaskcolMethod->get_active_row_number() == 4)) {
        showmaskexpMethod->set_active(0);
        expexpose->setEnabled(false); 
    }
    
    if (multiImage && showmaskexpMethod->get_active_text() == M("GENERAL_UNCHANGED")) {
        showmaskcolMethod->set_active(0);
    } else if((showmaskexpMethod->get_active_row_number() == 1 || showmaskexpMethod->get_active_row_number() == 2 || showmaskexpMethod->get_active_row_number() == 4)) {
        showmaskcolMethod->set_active(0);
        expcolor->setEnabled(false); 
    }
   

    // Update Retinex GUI according to inversret button state (to be compliant with inversretChanged function)
    if (multiImage && inversret->get_inconsistent()) {
        sensih->show();
    } else if (inversret->get_active()) {
        sensih->hide();
    } else {
        sensih->show();
    }

    // Update Sharpening GUI according to inverssha button state (to be compliant with inversshaChanged function)
    if (multiImage && inverssha->get_inconsistent()) {
        sensisha->show();
    } else if (inverssha->get_active()) {
        sensisha->hide();
    } else {
        sensisha->show();
    }
}

void Locallab::autoOpenCurve()
{
    // printf("autoOpenCurve\n");

    // TODO autoOpenCurve only considers linearity state of selected spot curve
    llshape->openIfNonlinear();
    ccshape->openIfNonlinear();
    LHshape->openIfNonlinear();
    HHshape->openIfNonlinear();
    CCmaskshape->openIfNonlinear();
    LLmaskshape->openIfNonlinear();
    CCmaskexpshape->openIfNonlinear();
    LLmaskexpshape->openIfNonlinear();
}
