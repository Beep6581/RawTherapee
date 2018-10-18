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
    expblur(new MyExpander(true, M("TP_LOCALLAB_BLUFR"))),
    exptonemap(new MyExpander(true, M("TP_LOCALLAB_TM"))),
    expreti(new MyExpander(true, M("TP_LOCALLAB_RETI"))),
    expsharp(new MyExpander(true, M("TP_LOCALLAB_SHARP"))),
    expcbdl(new MyExpander(true, M("TP_LOCALLAB_CBDL"))),
    expdenoi(new MyExpander(true, M("TP_LOCALLAB_DENOIS"))),

    // CurveEditorGroup widgets
    // Color & Light
    llCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_LUM"))),
    // Exposure
    curveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVEEDITOR_TONES_LABEL"))),
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
    hlcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 20))),
    hlcomprthresh(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 33))),
    black(Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0))),
    shcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50))),
    warm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WARM"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-orange-small.png"))))),
    sensiex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    // Vibrance
    saturated(Gtk::manage(new Adjuster(M("TP_VIBRANCE_SATURATED"), -100., 100., 1., 0.))),
    pastels(Gtk::manage(new Adjuster(M("TP_VIBRANCE_PASTELS"), -100., 100., 1., 0.))),
    sensiv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
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
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIH"), 0, 100, 1, 19))),
    // Sharpening
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 42, 500, 1, 4))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 75))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 75))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
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
    lumacontrastPlusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS"))))
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

    curvactivConn = curvactiv->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::curvactivChanged));

    lightness->setAdjusterListener(this);

    contrast->setAdjusterListener(this);

    chroma->setAdjusterListener(this);

    sensi->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensi->setAdjusterListener(this);

    qualitycurveMethod->append(M("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVENH"));
    qualitycurveMethod->set_active(0);
    qualitycurveMethod->set_tooltip_markup(M("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));
    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::qualitycurveMethodChanged));

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

    inversConn  = invers->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversChanged));

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

    ToolParamBlock* const exposeBox = Gtk::manage(new ToolParamBlock());
    exposeBox->pack_start(*expcomp);
    exposeBox->pack_start(*hlcompr);
    exposeBox->pack_start(*hlcomprthresh);
    exposeBox->pack_start(*black);
    exposeBox->pack_start(*shcompr);
    exposeBox->pack_start(*warm);
    exposeBox->pack_start(*sensiex);
    exposeBox->pack_start(*curveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
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

    panel->pack_start(*exptonemap, false, false);

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
    retiBox->pack_start(*sensih);
    retiBox->pack_start(*LocalcurveEditorgainT, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    retiBox->pack_start(*inversret);
    expreti->add(*retiBox);
    expreti->setLevel(2);

    panel->pack_start(*expreti, false, false);

    // Sharpening
    expsharp->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsharp));
    enablesharpConn = expsharp->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expsharp));

    sharradius->setAdjusterListener(this);

    sharamount->setAdjusterListener(this);

    shardamping->setAdjusterListener(this);

    shariter->setAdjusterListener(this);

    sensisha->set_tooltip_text(M("TP_LOCALLAB_SENSIS_TOOLTIP"));
    sensisha->setAdjusterListener(this);

    inversshaConn  = inverssha->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversshaChanged));

    ToolParamBlock* const sharpBox = Gtk::manage(new ToolParamBlock());
    sharpBox->pack_start(*sharradius);
    sharpBox->pack_start(*sharamount);
    sharpBox->pack_start(*shardamping);
    sharpBox->pack_start(*shariter);
    sharpBox->pack_start(*sensisha);
    sharpBox->pack_start(*inverssha);
    expsharp->add(*sharpBox);
    expsharp->setLevel(2);

    panel->pack_start(*expsharp, false, false);

    // Contrast by detail levels
    expcbdl->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcbdl));
    enablecbdlConn = expcbdl->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcbdl));

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
        expblur->set_expanded(expblur == expander);
        exptonemap->set_expanded(exptonemap == expander);
        expreti->set_expanded(expreti == expander);
        expsharp->set_expanded(expsharp == expander);
        expcbdl->set_expanded(expcbdl == expander);
        expdenoi->set_expanded(expdenoi == expander);

    }
}

void Locallab::enableToggled(MyExpander *expander)
{
    // TODO Locallab printf
    printf("enableToggled\n");

    if (getEnabled() && listener) {
        rtengine::ProcEvent event = NUMOFEVENTS;

        if (expander == expcolor) {
            event = EvLocenacolor;
        } else if (expander == expexpose) {
            event = EvLocenaexpose;
        } else if (expander == expvibrance) {
            event = EvLocenavibrance;
        } else if (expander == expblur) {
            event = EvLocenablur;
        } else if (expander == exptonemap) {
            event = EvLocenatonemap;
        } else if (expander == expreti) {
            event = EvLocenareti;
        } else if (expander == expsharp) {
            event = EvLocenasharp;
        } else if (expander == expcbdl) {
            event = EvLocenacbdl;
        } else if (expander == expdenoi) {
            event = EvLocenadenoi;
        } else {
            return;
        }

        if (expander->get_inconsistent()) {
            listener->panelChanged(event, M("GENERAL_UNCHANGED"));
        } else if (expander->getEnabled()) {
            listener->panelChanged(event, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(event, M("GENERAL_DISABLED"));
        }

    }
}

void Locallab::writeOptions(std::vector<int> &tpOpen)
{
    tpOpen.push_back(expsettings->getExpanded());
    tpOpen.push_back(expcolor->get_expanded());
    tpOpen.push_back(expexpose->get_expanded());
    tpOpen.push_back(expvibrance->get_expanded());
    tpOpen.push_back(expblur->get_expanded());
    tpOpen.push_back(exptonemap->get_expanded());
    tpOpen.push_back(expreti->get_expanded());
    tpOpen.push_back(expsharp->get_expanded());
    tpOpen.push_back(expcbdl->get_expanded());
    tpOpen.push_back(expdenoi->get_expanded());

}

void Locallab::updateToolState(std::vector<int> &tpOpen)
{
    if (tpOpen.size() >= 10) {
        expsettings->setExpanded(tpOpen.at(0));
        expcolor->set_expanded(tpOpen.at(1));
        expexpose->set_expanded(tpOpen.at(2));
        expvibrance->set_expanded(tpOpen.at(3));
        expblur->set_expanded(tpOpen.at(4));
        exptonemap->set_expanded(tpOpen.at(5));
        expreti->set_expanded(tpOpen.at(6));
        expsharp->set_expanded(tpOpen.at(7));
        expcbdl->set_expanded(tpOpen.at(8));
        expdenoi->set_expanded(tpOpen.at(9));
    }
}

void Locallab::lumaneutralPressed()
{
    // TODO Locallab printf
    printf("lumaneutralPressed\n");

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue(100);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void Locallab::lumacontrastPlusPressed()
{
    // TODO Locallab printf
    printf("lumacontrastPlusPressed\n");

    for (int i = 0; i < 5; i++) {
        float inc = (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void Locallab::lumacontrastMinusPressed()
{
    // TODO Locallab printf
    printf("lumacontrastMinusPressed\n");

    for (int i = 0; i < 5; i++) {
        float inc = - (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void Locallab::read(const ProcParams* pp, const ParamsEdited* pedited)
{
    printf("Locallab read\n");

    // Disable all listeners
    disableListener();

    if (pedited) {
        set_inconsistent(multiImage && !pedited->locallab.enabled);

        // Control spot settings
        ControlSpotPanel::SpotEdited* const se = new ControlSpotPanel::SpotEdited();

        if (pedited->locallab.nbspot) {
            se->addbutton = true;
            se->deletebutton = true;
        } else {
            se->addbutton = false;
            se->deletebutton = false;
        }

        se->treeview = pedited->locallab.nbspot || pedited->locallab.selspot;
        se->name = pedited->locallab.name;
        se->isvisible = pedited->locallab.isvisible;
        se->shape = pedited->locallab.shape;
        se->spotMethod = pedited->locallab.spotMethod;
        se->sensiexclu = pedited->locallab.sensiexclu;
        se->struc = pedited->locallab.struc;
        se->shapeMethod = pedited->locallab.shapeMethod;
        se->locX = pedited->locallab.locX;
        se->locXL = pedited->locallab.locXL;
        se->locY = pedited->locallab.locY;
        se->locYT = pedited->locallab.locYT;
        se->centerX = pedited->locallab.centerX;
        se->centerY = pedited->locallab.centerY;
        se->circrad = pedited->locallab.circrad;
        se->qualityMethod = pedited->locallab.qualityMethod;
        se->transit = pedited->locallab.transit;
        se->thresh = pedited->locallab.thresh;
        se->iter = pedited->locallab.iter;
        expsettings->setEditedStates(se);

        // Color & Light
        expcolor->set_inconsistent(!pedited->locallab.expcolor);
        curvactiv->set_inconsistent(multiImage && !pedited->locallab.curvactiv);
        lightness->setEditedState(pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setEditedState(pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setEditedState(pedited->locallab.chroma ? Edited : UnEdited);
        sensi->setEditedState(pedited->locallab.sensi ? Edited : UnEdited);

        if (!pedited->locallab.qualitycurveMethod) {
            qualitycurveMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        llshape->setUnChanged(!pedited->locallab.llcurve);
        ccshape->setUnChanged(!pedited->locallab.cccurve);
        LHshape->setUnChanged(!pedited->locallab.LHcurve);
        HHshape->setUnChanged(!pedited->locallab.HHcurve);
        invers->set_inconsistent(multiImage && !pedited->locallab.invers);

        // Exposure
        expexpose->set_inconsistent(!pedited->locallab.expexpose);
        expcomp->setEditedState(pedited->locallab.expcomp ? Edited : UnEdited);
        hlcompr->setEditedState(pedited->locallab.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setEditedState(pedited->locallab.hlcomprthresh ? Edited : UnEdited);
        black->setEditedState(pedited->locallab.black ? Edited : UnEdited);
        warm->setEditedState(pedited->locallab.warm ? Edited : UnEdited);
        shcompr->setEditedState(pedited->locallab.shcompr ? Edited : UnEdited);
        sensiex->setEditedState(pedited->locallab.sensiex ? Edited : UnEdited);
        shapeexpos->setUnChanged(!pedited->locallab.excurve);

        // Vibrance
        expvibrance->set_inconsistent(!pedited->locallab.expvibrance);
        saturated->setEditedState(pedited->locallab.saturated ? Edited : UnEdited);
        pastels->setEditedState(pedited->locallab.pastels ? Edited : UnEdited);
        psThreshold->setEditedState(pedited->locallab.psthreshold ? Edited : UnEdited);
        protectSkins->set_inconsistent(!pedited->locallab.protectskins);
        avoidColorShift->set_inconsistent(!pedited->locallab.avoidcolorshift);
        pastSatTog->set_inconsistent(!pedited->locallab.pastsattog);
        sensiv->setEditedState(pedited->locallab.sensiv ? Edited : UnEdited);
        skinTonesCurve->setUnChanged(!pedited->locallab.skintonescurve);

        // Blur & Noise
        expblur->set_inconsistent(!pedited->locallab.expblur);
        radius->setEditedState(pedited->locallab.radius ? Edited : UnEdited);
        strength->setEditedState(pedited->locallab.strength ? Edited : UnEdited);
        sensibn->setEditedState(pedited->locallab.sensibn ? Edited : UnEdited);

        if (!pedited->locallab.blurMethod) {
            blurMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        activlum->set_inconsistent(multiImage && !pedited->locallab.activlum);

        // Tone Mapping
        exptonemap->set_inconsistent(!pedited->locallab.exptonemap);
        stren->setEditedState(pedited->locallab.stren ? Edited : UnEdited);
        gamma->setEditedState(pedited->locallab.gamma ? Edited : UnEdited);
        estop->setEditedState(pedited->locallab.estop ? Edited : UnEdited);
        scaltm->setEditedState(pedited->locallab.scaltm ? Edited : UnEdited);
        rewei->setEditedState(pedited->locallab.rewei ? Edited : UnEdited);
        sensitm->setEditedState(pedited->locallab.sensitm ? Edited : UnEdited);

        // Retinex
        expreti->set_inconsistent(!pedited->locallab.expreti);

        if (!pedited->locallab.retinexMethod) {
            retinexMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        str->setEditedState(pedited->locallab.str ? Edited : UnEdited);
        chrrt->setEditedState(pedited->locallab.chrrt ? Edited : UnEdited);
        neigh->setEditedState(pedited->locallab.neigh ? Edited : UnEdited);
        vart->setEditedState(pedited->locallab.vart ? Edited : UnEdited);
        sensih->setEditedState(pedited->locallab.sensih ? Edited : UnEdited);
        cTgainshape->setUnChanged(!pedited->locallab.localTgaincurve);
        inversret->set_inconsistent(multiImage && !pedited->locallab.inversret);

        // Sharpening
        expsharp->set_inconsistent(!pedited->locallab.expsharp);
        sharradius->setEditedState(pedited->locallab.sharradius ? Edited : UnEdited);
        sharamount->setEditedState(pedited->locallab.sharamount ? Edited : UnEdited);
        shardamping->setEditedState(pedited->locallab.shardamping ? Edited : UnEdited);
        shariter->setEditedState(pedited->locallab.shariter ? Edited : UnEdited);
        sensisha->setEditedState(pedited->locallab.sensisha ? Edited : UnEdited);
        inverssha->set_inconsistent(multiImage && !pedited->locallab.inverssha);

        // Contrast by detail levels
        expcbdl->set_inconsistent(!pedited->locallab.expcbdl);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setEditedState(pedited->locallab.mult[i] ? Edited : UnEdited);
        }

        chromacbdl->setEditedState(pedited->locallab.chromacbdl ? Edited : UnEdited);
        threshold->setEditedState(pedited->locallab.threshold ? Edited : UnEdited);
        sensicb->setEditedState(pedited->locallab.sensicb ? Edited : UnEdited);

        // Denoise
        expdenoi->set_inconsistent(!pedited->locallab.expdenoi);
        noiselumf->setEditedState(pedited->locallab.noiselumf ? Edited : UnEdited);
        noiselumc->setEditedState(pedited->locallab.noiselumc ? Edited : UnEdited);
        noiselumdetail->setEditedState(pedited->locallab.noiselumdetail ? Edited : UnEdited);
        noiselequal->setEditedState(pedited->locallab.noiselequal ? Edited : UnEdited);
        noisechrof->setEditedState(pedited->locallab.noisechrof ? Edited : UnEdited);
        noisechroc->setEditedState(pedited->locallab.noisechroc ? Edited : UnEdited);
        noisechrodetail->setEditedState(pedited->locallab.noisechrodetail ? Edited : UnEdited);
        adjblur->setEditedState(pedited->locallab.adjblur ? Edited : UnEdited);
        bilateral->setEditedState(pedited->locallab.bilateral ? Edited : UnEdited);
        sensiden->setEditedState(pedited->locallab.sensiden ? Edited : UnEdited);

        // Others
        avoid->set_inconsistent(multiImage && !pedited->locallab.avoid);
    }

    setEnabled(pp->locallab.enabled);

    // Add non existent spots and update existent ones
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

        if (!expsettings->updateControlSpot(r)) {
            expsettings->addControlSpot(r);
        }
    }

    // Delete not anymore existent spots
    std::vector<int>* const list = expsettings->getSpotIdList();
    bool ispresent;

    for (int i = 0; i < (int)list->size(); i++) {
        ispresent = false;

        for (int j = 0; j < pp->locallab.nbspot && j < (int)pp->locallab.spots.size(); j++) {
            if (list->at(i) == pp->locallab.spots.at(j).id) {
                ispresent = true;
                break;
            }
        }

        if (!ispresent) {
            expsettings->deleteControlSpot(list->at(i));
        }
    }

    // Select active spot
    if (pp->locallab.nbspot > 0) {
        expsettings->setSelectedSpot(pp->locallab.spots.at(pp->locallab.selspot).id);
    }

    // Update Locallab tools GUI
    updateLocallabGUI(pp, pp->locallab.selspot);

    // Enable all listeners
    enableListener();
}

void Locallab::write(ProcParams* pp, ParamsEdited* pedited)
{
    printf("Locallab write\n");

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
            r->name = newSpot->name = "Control Spot #" + std::to_string(spotId);
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
            updateLocallabGUI(pp, pp->locallab.selspot);
            enableListener();

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
                    updateLocallabGUI(pp, pp->locallab.selspot);
                    enableListener();

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
            updateLocallabGUI(pp, pp->locallab.selspot);
            enableListener();

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

                    pp->locallab.spots.at(pp->locallab.selspot).llcurve = llshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).cccurve = ccshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).LHcurve = LHshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).HHcurve = HHshape->getCurve();
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
                    pp->locallab.spots.at(pp->locallab.selspot).sensih = sensih->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).localTgaincurve = cTgainshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).inversret = inversret->get_active();
                    // Sharpening
                    pp->locallab.spots.at(pp->locallab.selspot).expsharp = expsharp->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).sharradius = sharradius->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sharamount = sharamount->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shardamping = shardamping->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shariter = shariter->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensisha = sensisha->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).inverssha = inverssha->get_active();
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
            }

            // Update Locallab tools GUI
            disableListener();
            updateSpecificGUIState();
            enableListener();
    }

    if (pedited) {
        pedited->locallab.enabled = !get_inconsistent();
        // Control spot settings
        ControlSpotPanel::SpotEdited* const se = expsettings->getEditedStates();
        pedited->locallab.nbspot = se->addbutton || se->deletebutton;
        pedited->locallab.selspot = se->treeview;
        pedited->locallab.id = se->addbutton || se->deletebutton;
        pedited->locallab.name = se->name;
        pedited->locallab.isvisible = se->isvisible;
        pedited->locallab.shape = se->shape;
        pedited->locallab.spotMethod = se->spotMethod;
        pedited->locallab.sensiexclu = se->sensiexclu;
        pedited->locallab.struc = se->struc;
        pedited->locallab.shapeMethod = se->shapeMethod;
        pedited->locallab.locX = se->locX;
        pedited->locallab.locXL = se->locXL;
        pedited->locallab.locY = se->locY;
        pedited->locallab.locYT = se->locYT;
        pedited->locallab.centerX = se->centerX;
        pedited->locallab.centerY = se->centerY;
        pedited->locallab.circrad = se->circrad;
        pedited->locallab.qualityMethod = se->qualityMethod;
        pedited->locallab.transit = se->transit;
        pedited->locallab.thresh = se->thresh;
        pedited->locallab.iter = se->iter;
        // Color & Light
        pedited->locallab.expcolor = !expcolor->get_inconsistent();
        pedited->locallab.curvactiv = !curvactiv->get_inconsistent();
        pedited->locallab.lightness = lightness->getEditedState();
        pedited->locallab.contrast = contrast->getEditedState();
        pedited->locallab.chroma = chroma->getEditedState();
        pedited->locallab.sensi = sensi->getEditedState();
        pedited->locallab.qualitycurveMethod = qualitycurveMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.llcurve = !llshape->isUnChanged();
        pedited->locallab.cccurve = !ccshape->isUnChanged();
        pedited->locallab.LHcurve = !LHshape->isUnChanged();
        pedited->locallab.HHcurve = !HHshape->isUnChanged();
        pedited->locallab.invers = !invers->get_inconsistent();
        // Exposure
        pedited->locallab.expexpose = !expexpose->get_inconsistent();
        pedited->locallab.expcomp = expcomp->getEditedState();
        pedited->locallab.hlcompr = hlcompr->getEditedState();
        pedited->locallab.hlcomprthresh = hlcomprthresh->getEditedState();
        pedited->locallab.black = black->getEditedState();
        pedited->locallab.shcompr = shcompr->getEditedState();
        pedited->locallab.warm = warm->getEditedState();
        pedited->locallab.sensiex = sensiex->getEditedState();
        pedited->locallab.excurve = !shapeexpos->isUnChanged();
        // Vibrance
        pedited->locallab.expvibrance = !expvibrance->get_inconsistent();
        pedited->locallab.saturated = saturated->getEditedState();
        pedited->locallab.pastels = pastels->getEditedState();
        pedited->locallab.psthreshold = psThreshold->getEditedState();
        pedited->locallab.protectskins = !protectSkins->get_inconsistent();
        pedited->locallab.avoidcolorshift = !avoidColorShift->get_inconsistent();
        pedited->locallab.pastsattog = !pastSatTog->get_inconsistent();
        pedited->locallab.sensiv = sensiv->getEditedState();
        pedited->locallab.skintonescurve = !skinTonesCurve->isUnChanged();
        // Blur & Noise
        pedited->locallab.expblur = !expblur->get_inconsistent();
        pedited->locallab.radius = radius->getEditedState();
        pedited->locallab.strength = strength->getEditedState();
        pedited->locallab.sensibn = sensibn->getEditedState();
        pedited->locallab.blurMethod = blurMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.activlum = !activlum->get_inconsistent();
        // Tone Mapping
        pedited->locallab.exptonemap = !exptonemap->get_inconsistent();
        pedited->locallab.stren = stren->getEditedState();
        pedited->locallab.gamma = gamma->getEditedState();
        pedited->locallab.estop = estop->getEditedState();
        pedited->locallab.scaltm = scaltm->getEditedState();
        pedited->locallab.rewei = rewei->getEditedState();
        pedited->locallab.sensitm = sensitm->getEditedState();
        // Retinex
        pedited->locallab.expreti = !expreti->get_inconsistent();
        pedited->locallab.retinexMethod = retinexMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.str = str->getEditedState();
        pedited->locallab.chrrt = chrrt->getEditedState();
        pedited->locallab.neigh = neigh->getEditedState();
        pedited->locallab.vart = vart->getEditedState();
        pedited->locallab.sensih = sensih->getEditedState();
        pedited->locallab.localTgaincurve = !cTgainshape->isUnChanged();
        pedited->locallab.inversret = !inversret->get_inconsistent();
        // Sharpening
        pedited->locallab.expsharp = !expsharp->get_inconsistent();
        pedited->locallab.sharradius = sharradius->getEditedState();
        pedited->locallab.sharamount = sharamount->getEditedState();
        pedited->locallab.shardamping = shardamping->getEditedState();
        pedited->locallab.shariter = shariter->getEditedState();
        pedited->locallab.sensisha = sensisha->getEditedState();
        pedited->locallab.inverssha = !inverssha->get_inconsistent();
        // Contrast by detail levels
        pedited->locallab.expcbdl = !expcbdl->get_inconsistent();

        for (int i = 0; i < 5; i++) {
            pedited->locallab.mult[i] = multiplier[i]->getEditedState();
        }

        pedited->locallab.chromacbdl = chromacbdl->getEditedState();
        pedited->locallab.threshold = threshold->getEditedState();
        pedited->locallab.sensicb = sensicb->getEditedState();
        // Denoise
        pedited->locallab.expdenoi = !expdenoi->get_inconsistent();
        pedited->locallab.noiselumf = noiselumf->getEditedState();
        pedited->locallab.noiselumc = noiselumc->getEditedState();
        pedited->locallab.noiselumdetail = noiselumdetail->getEditedState();
        pedited->locallab.noiselequal = noiselequal->getEditedState();
        pedited->locallab.noisechrof = noisechrof->getEditedState();
        pedited->locallab.noisechroc = noisechroc->getEditedState();
        pedited->locallab.noisechrodetail = noisechrodetail->getEditedState();
        pedited->locallab.adjblur = adjblur->getEditedState();
        pedited->locallab.bilateral = bilateral->getEditedState();
        pedited->locallab.sensiden = sensiden->getEditedState();
        // Others
        pedited->locallab.avoid = !avoid->get_inconsistent();
    }
}

void Locallab::protectskins_toggled()
{
    printf("protectskins_toggled\n");

    if (batchMode) {
        /*
        if (protectSkins->get_inconsistent()) {
            protectSkins->set_inconsistent(false);
            pskinsconn.block(true);
            protectSkins->set_active(false);
            pskinsconn.block(false);
        } else if (lastProtectSkins) {
            protectSkins->set_inconsistent(true);
        }

        lastProtectSkins = protectSkins->get_active();
        */
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
    printf("avoidcolorshift_toggled\n");

    if (batchMode) {
        /*
        if (avoidColorShift->get_inconsistent()) {
            avoidColorShift->set_inconsistent(false);
            ashiftconn.block(true);
            avoidColorShift->set_active(false);
            ashiftconn.block(false);
        } else if (lastAvoidColorShift) {
            avoidColorShift->set_inconsistent(true);
        }

        lastAvoidColorShift = avoidColorShift->get_active();
        */
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
    printf("pastsattog_toggled\n");

    if (batchMode) {
        /*
        if (pastSatTog->get_inconsistent()) {
            pastSatTog->set_inconsistent(false);
            pastsattogconn.block(true);
            pastSatTog->set_active(false);
            pastsattogconn.block(false);
        } else if (lastPastSatTog) {
            pastSatTog->set_inconsistent(true);
        }

        lastPastSatTog = pastSatTog->get_active();
        */
    }

    // Update Vibrance GUI according to pastsattog button state (to be compliant with updateSpecificGUIState function)
    if (pastSatTog->get_active()) {
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
    }

    // Exposure
    if (getEnabled() && expexpose->getEnabled()) {
        if (ce == shapeexpos) {
            if (listener) {
                listener->panelChanged(Evlocallabshapeexpos, M("HISTORY_CUSTOMCURVE"));
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

    /*
    if (ce == skinTonesCurve) {
        listener->panelChanged(EvlocallabSkinTonesCurve, M("HISTORY_CUSTOMCURVE"));
        int strval = retrab->getValue();
        //update MIP
        retrab->setValue(strval + 1);
        adjusterChanged(retrab, strval + 1);
        usleep(10000);  //to test
        retrab->setValue(strval);
        adjusterChanged(retrab, strval);
    }
    */
}

void Locallab::retinexMethodChanged()
{
    printf("retinexMethodChanged\n");

    if (!batchMode) {
        /*
        retrab->hide();
        LocalcurveEditorgainTrab->hide();
        */
    }

    if (getEnabled() && expreti->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabretinexMethod, retinexMethod->get_active_text());
        }
    }
}

void Locallab::blurMethodChanged()
{
    printf("blurMethodChanged\n");

    // Update Blur & Noise GUI according to blurMethod combobox (to be compliant with updateSpecificGUIState function)
    if (blurMethod->get_active_row_number() == 0 || blurMethod->get_active_row_number() == 2) {
        sensibn->show();
    } else {
        sensibn->hide();
    }

    if (getEnabled() && expblur->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabblurMethod, blurMethod->get_active_text());
        }
    }
}

void Locallab::qualitycurveMethodChanged()
{
    printf("qualitycurveMethodChanged\n");

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text());
        }
    }
}

void Locallab::inversChanged()
{
    printf("inversChanged\n");

    /*
    if (batchMode) {
        if (invers->get_inconsistent()) {
            invers->set_inconsistent(false);
            inversConn.block(true);
            invers->set_active(false);
            inversConn.block(false);
        } else if (lastinvers) {
            invers->set_inconsistent(true);
        }

        lastinvers = invers->get_active();
    }
    */

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
    printf("curvactivChanged\n");

    // TODO Batch mode
    /*
    if (batchMode) {
        if (curvactiv->get_inconsistent()) {
            curvactiv->set_inconsistent(false);
            curvactivConn.block(true);
            curvactiv->set_active(false);
            curvactivConn.block(false);
        } else if (lastcurvactiv) {
            curvactiv->set_inconsistent(true);
        }

        lastcurvactiv = curvactiv->get_active();
    }
    */

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
    printf("activlumChanged\n");

    if (batchMode) {
        /*
        if (activlum->get_inconsistent()) {
            activlum->set_inconsistent(false);
            activlumConn.block(true);
            activlum->set_active(false);
            activlumConn.block(false);
        } else if (lastactivlum) {
            activlum->set_inconsistent(true);
        }

        lastactivlum = activlum->get_active();
        */
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
    printf("inversshaChanged\n");

    if (batchMode) {
        /*
        if (inverssha->get_inconsistent()) {
            inverssha->set_inconsistent(false);
            inversshaConn.block(true);
            inverssha->set_active(false);
            inversshaConn.block(false);
        } else if (lastinverssha) {
            inverssha->set_inconsistent(true);
        }

        lastinverssha = inverssha->get_active();
        */
    }

    // Update Sharpening GUI according to inverssha button state (to be compliant with updateSpecificGUIState function)
    if (inverssha->get_active()) {
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
    printf("inversretChanged\n");

    if (batchMode) {
        /*
        if (inversret->get_inconsistent()) {
            inversret->set_inconsistent(false);
            inversretConn.block(true);
            inversret->set_active(false);
            inversretConn.block(false);
        } else if (lastinversret) {
            inversret->set_inconsistent(true);
        }

        lastinversret = inversret->get_active();
        */
    }

    // Update Retinex GUI according to inversret button state (to be compliant with updateSpecificGUIState function)
    if (inversret->get_active()) {
        sensih->hide();
    } else {
        sensih->show();
    }

    if (getEnabled() && expsharp->getEnabled()) {
        if (listener) {
            if (inversret->get_active()) {
                listener->panelChanged(Evlocallabinversret, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinversret, M("GENERAL_DISABLED"));
            }
        }
    }
}

// TODO
void Locallab::setDefaults(const ProcParams * defParams, const ParamsEdited * pedited)
{
    /*
    degree->setDefault(defParams->locallab.degree);
    locY->setDefault(defParams->locallab.locY);
    locX->setDefault(defParams->locallab.locX);
    locYT->setDefault(defParams->locallab.locYT);
    locXL->setDefault(defParams->locallab.locXL);
    centerX->setDefault(defParams->locallab.centerX);
    centerY->setDefault(defParams->locallab.centerY);
    circrad->setDefault(defParams->locallab.circrad);
    adjblur->setDefault(defParams->locallab.adjblur);
    thres->setDefault(defParams->locallab.thres);
    proxi->setDefault(defParams->locallab.proxi);
    lightness->setDefault(defParams->locallab.lightness);
    contrast->setDefault(defParams->locallab.contrast);
    chroma->setDefault(defParams->locallab.chroma);
    warm->setDefault(defParams->locallab.warm);
    expcomp->setDefault(defParams->locallab.expcomp);
    black->setDefault(defParams->locallab.black);
    hlcompr->setDefault(defParams->locallab.hlcompr);
    hlcomprthresh->setDefault(defParams->locallab.hlcomprthresh);
    shcompr->setDefault(defParams->locallab.shcompr);

    noiselumf->setDefault(defParams->locallab.noiselumf);
    noiselumc->setDefault(defParams->locallab.noiselumc);
    noiselumdetail->setDefault(defParams->locallab.noiselumdetail);
    noiselequal->setDefault(defParams->locallab.noiselequal);
    noisechrodetail->setDefault(defParams->locallab.noisechrodetail);
    bilateral->setDefault(defParams->locallab.bilateral);
    sensiden->setDefault(defParams->locallab.sensiden);
    noisechrof->setDefault(defParams->locallab.noisechrof);
    noisechroc->setDefault(defParams->locallab.noisechroc);
    sharradius->setDefault(defParams->locallab.sharradius);
    sharamount->setDefault(defParams->locallab.sharamount);
    shardamping->setDefault(defParams->locallab.shardamping);
    shariter->setDefault(defParams->locallab.shariter);
    sensisha->setDefault(defParams->locallab.sensisha);
    sensi->setDefault(defParams->locallab.sensi);
    sensiex->setDefault(defParams->locallab.sensiex);
    sensih->setDefault(defParams->locallab.sensih);
    retrab->setDefault(defParams->locallab.retrab);
    sensiexclu->setDefault(defParams->locallab.sensiexclu);
    struc->setDefault(defParams->locallab.struc);
    sensicb->setDefault(defParams->locallab.sensicb);
    sensibn->setDefault(defParams->locallab.sensibn);
    sensitm->setDefault(defParams->locallab.sensitm);
    transit->setDefault(defParams->locallab.transit);
    radius->setDefault(defParams->locallab.radius);
    strength->setDefault(defParams->locallab.strength);
    stren->setDefault(defParams->locallab.stren);
    gamma->setDefault(defParams->locallab.gamma);
    estop->setDefault(defParams->locallab.estop);
    gamma->setDefault(defParams->locallab.gamma);
    scaltm->setDefault(defParams->locallab.scaltm);
    rewei->setDefault(defParams->locallab.rewei);
    neigh->setDefault(defParams->locallab.neigh);
    vart->setDefault(defParams->locallab.vart);
    chrrt->setDefault(defParams->locallab.chrrt);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setDefault(defParams->locallab.mult[i]);
    }

    threshold->setDefault(defParams->locallab.threshold);
    chromacbdl->setDefault(defParams->locallab.chromacbdl);
    pastels->setDefault(defParams->locallab.pastels);
    saturated->setDefault(defParams->locallab.saturated);
    psThreshold->setDefault<int> (defParams->locallab.psthreshold);
    sensiv->setDefault(defParams->locallab.sensiv);


    if (pedited) {
        degree->setDefaultEditedState(pedited->locallab.degree ? Edited : UnEdited);
        locY->setDefaultEditedState(pedited->locallab.locY ? Edited : UnEdited);
        locX->setDefaultEditedState(pedited->locallab.locX ? Edited : UnEdited);
        locYT->setDefaultEditedState(pedited->locallab.locYT ? Edited : UnEdited);
        locXL->setDefaultEditedState(pedited->locallab.locXL ? Edited : UnEdited);
        centerX->setDefaultEditedState(pedited->locallab.centerX ? Edited : UnEdited);
        centerY->setDefaultEditedState(pedited->locallab.centerY ? Edited : UnEdited);
        circrad->setDefaultEditedState(pedited->locallab.circrad ? Edited : UnEdited);
        adjblur->setDefaultEditedState(pedited->locallab.adjblur ? Edited : UnEdited);
        thres->setDefaultEditedState(pedited->locallab.thres ? Edited : UnEdited);
        proxi->setDefaultEditedState(pedited->locallab.proxi ? Edited : UnEdited);
        lightness->setDefaultEditedState(pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setDefaultEditedState(pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setDefaultEditedState(pedited->locallab.chroma ? Edited : UnEdited);
        warm->setDefaultEditedState(pedited->locallab.warm ? Edited : UnEdited);
        expcomp->setDefaultEditedState(pedited->locallab.expcomp ? Edited : UnEdited);
        black->setDefaultEditedState(pedited->locallab.black ? Edited : UnEdited);
        hlcompr->setDefaultEditedState(pedited->locallab.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setDefaultEditedState(pedited->locallab.hlcomprthresh ? Edited : UnEdited);
        shcompr->setDefaultEditedState(pedited->locallab.shcompr ? Edited : UnEdited);

        noiselumf->setDefaultEditedState(pedited->locallab.noiselumf ? Edited : UnEdited);
        noiselumc->setDefaultEditedState(pedited->locallab.noiselumc ? Edited : UnEdited);
        noiselumdetail->setDefaultEditedState(pedited->locallab.noiselumdetail ? Edited : UnEdited);
        noiselequal->setDefaultEditedState(pedited->locallab.noiselequal ? Edited : UnEdited);
        noisechrodetail->setDefaultEditedState(pedited->locallab.noisechrodetail ? Edited : UnEdited);
        bilateral->setDefaultEditedState(pedited->locallab.bilateral ? Edited : UnEdited);
        sensiden->setDefaultEditedState(pedited->locallab.sensiden ? Edited : UnEdited);
        noisechrof->setDefaultEditedState(pedited->locallab.noisechrof ? Edited : UnEdited);
        noisechroc->setDefaultEditedState(pedited->locallab.noisechroc ? Edited : UnEdited);
        sharradius->setDefaultEditedState(pedited->locallab.sharradius ? Edited : UnEdited);
        sharamount->setDefaultEditedState(pedited->locallab.sharamount ? Edited : UnEdited);
        shardamping->setDefaultEditedState(pedited->locallab.shardamping ? Edited : UnEdited);
        shariter->setDefaultEditedState(pedited->locallab.shariter ? Edited : UnEdited);
        sensisha->setDefaultEditedState(pedited->locallab.sensisha ? Edited : UnEdited);
        sensi->setDefaultEditedState(pedited->locallab.sensi ? Edited : UnEdited);
        sensiex->setDefaultEditedState(pedited->locallab.sensiex ? Edited : UnEdited);
        sensih->setDefaultEditedState(pedited->locallab.sensih ? Edited : UnEdited);
        retrab->setDefaultEditedState(pedited->locallab.retrab ? Edited : UnEdited);
        sensiexclu->setDefaultEditedState(pedited->locallab.sensiexclu ? Edited : UnEdited);
        struc->setDefaultEditedState(pedited->locallab.struc ? Edited : UnEdited);
        sensicb->setDefaultEditedState(pedited->locallab.sensicb ? Edited : UnEdited);
        sensibn->setDefaultEditedState(pedited->locallab.sensibn ? Edited : UnEdited);
        sensitm->setDefaultEditedState(pedited->locallab.sensitm ? Edited : UnEdited);
        radius->setDefaultEditedState(pedited->locallab.radius ? Edited : UnEdited);
        strength->setDefaultEditedState(pedited->locallab.strength ? Edited : UnEdited);
        stren->setDefaultEditedState(pedited->locallab.stren ? Edited : UnEdited);
        gamma->setDefaultEditedState(pedited->locallab.gamma ? Edited : UnEdited);
        estop->setDefaultEditedState(pedited->locallab.estop ? Edited : UnEdited);
        scaltm->setDefaultEditedState(pedited->locallab.scaltm ? Edited : UnEdited);
        rewei->setDefaultEditedState(pedited->locallab.rewei ? Edited : UnEdited);
        transit->setDefaultEditedState(pedited->locallab.transit ? Edited : UnEdited);
        str->setDefaultEditedState(pedited->locallab.str ? Edited : UnEdited);
        neigh->setDefaultEditedState(pedited->locallab.neigh ? Edited : UnEdited);
        vart->setDefaultEditedState(pedited->locallab.vart ? Edited : UnEdited);
        chrrt->setDefaultEditedState(pedited->locallab.chrrt ? Edited : UnEdited);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState(pedited->locallab.mult[i] ? Edited : UnEdited);
        }

        threshold->setDefaultEditedState(pedited->locallab.threshold ? Edited : UnEdited);
        chromacbdl->setDefaultEditedState(pedited->locallab.chromacbdl ? Edited : UnEdited);

        pastels->setDefaultEditedState(pedited->locallab.pastels ? Edited : UnEdited);
        saturated->setDefaultEditedState(pedited->locallab.saturated ? Edited : UnEdited);
        psThreshold->setDefaultEditedState(pedited->locallab.psthreshold ? Edited : UnEdited);
        sensiv->setDefaultEditedState(pedited->locallab.sensiv ? Edited : UnEdited);

    } else {
        degree->setDefaultEditedState(Irrelevant);
        locY->setDefaultEditedState(Irrelevant);
        locX->setDefaultEditedState(Irrelevant);
        locYT->setDefaultEditedState(Irrelevant);
        locXL->setDefaultEditedState(Irrelevant);
        centerX->setDefaultEditedState(Irrelevant);
        centerY->setDefaultEditedState(Irrelevant);
        circrad->setDefaultEditedState(Irrelevant);
        adjblur->setDefaultEditedState(Irrelevant);
        thres->setDefaultEditedState(Irrelevant);
        proxi->setDefaultEditedState(Irrelevant);
        lightness->setDefaultEditedState(Irrelevant);
        contrast->setDefaultEditedState(Irrelevant);
        chroma->setDefaultEditedState(Irrelevant);
        warm->setDefaultEditedState(Irrelevant);
        expcomp->setDefaultEditedState(Irrelevant);
        black->setDefaultEditedState(Irrelevant);
        hlcompr->setDefaultEditedState(Irrelevant);
        hlcomprthresh->setDefaultEditedState(Irrelevant);
        shcompr->setDefaultEditedState(Irrelevant);

        noiselumf->setDefaultEditedState(Irrelevant);
        noiselumc->setDefaultEditedState(Irrelevant);
        noiselumdetail->setDefaultEditedState(Irrelevant);
        noiselequal->setDefaultEditedState(Irrelevant);
        noisechrodetail->setDefaultEditedState(Irrelevant);
        bilateral->setDefaultEditedState(Irrelevant);
        sensiden->setDefaultEditedState(Irrelevant);
        noisechrof->setDefaultEditedState(Irrelevant);
        noisechroc->setDefaultEditedState(Irrelevant);
        sharradius->setDefaultEditedState(Irrelevant);
        sharamount->setDefaultEditedState(Irrelevant);
        shardamping->setDefaultEditedState(Irrelevant);
        shariter->setDefaultEditedState(Irrelevant);
        sensisha->setDefaultEditedState(Irrelevant);
        sensi->setDefaultEditedState(Irrelevant);
        sensiex->setDefaultEditedState(Irrelevant);
        sensih->setDefaultEditedState(Irrelevant);
        retrab->setDefaultEditedState(Irrelevant);
        sensicb->setDefaultEditedState(Irrelevant);
        sensiexclu->setDefaultEditedState(Irrelevant);
        sensicb->setDefaultEditedState(Irrelevant);
        struc->setDefaultEditedState(Irrelevant);
        sensitm->setDefaultEditedState(Irrelevant);
        radius->setDefaultEditedState(Irrelevant);
        strength->setDefaultEditedState(Irrelevant);
        stren->setDefaultEditedState(Irrelevant);
        gamma->setDefaultEditedState(Irrelevant);
        estop->setDefaultEditedState(Irrelevant);
        scaltm->setDefaultEditedState(Irrelevant);
        rewei->setDefaultEditedState(Irrelevant);
        transit->setDefaultEditedState(Irrelevant);
        str->setDefaultEditedState(Irrelevant);
        neigh->setDefaultEditedState(Irrelevant);
        vart->setDefaultEditedState(Irrelevant);
        chrrt->setDefaultEditedState(Irrelevant);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState(Irrelevant);
        }

        threshold->setDefaultEditedState(Irrelevant);
        chromacbdl->setDefaultEditedState(Irrelevant);

        pastels->setDefaultEditedState(Irrelevant);
        saturated->setDefaultEditedState(Irrelevant);
        psThreshold->setDefaultEditedState(Irrelevant);
        sensiv->setDefaultEditedState(Irrelevant);

    }
    */
}

void Locallab::adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop)
{
    printf("adjusterChangedTS\n");

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabPastSatThreshold, psThreshold->getHistoryString());
        }
    }
}

void Locallab::adjusterChanged(Adjuster * a, double newval)
{
    printf("adjusterChanged\n");

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
        shcompr->set_sensitive(!((int)black->getValue() == 0)); // At black = 0, shcompr value has no effect
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
    if (a == pastels && pastSatTog->get_active()) {
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

        if (a == sensih) {
            if (listener) {
                listener->panelChanged(Evlocallabsensih, sensih->getTextValue());
            }
        }
    }

    // Sharpening
    if (getEnabled() && expsharp->getEnabled()) {
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

        if (a == sensisha) {
            if (listener) {
                listener->panelChanged(Evlocallabsensis, sensisha->getTextValue());
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
    printf("avoidChanged\n");

    if (batchMode) {
        /*
        if (avoid->get_inconsistent()) {
            avoid->set_inconsistent(false);
            avoidConn.block(true);
            avoid->set_active(false);
            avoidConn.block(false);
        } else if (lastavoid) {
            avoid->set_inconsistent(true);
        }

        lastavoid = avoid->get_active();
        */
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

// TODO
void Locallab::trimValues(rtengine::procparams::ProcParams * pp)
{
    /*
    degree->trimValue(pp->locallab.degree);
    locY->trimValue(pp->locallab.locY);
    locX->trimValue(pp->locallab.locX);
    locYT->trimValue(pp->locallab.locYT);
    locXL->trimValue(pp->locallab.locXL);
    centerX->trimValue(pp->locallab.centerX);
    centerY->trimValue(pp->locallab.centerY);
    circrad->trimValue(pp->locallab.circrad);
    adjblur->trimValue(pp->locallab.adjblur);
    thres->trimValue(pp->locallab.thres);
    proxi->trimValue(pp->locallab.proxi);
    lightness->trimValue(pp->locallab.lightness);
    contrast->trimValue(pp->locallab.contrast);
    chroma->trimValue(pp->locallab.chroma);
    warm->trimValue(pp->locallab.warm);
    expcomp->trimValue(pp->locallab.expcomp);
    hlcompr->trimValue(pp->locallab.hlcompr);
    hlcomprthresh->trimValue(pp->locallab.hlcomprthresh);
    black->trimValue(pp->locallab.black);
    shcompr->trimValue(pp->locallab.shcompr);
    noiselumf->trimValue(pp->locallab.noiselumf);
    noiselumc->trimValue(pp->locallab.noiselumc);
    noiselumdetail->trimValue(pp->locallab.noiselumdetail);
    noiselequal->trimValue(pp->locallab.noiselequal);
    noisechrodetail->trimValue(pp->locallab.noisechrodetail);
    bilateral->trimValue(pp->locallab.bilateral);
    sensiden->trimValue(pp->locallab.sensiden);
    noisechrof->trimValue(pp->locallab.noisechrof);
    noisechroc->trimValue(pp->locallab.noisechroc);
    sharradius->trimValue(pp->locallab.sharradius);
    sharamount->trimValue(pp->locallab.sharamount);
    shardamping->trimValue(pp->locallab.shardamping);
    shariter->trimValue(pp->locallab.shariter);
    sensisha->trimValue(pp->locallab.sensisha);
    sensi->trimValue(pp->locallab.sensi);
    sensiex->trimValue(pp->locallab.sensiex);
    sensih->trimValue(pp->locallab.sensih);
    sensiexclu->trimValue(pp->locallab.sensiexclu);
    struc->trimValue(pp->locallab.struc);
    sensicb->trimValue(pp->locallab.sensicb);
    sensibn->trimValue(pp->locallab.sensibn);
    sensitm->trimValue(pp->locallab.sensitm);
    radius->trimValue(pp->locallab.radius);
    strength->trimValue(pp->locallab.strength);
    stren->trimValue(pp->locallab.stren);
    gamma->trimValue(pp->locallab.gamma);
    estop->trimValue(pp->locallab.estop);
    scaltm->trimValue(pp->locallab.scaltm);
    rewei->trimValue(pp->locallab.rewei);
    transit->trimValue(pp->locallab.transit);
    str->trimValue(pp->locallab.str);
    neigh->trimValue(pp->locallab.neigh);
    vart->trimValue(pp->locallab.vart);
    chrrt->trimValue(pp->locallab.chrrt);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->trimValue(pp->locallab.mult[i]);
    }

    chromacbdl->trimValue(pp->locallab.chromacbdl);

    threshold->trimValue(pp->locallab.threshold);
    pastels->trimValue(pp->locallab.pastels);
    saturated->trimValue(pp->locallab.saturated);
    sensiv->trimValue(pp->locallab.sensiv);
     */
}

void Locallab::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);
    printf("BatchMode : %d\n", batchMode);

    adjblur->showEditedCB();
    lightness->showEditedCB();
    contrast->showEditedCB();
    chroma->showEditedCB();
    warm->showEditedCB();
    expcomp->showEditedCB();
    black->showEditedCB();
    hlcompr->showEditedCB();
    hlcomprthresh->showEditedCB();
    shcompr->showEditedCB();
    noiselumf->showEditedCB();
    noiselumc->showEditedCB();
    noiselumdetail->showEditedCB();
    noiselequal->showEditedCB();
    noisechrodetail->showEditedCB();
    bilateral->showEditedCB();
    sensiden->showEditedCB();
    noisechroc->showEditedCB();
    noiselumf->showEditedCB();
    sharradius->showEditedCB();
    sharamount->showEditedCB();
    shardamping->showEditedCB();
    shariter->showEditedCB();
    sensisha->showEditedCB();
    sensi->showEditedCB();
    sensiex->showEditedCB();
    sensih->showEditedCB();
    sensicb->showEditedCB();
    sensibn->showEditedCB();
    sensitm->showEditedCB();
    radius->showEditedCB();
    strength->showEditedCB();
    stren->showEditedCB();
    gamma->showEditedCB();
    estop->showEditedCB();
    scaltm->showEditedCB();
    rewei->showEditedCB();
    str->showEditedCB();
    neigh->showEditedCB();
    vart->showEditedCB();
    LocalcurveEditorgainT->setBatchMode(batchMode);
    llCurveEditorG->setBatchMode(batchMode);
    chrrt->showEditedCB();

    for (int i = 0; i < 5; i++) {
        multiplier[i]->showEditedCB();
    }

    threshold->showEditedCB();
    chromacbdl->showEditedCB();
    pastels->showEditedCB();
    saturated->showEditedCB();
    psThreshold->showEditedCB();
    sensiv->showEditedCB();

    curveEditorGG->setBatchMode(batchMode);

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
    // cTgainshape->setEditProvider(provider);
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
    printf("enableListener\n");

    FoldableToolPanel::enableListener();
    // Color & Light
    enablecolorConn.block(false);
    curvactivConn.block(false);
    qualitycurveMethodConn.block(false);
    inversConn.block(false);
    // Exposure
    enableexposeConn.block(false);
    // Vibrance
    enablevibranceConn.block(false);
    pskinsconn.block(false);
    ashiftconn.block(false);
    pastsattogconn.block(false);
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
    // Contrast by detail levels
    enablecbdlConn.block(false);
    // Denoise
    enabledenoiConn.block(false);
    avoidConn.block(false);
}

void Locallab::disableListener()
{
    printf("disableListener\n");

    FoldableToolPanel::disableListener();
    // Color & Light
    enablecolorConn.block(true);
    curvactivConn.block(true);
    qualitycurveMethodConn.block(true);
    inversConn.block(true);
    // Exposure
    enableexposeConn.block(true);
    // Vibrance
    enablevibranceConn.block(true);
    pskinsconn.block(true);
    ashiftconn.block(true);
    pastsattogconn.block(true);
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
    // Contrast by detail levels
    enablecbdlConn.block(true);
    // Denoise
    enabledenoiConn.block(true);
    avoidConn.block(true);
}

void Locallab::updateLocallabGUI(const rtengine::procparams::ProcParams* pp, int index)
{
    printf("updateLocallabGUI\n");

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
        } else {
            qualitycurveMethod->set_active(2);
        }

        llshape->setCurve(pp->locallab.spots.at(index).llcurve);
        ccshape->setCurve(pp->locallab.spots.at(index).cccurve);
        LHshape->setCurve(pp->locallab.spots.at(index).LHcurve);
        HHshape->setCurve(pp->locallab.spots.at(index).HHcurve);
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
        sensih->setValue(pp->locallab.spots.at(index).sensih);
        cTgainshape->setCurve(pp->locallab.spots.at(index).localTgaincurve);
        inversret->set_active(pp->locallab.spots.at(index).inversret);

        // Sharpening
        expsharp->setEnabled(pp->locallab.spots.at(index).expsharp);
        sharradius->setValue(pp->locallab.spots.at(index).sharradius);
        sharamount->setValue(pp->locallab.spots.at(index).sharamount);
        shardamping->setValue(pp->locallab.spots.at(index).shardamping);
        shariter->setValue(pp->locallab.spots.at(index).shariter);
        sensisha->setValue(pp->locallab.spots.at(index).sensisha);
        inverssha->set_active(pp->locallab.spots.at(index).inverssha);

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
    }
}

void Locallab::updateSpecificGUIState()
{
    // Update Color & Light GUI according to invers button state
    if (invers->get_active()) {
        sensi->show();
        llCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        labqualcurv->hide();

    } else {
        sensi->show();
        llCurveEditorG->show();
        curvactiv->show();
        qualitycurveMethod->show();
        labqualcurv->show();

    }

    // Update Exposure GUI according to black adjuster state (to be compliant with adjusterChanged function)
    shcompr->set_sensitive(!((int)black->getValue() == 0)); // At black = 0, shcompr value has no effect

    // Update Vibrance GUI according to pastsattog button state (to be compliant with pastsattog_toggled function)
    if (pastSatTog->get_active()) {
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
    if (blurMethod->get_active_row_number() == 0 || blurMethod->get_active_row_number() == 2) {
        sensibn->show();
    } else {
        sensibn->hide();
    }

    // Update Retinex GUI according to inversret button state (to be compliant with inversretChanged function)
    if (inversret->get_active()) {
        sensih->hide();
    } else {
        sensih->show();
    }

    // Update Sharpening GUI according to inverssha button state (to be compliant with inversshaChanged function)
    if (inverssha->get_active()) {
        sensisha->hide();
    } else {
        sensisha->show();
    }
}

void Locallab::autoOpenCurve()
{
    printf("autoOpenCurve\n");

    // TODO autoOpenCurve only considers linearity state of selected spot curve
    llshape->openIfNonlinear();
    ccshape->openIfNonlinear();
    LHshape->openIfNonlinear();
    HHshape->openIfNonlinear();
}
