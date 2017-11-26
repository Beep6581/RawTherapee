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

using namespace rtengine;
using namespace rtengine::procparams;
extern Options options;



Locallab::Locallab():
    FoldableToolPanel(this, "locallab", M("TP_LOCALLAB_LABEL"), true, true),
    EditSubscriber(ET_OBJECTS), lastObject(-1),
    expcolor(new MyExpander(true, M("TP_LOCALLAB_COFR"))),
    expexpose(new MyExpander(true, M("TP_LOCALLAB_EXPOSE"))),
    expvibrance(new MyExpander(true, M("TP_LOCALLAB_VIBRANCE"))),
    expblur(new MyExpander(true, M("TP_LOCALLAB_BLUFR"))),
    exptonemap(new MyExpander(true, M("TP_LOCALLAB_TM"))),
    expreti(new MyExpander(true, M("TP_LOCALLAB_RETI"))),
    expsharp(new MyExpander(true, M("TP_LOCALLAB_SHARP"))),
    expcbdl(new MyExpander(true, M("TP_LOCALLAB_CBDL"))),
    expdenoi(new MyExpander(true, M("TP_LOCALLAB_DENOIS"))),
    expsettings(new MyExpander(false, M("TP_LOCALLAB_SETTINGS"))),

    LocalcurveEditorgainT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAIN"))),
    LocalcurveEditorgainTrab(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAINRAB"))),
    llCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_LUM"))),


    anbspot(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ANBSPOT"), 0, 1, 1, 0))),
    locX(Gtk::manage(new Adjuster(M("TP_LOCAL_WIDTH"), 0, 2250, 1, 250))),
    locXL(Gtk::manage(new Adjuster(M("TP_LOCAL_WIDTH_L"), 0, 2250, 1, 250))),
    degree(Gtk::manage(new Adjuster(M("TP_LOCAL_DEGREE"), -180, 180, 1, 0))),
    locY(Gtk::manage(new Adjuster(M("TP_LOCAL_HEIGHT"), 0, 2250, 1, 250))),
    locYT(Gtk::manage(new Adjuster(M("TP_LOCAL_HEIGHT_T"), 0, 2250, 1, 250))),
    centerX(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTER_X"), -1000, 1000, 1, 0))),
    centerY(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTER_Y"), -1000, 1000, 1, 0))),
    circrad(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CIRCRADIUS"), 2, 150, 1, 18))),
    sensiexclu(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIEXCLU"), 0, 100, 1, 19))),
    struc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUC"), 0, 5, 1, 0))),
    thres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_THRES"), 1, 35, 1, 18))),
    proxi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_PROXI"), 0, 60, 1, 0))),
    lightness(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LIGHTNESS"), -100, 100, 1, 0))),
    contrast(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0))),
    chroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMA"), -100, 150, 1, 0))),
    sensi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    expcomp(Gtk::manage(new Adjuster(M("TP_EXPOSURE_EXPCOMP"), -200, 200, 5, 0))),
    hlcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 20))),
    hlcomprthresh(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 33))),
    black(Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0))),
    shcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50))),
    sensiex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    radius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADIUS"), 1, 100, 1, 1))),
    strength(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENGTH"), 0, 100, 1, 0))),
    sensibn(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIBN"), 0, 100, 1, 40))),
    transit(Gtk::manage(new Adjuster(M("TP_LOCALLAB_TRANSIT"), 5, 95, 1, 60))),
    stren(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STREN"), -50, 100, 1, 0))),
    gamma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAM"), 80, 150, 1, 100))),
    estop(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ESTOP"), 10, 400, 1, 140))),
    scaltm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALTM"), 1, 100, 1, 10))),
    rewei(Gtk::manage(new Adjuster(M("TP_LOCALLAB_REWEI"), 0, 9, 1, 0))),
    sensitm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    str(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STR"), 0, 100, 1, 0))),
    neigh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NEIGH"), 14, 150, 1, 50))),
    vart(Gtk::manage(new Adjuster(M("TP_LOCALLAB_VART"), 50, 500, 1, 200))),
    chrrt(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHRRT"), 0, 100, 1, 0))),
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIH"), 0, 100, 1, 19))),
    retrab(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RETRAB"), 0, 10000, 1, 500))),
    chromacbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMACBDL"), 0, 300, 1, 0))),
    threshold(Gtk::manage(new Adjuster(M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 100, 1, 20))),
    sensicb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSICB"), 0, 100, 1, 19))),
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 42, 500, 1, 4))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 75))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 75))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    noiselumf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINE"), 0, 100, 1, 0))),
    noiselumc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMCOARSE"), 0, 100, 1, 0))),
    noisechrof(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROFINE"), 0, 100, 1, 0))),
    noisechroc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROCOARSE"), 0, 100, 1, 0))),
    hueref(Gtk::manage(new Adjuster(M("TP_LOCALLAB_HUEREF"), -3.15, 3.15, 0.01, 0))),
    chromaref(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMAREF"), 0, 200, 0.01, 0))),
    lumaref(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LUMAMAREF"), 0, 100, 0.01, 0))),
    sobelref(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOBELREF"), 0, 100, 0.01, 0))),
    centerXbuf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTERBUF_X"), -1000, 1000, 1, 0))),
    centerYbuf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTERBUF_Y"), -1000, 1000, 1, 0))),
    adjblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ADJBLUR"), 0, 100, 1, 0))),

    Smethod(Gtk::manage(new MyComboBoxText())),
    Exclumethod(Gtk::manage(new MyComboBoxText())),

    retinexMethod(Gtk::manage(new MyComboBoxText())),
    qualityMethod(Gtk::manage(new MyComboBoxText())),
    qualitycurveMethod(Gtk::manage(new MyComboBoxText())),
    blurMethod(Gtk::manage(new MyComboBoxText())),
    dustMethod(Gtk::manage(new MyComboBoxText())),

    excluFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_EXCLUF")))),
    artifFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_ARTIF")))),
    shapeFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHFR")))),
    superFrame(Gtk::manage(new Gtk::Frame())),
    dustFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_DUST")))),

    //  artifVBox (Gtk::manage (new Gtk::VBox ())),
//   shapeVBox (Gtk::manage (new Gtk::VBox ())),
//   tmBox (Gtk::manage (new Gtk::VBox())),
    //  retiBox (Gtk::manage (new Gtk::VBox())),
    //  colorVBox (Gtk::manage ( new Gtk::VBox())),
    // blurrVBox (Gtk::manage ( new Gtk::VBox())),
//   sharpVBox (Gtk::manage ( new Gtk::VBox())),
//   cbdlVBox (Gtk::manage ( new Gtk::VBox())),
//   denoisVBox (Gtk::manage ( new Gtk::VBox())),
//   superVBox (Gtk::manage (new Gtk::VBox ())),


    labmdh(Gtk::manage(new Gtk::Label(M("TP_LOCRETI_METHOD") + ":"))),
    labqual(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_QUAL_METHOD") + ":"))),
    labqualcurv(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_QUALCURV_METHOD") + ":"))),
    labmS(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_STYPE") + ":"))),
    labmEx(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_EXCLUTYPE") + ":"))),

    ctboxS(Gtk::manage(new Gtk::HBox())),
    ctboxEx(Gtk::manage(new Gtk::HBox())),
    dhbox(Gtk::manage(new Gtk::HBox())),
    qualbox(Gtk::manage(new Gtk::HBox())),
    qualcurvbox(Gtk::manage(new Gtk::HBox())),

    avoid(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_AVOID")))),
    activlum(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ACTIV")))),
    invers(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    curvactiv(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CURV")))),
    inversrad(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    inversret(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    inverssha(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    cutpast(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CUTPAST")))),
    lastdust(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_LASTDUST")))),
    draggedPointOldAngle(-1000.)

{
    CurveListener::setMulti(true);
    ProcParams params;
    editHBox = Gtk::manage(new Gtk::HBox());
    edit = Gtk::manage(new Gtk::ToggleButton());
    edit->add(*Gtk::manage(new RTImage("editmodehand.png")));
    edit->set_tooltip_text(M("EDIT_OBJECT_TOOLTIP"));
    editConn = edit->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::editToggled));
    editHBox->pack_start(*edit, Gtk::PACK_SHRINK, 0);
    pack_start(*editHBox, Gtk::PACK_SHRINK, 0);
    int realnbspot;


    realnbspot = options.rtSettings.nspot;
    nbspot = Gtk::manage(new Adjuster(M("TP_LOCALLAB_NBSPOT"), 1, realnbspot, 1, 1));

    if (options.rtSettings.locdelay) {

        if (nbspot->delay < 200) {
            nbspot->delay = 200;
        }
    }


    nbspot->setAdjusterListener(this);
    nbspot->set_tooltip_text(M("TP_LOCALLAB_NBSPOT_TOOLTIP"));


    anbspot->setAdjusterListener(this);
    anbspot->set_tooltip_text(M("TP_LOCALLAB_ANBSPOT_TOOLTIP"));

    shapeFrame->set_label_align(0.025, 0.5);

    expsettings->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsettings));


    expcolor->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcolor));
    enablecolorConn = expcolor->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcolor));

    expexpose->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expexpose));
    enableexposeConn = expexpose->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expexpose));

    expvibrance->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expvibrance));
    enablevibranceConn = expvibrance->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expvibrance));

    expblur->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expblur));
    enableblurConn = expblur->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expblur));

    exptonemap->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), exptonemap));
    enabletonemapConn = exptonemap->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), exptonemap));

    expreti->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expreti));
    enableretiConn = expreti->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expreti));

    expsharp->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsharp));
    enablesharpConn = expsharp->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expsharp));

    expcbdl->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcbdl));
    enablecbdlConn = expcbdl->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcbdl));

    expdenoi->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expdenoi));
    enabledenoiConn = expdenoi->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expdenoi));

    ctboxEx->pack_start(*labmEx, Gtk::PACK_SHRINK, 4);
    ctboxEx->set_tooltip_markup(M("TP_LOCALLAB_EXCLUTYPE_TOOLTIP"));

    Exclumethod->append(M("TP_LOCALLAB_EXNORM"));
    Exclumethod->append(M("TP_LOCALLAB_EXECLU"));
    Exclumethod->set_active(0);
    Exclumethodconn = Exclumethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::ExclumethodChanged));

    sensiexclu->set_tooltip_text(M("TP_LOCALLAB_SENSIEXCLU_TOOLTIP"));
    sensiexclu->setAdjusterListener(this);

    struc->set_tooltip_text(M("TP_LOCALLAB_STRUC_TOOLTIP"));
    struc->setAdjusterListener(this);

    ctboxS->pack_start(*labmS, Gtk::PACK_SHRINK, 4);
    ctboxS->set_tooltip_markup(M("TP_LOCALLAB_STYPE_TOOLTIP"));

    Smethod->append(M("TP_LOCALLAB_IND"));
    Smethod->append(M("TP_LOCALLAB_SYM"));
    Smethod->append(M("TP_LOCALLAB_INDSL"));
    Smethod->append(M("TP_LOCALLAB_SYMSL"));
    Smethod->set_active(0);
    Smethodconn = Smethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::SmethodChanged));


    //locX->set_tooltip_text (M("TP_LOCAL_WIDTH_TOOLTIP"));
    locX->setAdjusterListener(this);

    //locX->set_tooltip_text (M("TP_LOCAL_WIDTH_TOOLTIP"));
    locXL->setAdjusterListener(this);

    //degree->set_tooltip_text (M("TP_LOCAL_DEGREE_TOOLTIP"));
    degree->setAdjusterListener(this);

    //locY->set_tooltip_text (M("TP_LOCAL_HEIGHT_TOOLTIP"));
    locY->setAdjusterListener(this);

    //locY->set_tooltip_text (M("TP_LOCAL_HEIGHT_TOOLTIP"));
    locYT->setAdjusterListener(this);

    //centerX->set_tooltip_text (M("TP_LOCALLAB_CENTER_X_TOOLTIP"));
    centerX->setAdjusterListener(this);

    //centerY->set_tooltip_text (M("TP_LOCALLAB_CENTER_Y_TOOLTIP"));
    centerY->setAdjusterListener(this);

    circrad->setAdjusterListener(this);


    qualityMethod->append(M("TP_LOCALLAB_STD"));
    qualityMethod->append(M("TP_LOCALLAB_ENH"));
    qualityMethod->append(M("TP_LOCALLAB_ENHDEN"));
    qualityMethod->set_active(0);
    qualityMethodConn = qualityMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::qualityMethodChanged));
    qualityMethod->set_tooltip_markup(M("TP_LOCALLAB_METHOD_TOOLTIP"));

    thres->setAdjusterListener(this);

    proxi->setAdjusterListener(this);
    std::vector<GradientMilestone> milestones;
    std::vector<double> defaultCurve;
    std::vector<double> defaultCurve2;
    std::vector<double> defaultCurve2rab;
    std::vector<double> defaultCurve3;
    std::vector<double> defaultCurve4;
    std::vector<double> defaultCurve5;

    irg   = Gtk::manage(new RTImage("Chanmixer-RG.png"));

    qualitycurveMethod->append(M("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVENH"));
    qualitycurveMethod->set_active(0);
    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::qualitycurveMethodChanged));
    qualitycurveMethod->set_tooltip_markup(M("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));



    llCurveEditorG->setCurveListener(this);

    rtengine::LocallabParams::getDefaultLLCurve(defaultCurve);
    llshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "L(L)"));
    llshape->setResetCurve(DiagonalCurveType(defaultCurve.at(0)), defaultCurve);
    llshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    milestones.push_back(GradientMilestone(0., 0., 0., 0.));
    milestones.push_back(GradientMilestone(1., 1., 1., 1.));
    llshape->setBottomBarBgGradient(milestones);
    llshape->setLeftBarBgGradient(milestones);

    rtengine::LocallabParams::getDefaultCCCurve(defaultCurve4);
    ccshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "C(C)"));
    ccshape->setResetCurve(DiagonalCurveType(defaultCurve4.at(0)), defaultCurve4);
    ccshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    milestones.push_back(GradientMilestone(0., 0., 0., 0.));
    milestones.push_back(GradientMilestone(1., 1., 1., 1.));
    ccshape->setBottomBarBgGradient(milestones);
    ccshape->setLeftBarBgGradient(milestones);

    rtengine::LocallabParams::getDefaultLHCurve(defaultCurve3);

    LHshape = static_cast<FlatCurveEditor*>(llCurveEditorG->addCurve(CT_Flat, "L(H)", nullptr, false, true));

    LHshape->setIdentityValue(0.);
    LHshape->setResetCurve(FlatCurveType(defaultCurve3.at(0)), defaultCurve3);
    LHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    LHshape->setCurveColorProvider(this, 1);
    milestones.clear();

    for (int i = 0; i < 7; i++) {
        float R, G, B;
        float x = float (i) * (1.0f / 6.0);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        milestones.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    LHshape->setBottomBarBgGradient(milestones);

    rtengine::LocallabParams::getDefaultHHCurve(defaultCurve5);

    HHshape = static_cast<FlatCurveEditor*>(llCurveEditorG->addCurve(CT_Flat, "H(H)", nullptr, false, true));

    HHshape->setIdentityValue(0.);
    HHshape->setResetCurve(FlatCurveType(defaultCurve5.at(0)), defaultCurve5);
    HHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    HHshape->setCurveColorProvider(this, 1);
    milestones.clear();

    for (int i = 0; i < 7; i++) {
        float R, G, B;
        float x = float (i) * (1.0f / 6.0);

        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        milestones.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    HHshape->setBottomBarBgGradient(milestones);


    llCurveEditorG->curveListComplete();



    //lightness->set_tooltip_text (M("TP_LOCALLAB_LIGHTNESS_TOOLTIP"));
    lightness->setAdjusterListener(this);

    //contrast->set_tooltip_text (M("TP_LOCALLAB_CONTRAST_TOOLTIP"));
    contrast->setAdjusterListener(this);

    //chroma->set_tooltip_text (M("TP_LOCALLAB_CHROMA_TOOLTIP"));
    chroma->setAdjusterListener(this);

    sensi->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensi->setAdjusterListener(this);

    centerXbuf->setAdjusterListener(this);;
    centerYbuf->setAdjusterListener(this);;
    adjblur->setAdjusterListener(this);;

//exposure

    expcomp->setAdjusterListener(this);
    hlcomprthresh->setAdjusterListener(this);
    black->setAdjusterListener(this);
    hlcompr->setAdjusterListener(this);
    shcompr->setAdjusterListener(this);
    sensiex->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensiex->setAdjusterListener(this);

    //radius->set_tooltip_text (M("TP_LOCALLAB_RADIUS_TOOLTIP"));
    radius->setAdjusterListener(this);
    //radius->set_tooltip_text (M("TP_LOCALLAB_RADIUS_TOOLTIP"));
    strength->setAdjusterListener(this);


    sensibn->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensibn->setAdjusterListener(this);

    activlum->set_active(false);
    activlumConn  = activlum->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::activlumChanged));

    transit->set_tooltip_text(M("TP_LOCALLAB_TRANSIT_TOOLTIP"));
    transit->setAdjusterListener(this);

    invers->set_active(false);
    inversConn  = invers->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversChanged));

    curvactiv->set_active(false);
    curvactivConn  = curvactiv->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::curvactivChanged));

    inversrad->set_active(false);
    inversradConn  = inversrad->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversradChanged));

    inversret->set_active(false);
    inversretConn  = inversret->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversretChanged));

    cutpast->set_active(false);
    cutpastConn  = cutpast->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::cutpastChanged));
    cutpast->set_tooltip_text(M("TP_LOCALLAB_CUTPAST_TOOLTIP"));
//tone mapping local

    lastdust->set_active(false);
    lastdustConn  = lastdust->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::lastdustChanged));
    lastdust->set_tooltip_text(M("TP_LOCALLAB_LASTDUST_TOOLTIP"));

    stren->setAdjusterListener(this);

    gamma->setAdjusterListener(this);

    estop->setAdjusterListener(this);

    scaltm->setAdjusterListener(this);

    rewei->setAdjusterListener(this);

    sensitm->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensitm->setAdjusterListener(this);

//end TM


//retinex local

    dhbox->pack_start(*labmdh, Gtk::PACK_SHRINK, 1);

    retinexMethod->append(M("TP_RETINEX_LOW"));
    retinexMethod->append(M("TP_RETINEX_UNIFORM"));
    retinexMethod->append(M("TP_RETINEX_HIGH"));
    retinexMethod->set_active(0);
    retinexMethodConn = retinexMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::retinexMethodChanged));
    retinexMethod->set_tooltip_markup(M("TP_LOCRETI_METHOD_TOOLTIP"));

    str->setAdjusterListener(this);
    neigh->setAdjusterListener(this);
    vart->setAdjusterListener(this);
    chrrt->setAdjusterListener(this);
    sensih->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensih->setAdjusterListener(this);
    retrab->setAdjusterListener(this);


    LocalcurveEditorgainT->setCurveListener(this);
    rtengine::LocallabParams::getDefaultLocalgainCurveT(defaultCurve2);


    cTgainshape = static_cast<FlatCurveEditor*>(LocalcurveEditorgainT->addCurve(CT_Flat, "", nullptr, false, false));

    cTgainshape->setIdentityValue(0.);
    cTgainshape->setResetCurve(FlatCurveType(defaultCurve2.at(0)), defaultCurve2);
    cTgainshape->setTooltip(M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));

    LocalcurveEditorgainTrab->setCurveListener(this);

    rtengine::LocallabParams::getDefaultLocalgainCurveTrab(defaultCurve2rab);


    cTgainshaperab = static_cast<FlatCurveEditor*>(LocalcurveEditorgainTrab->addCurve(CT_Flat, "", nullptr, false, false));


    cTgainshaperab->setIdentityValue(0.);
    cTgainshaperab->setResetCurve(FlatCurveType(defaultCurve2rab.at(0)), defaultCurve2rab);
    cTgainshaperab->setTooltip(M("TP_RETINEX_GAINTRANSMISSIONRAB_TOOLTIP"));

    LocalcurveEditorgainT->curveListComplete();
    LocalcurveEditorgainT->show();
    LocalcurveEditorgainTrab->curveListComplete();
    LocalcurveEditorgainTrab->show();


// end reti
    ToolParamBlock* const shapeBox = Gtk::manage(new ToolParamBlock());
    avoid->set_active(false);
    avoidConn  = avoid->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::avoidChanged));
    shapeBox->pack_start(*nbspot);
    pack_start(*anbspot);

    hueref->setAdjusterListener(this);
    chromaref->setAdjusterListener(this);
    lumaref->setAdjusterListener(this);
    sobelref->setAdjusterListener(this);

    pack_start(*hueref);
    pack_start(*chromaref);
    pack_start(*lumaref);
    pack_start(*sobelref);

    anbspot->hide();//keep anbspot  - i used it to test diffrent algo...
    hueref->hide();
    chromaref->hide();
    lumaref->hide();
    sobelref->hide();
    ctboxEx->pack_start(*Exclumethod);
    shapeBox->pack_start(*ctboxEx);

    excluFrame->set_label_align(0.025, 0.5);
    excluFrame->set_tooltip_text(M("TP_LOCALLAB_EXCLUF_TOOLTIP"));
    ToolParamBlock* const excluBox = Gtk::manage(new ToolParamBlock());

    excluBox->pack_start(*sensiexclu);
    //excluBox->pack_start (*struc);
    excluFrame->add(*excluBox);
    shapeBox->pack_start(*excluFrame);


    ctboxS->pack_start(*Smethod);
    shapeBox->pack_start(*ctboxS);


    shapeBox->pack_start(*locX);
    shapeBox->pack_start(*locXL);
    //pack_start (*degree);
    shapeBox->pack_start(*locY);
    shapeBox->pack_start(*locYT);
    shapeBox->pack_start(*centerX);
    shapeBox->pack_start(*centerY);
    shapeBox->pack_start(*circrad);
    qualbox->pack_start(*labqual, Gtk::PACK_SHRINK, 4);
    qualbox->pack_start(*qualityMethod);
    shapeBox->pack_start(*qualbox);
    shapeBox->pack_start(*transit);

    artifFrame->set_label_align(0.025, 0.5);
    artifFrame->set_tooltip_text(M("TP_LOCALLAB_ARTIF_TOOLTIP"));

    ToolParamBlock* const artifBox = Gtk::manage(new ToolParamBlock());

    artifBox->pack_start(*thres);
    artifBox->pack_start(*proxi);
    artifFrame->add(*artifBox);
    shapeBox->pack_start(*artifFrame);

    expsettings->add(*shapeBox);
    expsettings->setLevel(2);
    pack_start(*expsettings);



    Gtk::HBox * buttonBox1 = Gtk::manage(new Gtk::HBox(true, 10));

    Gtk::Button * lumacontrastMinusButton = Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")));
    buttonBox1->pack_start(*lumacontrastMinusButton);
    lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumacontrastMinusPressed));

    Gtk::Button * lumaneutralButton = Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")));
    buttonBox1->pack_start(*lumaneutralButton);
    lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumaneutralPressed));

    Gtk::Button * lumacontrastPlusButton = Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS")));
    buttonBox1->pack_start(*lumacontrastPlusButton);
    lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumacontrastPlusPressed));
    ToolParamBlock* const cbdlBox = Gtk::manage(new ToolParamBlock());

    cbdlBox->pack_start(*buttonBox1);

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
        cbdlBox->pack_start(*multiplier[i]);
    }

    Gtk::HSeparator *separator3 = Gtk::manage(new  Gtk::HSeparator());
    cbdlBox->pack_start(*separator3, Gtk::PACK_SHRINK, 2);

    chromacbdl->set_tooltip_text(M("TP_LOCALLAB_CHROMACB_TOOLTIP"));

    chromacbdl->setAdjusterListener(this);
    cbdlBox->pack_start(*chromacbdl);

    threshold->setAdjusterListener(this);
    cbdlBox->pack_start(*threshold);

    sensicb->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensicb->setAdjusterListener(this);
    cbdlBox->pack_start(*sensicb);

    sharradius->setAdjusterListener(this);

    sharamount->setAdjusterListener(this);

    shardamping->setAdjusterListener(this);

    shariter->setAdjusterListener(this);


    sensisha->set_tooltip_text(M("TP_LOCALLAB_SENSIS_TOOLTIP"));
    sensisha->setAdjusterListener(this);

    inverssha->set_active(false);
    inversshaConn  = inverssha->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversshaChanged));
    ToolParamBlock* const sharpBox = Gtk::manage(new ToolParamBlock());

    sharpBox->pack_start(*sharradius);
    sharpBox->pack_start(*sharamount);
    sharpBox->pack_start(*shardamping);
    sharpBox->pack_start(*shariter);
    sharpBox->pack_start(*sensisha);
    sharpBox->pack_start(*inverssha);


    noiselumf->setAdjusterListener(this);

    noiselumc->setAdjusterListener(this);

    noisechrof->setAdjusterListener(this);

    noisechroc->setAdjusterListener(this);
    ToolParamBlock* const denoisBox = Gtk::manage(new ToolParamBlock());

    denoisBox->pack_start(*noiselumf);
    denoisBox->pack_start(*noiselumc);
    denoisBox->pack_start(*noisechrof);
    denoisBox->pack_start(*noisechroc);

    neutrHBox1 = Gtk::manage(new Gtk::HBox());

    neutral1 = Gtk::manage(new Gtk::Button(M("TP_LOCALLAB_NEUTRAL")));
    RTImage *resetImg1 = Gtk::manage(new RTImage("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
    neutral1->set_image(*resetImg1);
    neutral1->set_tooltip_text(M("TP_LOCALLAB_NEUTRAL_TIP"));
    neutralconn1 = neutral1->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::neutral_pressed));
    neutral1->show();
    neutrHBox1->pack_start(*neutral1);
    pack_start(*neutrHBox1);

    superFrame->set_label_align(0.025, 0.5);
    //  Gtk::VBox *superVBox = Gtk::manage ( new Gtk::VBox());
    ToolParamBlock* const superBox = Gtk::manage(new ToolParamBlock());

    superFrame->set_label_widget(*curvactiv);
    ToolParamBlock* const colorBox = Gtk::manage(new ToolParamBlock());

    //  ToolParamBlock* const dustBox = Gtk::manage (new ToolParamBlock());
    dustMethod->append(M("TP_LOCALLAB_DSCOP"));
    dustMethod->append(M("TP_LOCALLAB_DSMOV"));
    dustMethod->append(M("TP_LOCALLAB_DSPAS"));
    dustMethod->set_active(0);
    dustMethodConn = dustMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::dustMethodChanged));
    dustMethod->set_tooltip_markup(M("TP_LOCALLAB_BLMETHOD_TOOLTIP"));

    superBox->pack_start(*lightness);
    superBox->pack_start(*contrast);
    superFrame->add(*superBox);
    colorBox->pack_start(*superFrame);

    colorBox->pack_start(*chroma);
    colorBox->pack_start(*sensi);

    dustFrame->set_label_align(0.025, 0.5);
//    dustBox->pack_start (*dustMethod);
//    dustBox->pack_start (*lastdust);

//    dustBox->pack_start (*cutpast);
//    dustBox->pack_start (*centerXbuf);
//    dustBox->pack_start (*centerYbuf);
//    dustBox->pack_start (*adjblur);
//    dustFrame->add (*dustBox);
//  colorBox->pack_start (*dustFrame);
    centerXbuf->hide();
    centerYbuf->hide();

    qualcurvbox->pack_start(*labqualcurv, Gtk::PACK_SHRINK, 4);
    qualcurvbox->pack_start(*qualitycurveMethod);

    colorBox->pack_start(*qualcurvbox);


    colorBox->pack_start(*llCurveEditorG, Gtk::PACK_SHRINK, 2);
    colorBox->pack_start(*invers);

    expcolor->add(*colorBox);
    expcolor->setLevel(2);
    pack_start(*expcolor);

    ToolParamBlock* const exposeBox = Gtk::manage(new ToolParamBlock());

    curveEditorG = new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVEEDITOR_TONES_LABEL"));
    curveEditorG->setCurveListener(this);

    shape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, ""));
    shape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_TONES_TOOLTIP"));
    curveEditorG->curveListComplete();


    exposeBox->pack_start(*expcomp);
    exposeBox->pack_start(*hlcompr);
    exposeBox->pack_start(*hlcomprthresh);
    exposeBox->pack_start(*black);
    exposeBox->pack_start(*shcompr);
    exposeBox->pack_start(*sensiex);
    exposeBox->pack_start(*curveEditorG);

    expexpose->add(*exposeBox);
    expexpose->setLevel(2);
    pack_start(*expexpose);

    ToolParamBlock* const vibranceBox = Gtk::manage(new ToolParamBlock());
    std::vector<GradientMilestone> milestonesvib;
    float R, G, B;
    // -0.1 rad < Hue < 1.6 rad
    Color::hsv2rgb01(0.92f, 0.45f, 0.6f, R, G, B);
    milestonesvib.push_back(GradientMilestone(0.0, double (R), double (G), double (B)));
    Color::hsv2rgb01(0.14056f, 0.45f, 0.6f, R, G, B);
    milestonesvib.push_back(GradientMilestone(1.0, double (R), double (G), double (B)));

    saturated = Gtk::manage(new Adjuster(M("TP_VIBRANCE_SATURATED"), -100., 100., 1., 0.));
    saturated->setAdjusterListener(this);
    saturated->set_sensitive(false);
    vibranceBox->pack_start(*saturated, Gtk::PACK_SHRINK, 0);

    pastels = Gtk::manage(new Adjuster(M("TP_VIBRANCE_PASTELS"), -100., 100., 1., 0.));
    pastels->setAdjusterListener(this);
    vibranceBox->pack_start(*pastels, Gtk::PACK_SHRINK, 0);

    psThreshold = Gtk::manage(new ThresholdAdjuster(M("TP_VIBRANCE_PSTHRESHOLD"), -100., 100., 0., M("TP_VIBRANCE_PSTHRESHOLD_WEIGTHING"), 0, 0., 100., 75., M("TP_VIBRANCE_PSTHRESHOLD_SATTHRESH"), 0, this, false));
    psThreshold->setAdjusterListener(this);
    psThreshold->set_tooltip_markup(M("TP_VIBRANCE_PSTHRESHOLD_TOOLTIP"));
    psThreshold->set_sensitive(false);
    vibranceBox->pack_start(*psThreshold, Gtk::PACK_SHRINK, 0);

    protectSkins = Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PROTECTSKINS")));
    protectSkins->set_active(true);
    vibranceBox->pack_start(*protectSkins, Gtk::PACK_SHRINK, 0);

    avoidColorShift = Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_AVOIDCOLORSHIFT")));
    avoidColorShift->set_active(true);
    vibranceBox->pack_start(*avoidColorShift, Gtk::PACK_SHRINK, 0);

    pastSatTog = Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PASTSATTOG")));
    pastSatTog->set_active(true);
    vibranceBox->pack_start(*pastSatTog, Gtk::PACK_SHRINK, 0);

    sensiv = Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19));
    sensiv->setAdjusterListener(this);

    vibranceBox->pack_start(*sensiv, Gtk::PACK_SHRINK, 0);

    curveEditorGG = new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_LABEL"));
    curveEditorGG->setCurveListener(this);

    skinTonesCurve = static_cast<DiagonalCurveEditor*>(curveEditorGG->addCurve(CT_Diagonal, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES")));
    skinTonesCurve->setTooltip(M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_TOOLTIP"));
    skinTonesCurve->setBottomBarBgGradient(milestonesvib);
    skinTonesCurve->setLeftBarBgGradient(milestonesvib);
    skinTonesCurve->setRangeLabels(
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE1"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE2"),
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE3"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE4")
    );
    skinTonesCurve->setRangeDefaultMilestones(0.1, 0.4, 0.85);
    curveEditorGG->curveListComplete();

    vibranceBox->pack_start(*curveEditorGG, Gtk::PACK_SHRINK, 4);

    pskinsconn = protectSkins->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::protectskins_toggled));
    ashiftconn = avoidColorShift->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::avoidcolorshift_toggled));
    pastsattogconn = pastSatTog->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::pastsattog_toggled));

    expvibrance->add(*vibranceBox);
    expvibrance->setLevel(2);
    pack_start(*expvibrance);




    ToolParamBlock* const blurrBox = Gtk::manage(new ToolParamBlock());
    blurMethod->append(M("TP_LOCALLAB_BLNORM"));
    blurMethod->append(M("TP_LOCALLAB_BLINV"));
    blurMethod->append(M("TP_LOCALLAB_BLSYM"));
    blurMethod->set_active(0);
    blurMethodConn = blurMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::blurMethodChanged));
    blurMethod->set_tooltip_markup(M("TP_LOCALLAB_BLMETHOD_TOOLTIP"));

    blurrBox->pack_start(*radius);
    blurrBox->pack_start(*strength);
    blurrBox->pack_start(*sensibn);
    blurrBox->pack_start(*blurMethod);


    blurrBox->pack_start(*activlum);

//    blurrBox->pack_start (*inversrad);
    expblur->add(*blurrBox);
    expblur->setLevel(2);
    pack_start(*expblur);
    ToolParamBlock* const tmBox = Gtk::manage(new ToolParamBlock());

    tmBox->pack_start(*stren);
    tmBox->pack_start(*gamma);
    tmBox->pack_start(*estop);
    tmBox->pack_start(*scaltm);
    tmBox->pack_start(*rewei);
    tmBox->pack_start(*sensitm);

    exptonemap->add(*tmBox);
    exptonemap->setLevel(2);
    pack_start(*exptonemap);

    ToolParamBlock* const retiBox = Gtk::manage(new ToolParamBlock());

    retiBox->pack_start(*retinexMethod);
    retiBox->pack_start(*str);
    retiBox->pack_start(*chrrt);
    retiBox->pack_start(*neigh);
    retiBox->pack_start(*vart);
    retiBox->pack_start(*sensih);
    retiBox->pack_start(*retrab);

    retiBox->pack_start(*LocalcurveEditorgainTrab, Gtk::PACK_SHRINK, 4);

    retiBox->pack_start(*LocalcurveEditorgainT, Gtk::PACK_SHRINK, 4);
    retiBox->pack_start(*inversret);

    expreti->add(*retiBox);
    expreti->setLevel(2);
    pack_start(*expreti);


    expsharp->add(*sharpBox);
    expsharp->setLevel(2);
    pack_start(*expsharp);

    expcbdl->add(*cbdlBox);
    expcbdl->setLevel(2);
    pack_start(*expcbdl);

    expdenoi->add(*denoisBox);
    expdenoi->setLevel(2);
    pack_start(*expdenoi);


//    pack_start (*transit);
    pack_start(*avoid); //keep avoid clor shift in case of

    neutrHBox = Gtk::manage(new Gtk::HBox());

    neutral = Gtk::manage(new Gtk::Button(M("TP_LOCALLAB_NEUTRAL")));
    RTImage *resetImg = Gtk::manage(new RTImage("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
    neutral->set_image(*resetImg);
    neutral->set_tooltip_text(M("TP_LOCALLAB_NEUTRAL_TIP"));
    neutralconn = neutral->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::neutral_pressed));
    neutral->show();
    neutrHBox->pack_start(*neutral);
    pack_start(*neutrHBox);


    // Instantiating the Editing geometry; positions will be initialized later
//   Line  *hLine, *vLine, *locYLine[2], *locXLine[2];
    Line  *locYLine[2], *locXLine[2];
    Circle *centerCircle;
//   Arcellipse *oneellipse;

    Beziers *onebeziers[4] = {};
    Beziers *twobeziers[4] = {};
    Beziers *thrbeziers[4] = {};
    Beziers *foubeziers[4] = {};
    float innw = 0.7f;
    // Visible geometry
    locXLine[0] = new Line();
    locXLine[0]->innerLineWidth = 2;
    locXLine[1] = new Line();
    locXLine[1]->innerLineWidth = 2;
    locXLine[0]->datum  = locXLine[1]->datum = Geometry::IMAGE;

    locYLine[0] = new Line();
    locYLine[0]->innerLineWidth = 2;
    locYLine[1] = new Line();
    locYLine[1]->innerLineWidth = 2;
    locYLine[0]->datum = locYLine[1]->datum = Geometry::IMAGE;

    centerCircle = new Circle();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = true;
    centerCircle->radius = circrad->getValue(); //19;
    centerCircle->filled = false;

    if (options.showdelimspot) {
        onebeziers[0] = new Beziers();
        onebeziers[0]->datum = Geometry::IMAGE;
        onebeziers[0]->innerLineWidth = innw;

        onebeziers[1] = new Beziers();
        onebeziers[1]->datum = Geometry::IMAGE;
        onebeziers[1]->innerLineWidth = innw;

        onebeziers[2] = new Beziers();
        onebeziers[2]->datum = Geometry::IMAGE;
        onebeziers[2]->innerLineWidth = innw;

        onebeziers[3] = new Beziers();
        onebeziers[3]->datum = Geometry::IMAGE;
        onebeziers[3]->innerLineWidth = innw;

        twobeziers[0] = new Beziers();
        twobeziers[0]->datum = Geometry::IMAGE;
        twobeziers[0]->innerLineWidth = innw;

        twobeziers[1] = new Beziers();
        twobeziers[1]->datum = Geometry::IMAGE;
        twobeziers[1]->innerLineWidth = innw;

        twobeziers[2] = new Beziers();
        twobeziers[2]->datum = Geometry::IMAGE;
        twobeziers[2]->innerLineWidth = innw;

        twobeziers[3] = new Beziers();
        twobeziers[3]->datum = Geometry::IMAGE;
        twobeziers[3]->innerLineWidth = innw;

        thrbeziers[0] = new Beziers();
        thrbeziers[0]->datum = Geometry::IMAGE;
        thrbeziers[0]->innerLineWidth = innw;

        thrbeziers[1] = new Beziers();
        thrbeziers[1]->datum = Geometry::IMAGE;
        thrbeziers[1]->innerLineWidth = innw;

        thrbeziers[2] = new Beziers();
        thrbeziers[2]->datum = Geometry::IMAGE;
        thrbeziers[2]->innerLineWidth = innw;

        thrbeziers[3] = new Beziers();
        thrbeziers[3]->datum = Geometry::IMAGE;
        thrbeziers[3]->innerLineWidth = innw;

        foubeziers[0] = new Beziers();
        foubeziers[0]->datum = Geometry::IMAGE;
        foubeziers[0]->innerLineWidth = innw;

        foubeziers[1] = new Beziers();
        foubeziers[1]->datum = Geometry::IMAGE;
        foubeziers[1]->innerLineWidth = innw;

        foubeziers[2] = new Beziers();
        foubeziers[2]->datum = Geometry::IMAGE;
        foubeziers[2]->innerLineWidth = innw;

        foubeziers[3] = new Beziers();
        foubeziers[3]->datum = Geometry::IMAGE;
        foubeziers[3]->innerLineWidth = innw;
    }

    // oneellipse->radiusInImageSpace = true;
    // oneellipse->radius = locX->getValue();
    // oneellipse->filled = false;

    EditSubscriber::visibleGeometry.push_back(locXLine[0]);
    EditSubscriber::visibleGeometry.push_back(locXLine[1]);
    EditSubscriber::visibleGeometry.push_back(locYLine[0]);
    EditSubscriber::visibleGeometry.push_back(locYLine[1]);
    EditSubscriber::visibleGeometry.push_back(centerCircle);

    if (options.showdelimspot) {
        EditSubscriber::visibleGeometry.push_back(onebeziers[0]);
        EditSubscriber::visibleGeometry.push_back(onebeziers[1]);
        EditSubscriber::visibleGeometry.push_back(onebeziers[2]);
        EditSubscriber::visibleGeometry.push_back(onebeziers[3]);
        EditSubscriber::visibleGeometry.push_back(twobeziers[0]);
        EditSubscriber::visibleGeometry.push_back(twobeziers[1]);
        EditSubscriber::visibleGeometry.push_back(twobeziers[2]);
        EditSubscriber::visibleGeometry.push_back(twobeziers[3]);
        EditSubscriber::visibleGeometry.push_back(thrbeziers[0]);
        EditSubscriber::visibleGeometry.push_back(thrbeziers[1]);
        EditSubscriber::visibleGeometry.push_back(thrbeziers[2]);
        EditSubscriber::visibleGeometry.push_back(thrbeziers[3]);
        EditSubscriber::visibleGeometry.push_back(foubeziers[0]);
        EditSubscriber::visibleGeometry.push_back(foubeziers[1]);
        EditSubscriber::visibleGeometry.push_back(foubeziers[2]);
        EditSubscriber::visibleGeometry.push_back(foubeziers[3]);
    }

    // MouseOver geometry
    locXLine[0] = new Line();
    locXLine[0]->innerLineWidth = 2;
    locXLine[1] = new Line();
    locXLine[1]->innerLineWidth = 2;
    locXLine[0]->datum  = locXLine[1]->datum = Geometry::IMAGE;

    locYLine[0] = new Line();
    locYLine[0]->innerLineWidth = 2;
    locYLine[1] = new Line();
    locYLine[1]->innerLineWidth = 2;
    locYLine[0]->datum = locYLine[1]->datum = Geometry::IMAGE;

    centerCircle = new Circle();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = true;
    centerCircle->radius = circrad->getValue();//19;
    centerCircle->filled = true;

    if (options.showdelimspot) {
        onebeziers[0]   = new Beziers();
        onebeziers[0]->datum = Geometry::IMAGE;
        onebeziers[0]->innerLineWidth = innw;

        onebeziers[1]   = new Beziers();
        onebeziers[1]->datum = Geometry::IMAGE;
        onebeziers[1]->innerLineWidth = innw;

        onebeziers[2]   = new Beziers();
        onebeziers[2]->datum = Geometry::IMAGE;
        onebeziers[2]->innerLineWidth = innw;

        onebeziers[3]   = new Beziers();
        onebeziers[3]->datum = Geometry::IMAGE;
        onebeziers[3]->innerLineWidth = innw;

        twobeziers[0] = new Beziers();
        twobeziers[0]->datum = Geometry::IMAGE;
        twobeziers[0]->innerLineWidth = innw;

        twobeziers[1] = new Beziers();
        twobeziers[1]->datum = Geometry::IMAGE;
        twobeziers[1]->innerLineWidth = innw;

        twobeziers[2] = new Beziers();
        twobeziers[2]->datum = Geometry::IMAGE;
        twobeziers[2]->innerLineWidth = innw;

        twobeziers[3] = new Beziers();
        twobeziers[3]->datum = Geometry::IMAGE;
        twobeziers[3]->innerLineWidth = innw;

        thrbeziers[0] = new Beziers();
        thrbeziers[0]->datum = Geometry::IMAGE;
        thrbeziers[0]->innerLineWidth = innw;

        thrbeziers[1] = new Beziers();
        thrbeziers[1]->datum = Geometry::IMAGE;
        thrbeziers[1]->innerLineWidth = innw;

        thrbeziers[2] = new Beziers();
        thrbeziers[2]->datum = Geometry::IMAGE;
        thrbeziers[2]->innerLineWidth = innw;

        thrbeziers[3] = new Beziers();
        thrbeziers[3]->datum = Geometry::IMAGE;
        thrbeziers[3]->innerLineWidth = innw;

        foubeziers[0] = new Beziers();
        foubeziers[0]->datum = Geometry::IMAGE;
        foubeziers[0]->innerLineWidth = innw;

        foubeziers[1] = new Beziers();
        foubeziers[1]->datum = Geometry::IMAGE;
        foubeziers[1]->innerLineWidth = innw;

        foubeziers[2] = new Beziers();
        foubeziers[2]->datum = Geometry::IMAGE;
        foubeziers[2]->innerLineWidth = innw;

        foubeziers[3] = new Beziers();
        foubeziers[3]->datum = Geometry::IMAGE;
        foubeziers[3]->innerLineWidth = innw;

    }

//   oneellipse->radiusInImageSpace = true;
//   oneellipse->radius = 10;//locX->getValue();
//    oneellipse->filled = false;

    EditSubscriber::mouseOverGeometry.push_back(locXLine[0]);
    EditSubscriber::mouseOverGeometry.push_back(locXLine[1]);

    EditSubscriber::mouseOverGeometry.push_back(locYLine[0]);
    EditSubscriber::mouseOverGeometry.push_back(locYLine[1]);

    EditSubscriber::mouseOverGeometry.push_back(centerCircle);

    if (options.showdelimspot) {
        EditSubscriber::mouseOverGeometry.push_back(onebeziers[0]);
        EditSubscriber::mouseOverGeometry.push_back(onebeziers[1]);
        EditSubscriber::mouseOverGeometry.push_back(onebeziers[2]);
        EditSubscriber::mouseOverGeometry.push_back(onebeziers[3]);
        EditSubscriber::mouseOverGeometry.push_back(twobeziers[0]);
        EditSubscriber::mouseOverGeometry.push_back(twobeziers[1]);
        EditSubscriber::mouseOverGeometry.push_back(twobeziers[2]);
        EditSubscriber::mouseOverGeometry.push_back(twobeziers[3]);
        EditSubscriber::mouseOverGeometry.push_back(thrbeziers[0]);
        EditSubscriber::mouseOverGeometry.push_back(thrbeziers[1]);
        EditSubscriber::mouseOverGeometry.push_back(thrbeziers[2]);
        EditSubscriber::mouseOverGeometry.push_back(thrbeziers[3]);
        EditSubscriber::mouseOverGeometry.push_back(foubeziers[0]);
        EditSubscriber::mouseOverGeometry.push_back(foubeziers[1]);
        EditSubscriber::mouseOverGeometry.push_back(foubeziers[2]);
        EditSubscriber::mouseOverGeometry.push_back(foubeziers[3]);
    }

    show_all();
}

Locallab::~Locallab()
{
    for (std::vector<Geometry*>::const_iterator i = visibleGeometry.begin(); i != visibleGeometry.end(); ++i) {
        delete *i;
    }

    for (std::vector<Geometry*>::const_iterator i = mouseOverGeometry.begin(); i != mouseOverGeometry.end(); ++i) {
        delete *i;
    }

    delete LocalcurveEditorgainT;
    delete LocalcurveEditorgainTrab;
    delete llCurveEditorG;

}
void Locallab::foldAllButMe(GdkEventButton* event, MyExpander *expander)
{
    if (event->button == 3) {
        expsettings->set_expanded(expsettings == expander);
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
    if (listener) {
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
    tpOpen.push_back(expsettings->get_expanded());
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
        expsettings->set_expanded(tpOpen.at(0));
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



void Locallab::neutral_pressed()
{
    Smethod->set_active(0);
    Exclumethod->set_active(0);
    locX->resetValue(false);
    locXL->resetValue(false);
    locY->resetValue(false);
    locYT->resetValue(false);
    centerX->resetValue(false);
    centerY->resetValue(false);
    circrad->resetValue(false);
    centerXbuf->resetValue(false);
    centerYbuf->resetValue(false);
    adjblur->resetValue(false);
    blurMethod->set_active(0);
    dustMethod->set_active(1);
    sensiexclu->resetValue(false);
    struc->resetValue(false);

    qualityMethod->set_active(0);
    qualitycurveMethod->set_active(0);
    thres->resetValue(false);
    proxi->resetValue(false);
    lightness->resetValue(false);
    chroma->resetValue(false);
    contrast->resetValue(false);
    sensi->resetValue(false);
    radius->resetValue(false);
    strength->resetValue(false);
    transit->resetValue(false);
    sensibn->resetValue(false);
    invers->set_active(false);
    cutpast->set_active(false);
    lastdust->set_active(false);
    curvactiv->set_active(false);
    inversrad->set_active(false);
    inversret->set_active(false);
    stren->resetValue(false);
    gamma->resetValue(false);
    estop->resetValue(false);
    scaltm->resetValue(false);
    rewei->resetValue(false);
    sensitm->resetValue(false);
    retinexMethod->set_active(2);
    str->resetValue(false);
    neigh->resetValue(false);
    vart->resetValue(false);
    chrrt->resetValue(false);
    sensih->resetValue(false);
    retrab->resetValue(false);
//    cTgainshape->reset();
//  cTgainshape->setCurve (creti);
    avoid->set_active(false);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->resetValue(false);
    }

    chromacbdl->resetValue(false);
    threshold->resetValue(false);
    sensicb->resetValue(false);
    sharradius->resetValue(false);
    sharamount->resetValue(false);
    shardamping->resetValue(false);
    shariter->resetValue(false);
    sensisha->resetValue(false);
    inverssha->set_active(false);
    noiselumf->resetValue(false);
    noiselumc->resetValue(false);
    noisechrof->resetValue(false);
    noisechroc->resetValue(false);


}


void Locallab::lumaneutralPressed()
{

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue(100);
        adjusterChanged(multiplier[i], 100);
    }
}


void Locallab::lumacontrastPlusPressed()
{

    for (int i = 0; i < 5; i++) {
        float inc = (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
        adjusterChanged(multiplier[i], multiplier[i]->getValue());
    }
}


void Locallab::lumacontrastMinusPressed()
{

    for (int i = 0; i < 5; i++) {
        float inc = - (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
        adjusterChanged(multiplier[i], multiplier[i]->getValue());
    }
}



void Locallab::autoOpenCurve()
{
    cTgainshape->openIfNonlinear();
    //  llshape->openIfNonlinear();
    //  LHshape->openIfNonlinear();

}


int localChangedUI(void* data)
{

    GThreadLock lock;
    (static_cast<Locallab*>(data))->localComputed_();

    return 0;
}

int localretChangedUI(void* data)
{

    GThreadLock lock;
    (static_cast<Locallab*>(data))->localretComputed_();

    return 0;
}

bool Locallab::localretComputed_()
{
    disableListener();

    //Reticurv
//update GUI and MIP specially for curve

    int *s_datc;
    s_datc = new int[70];
    int siz;
    //printf("nexts=%s\n", nextstr2.c_str());
    ImProcFunctions::strcurv_data(nextstr2, s_datc, siz);
    std::vector<double>   creti;

    for (int j = 0; j < siz; j++) {
        creti.push_back((double)(s_datc[j]) / 1000.);
    }

    delete [] s_datc;

    cTgainshape->setCurve(creti);

    int *s_datcl;
    s_datcl = new int[70];
    int sizl;
    ImProcFunctions::strcurv_data(nextll_str2, s_datcl, sizl);
    std::vector<double>   cll;

    for (int j = 0; j < sizl; j++) {
        cll.push_back((double)(s_datcl[j]) / 1000.);
    }

    delete [] s_datcl;

    llshape->setCurve(cll);


    int *s_datcc;
    s_datcc = new int[70];
    int sizc;
    ImProcFunctions::strcurv_data(nextcc_str2, s_datcc, sizc);
    std::vector<double>   ccc;

    for (int j = 0; j < sizc; j++) {
        ccc.push_back((double)(s_datcc[j]) / 1000.);
    }

    delete [] s_datcc;

    ccshape->setCurve(ccc);

    int *s_datch;
    s_datch = new int[70];
    int sizh;
    ImProcFunctions::strcurv_data(nextlh_str2, s_datch, sizh);
    std::vector<double>   clh;

    for (int j = 0; j < sizh; j++) {
        clh.push_back((double)(s_datch[j]) / 1000.);
    }

    delete [] s_datch;

    LHshape->setCurve(clh);


    int *s_datchh;
    s_datchh = new int[70];
    int sizhh;
    ImProcFunctions::strcurv_data(nexthh_str2, s_datchh, sizhh);
    std::vector<double>   chh;

    for (int j = 0; j < sizhh; j++) {
        chh.push_back((double)(s_datchh[j]) / 1000.);
    }

    delete [] s_datchh;

    HHshape->setCurve(chh);

    //skinTonesCurve
    int *s_datcsk;
    s_datcsk = new int[70];
    int sizsk;
    ImProcFunctions::strcurv_data(nextsk_str2, s_datcsk, sizsk);
    std::vector<double>   csk;

    for (int j = 0; j < sizsk; j++) {
        csk.push_back((double)(s_datcsk[j]) / 1000.);
    }

    delete [] s_datcsk;

    skinTonesCurve->setCurve(csk);

    //PSthreshold
    int sizps = 2;
    int s_datcps[sizps + 1];
    ImProcFunctions::strcurv_data(nextps_str2, s_datcps, sizps);
    psThreshold->setValue(s_datcps[0], s_datcps[1]);


    //exCurve
    int *s_datcex;
    s_datcex = new int[70];
    int sizex;
    ImProcFunctions::strcurv_data(nextex_str2, s_datcex, sizex);
    std::vector<double>   cex;

    for (int j = 0; j < sizex; j++) {
        cex.push_back((double)(s_datcex[j]) / 1000.);
    }

    delete [] s_datcex;

    shape->setCurve(cex);

    enableListener();

    //update all sliders by this strange process!
    if (anbspot->getValue() == 0) {
        anbspot->setValue(1);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged(anbspot, 1);

    } else if (anbspot->getValue() == 1) {
        anbspot->setValue(0);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged(anbspot, 0);

    }

    //update all curves
    std::vector<double>   cretirab;
    cretirab = cTgainshaperab->getCurve();

    if (cretirab.at(5) == 0.70) {
        cretirab.at(5) = 0.9;
        cTgainshaperab->setCurve(cretirab);

        curveChanged(cTgainshaperab);
    } else if (cretirab.at(5) == 0.90) {
        cretirab.at(5) = 0.7;
        cTgainshaperab->setCurve(cretirab);
        curveChanged(cTgainshaperab);

    }


//    printf("G2 anbspot=%i\n", anbspot->getValue());

    if (listener) { //for all sliders
        listener->panelChanged(Evlocallabanbspot, ""); //anbspot->getTextValue());
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabCTgainCurverab, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabCTgainCurve, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(Evlocallabllshape, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(Evlocallabccshape, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabLHshape, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabHHshape, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabSkinTonesCurve, M(""));
    }

    if (listener) {//for PSthreshold
        listener->panelChanged(EvlocallabPastSatThreshold, M(""));
    }

    if (listener) {//for PSthreshold
        listener->panelChanged(Evlocallabshape, M(""));
    }


    return false;

}

bool Locallab::localComputed_()
{
//update GUI and MIP
    disableListener();

    //size spot
    circrad->setValue(nextdatasp[2]);
    //center and cursor
    locX->setValue(nextdatasp[3]);
    locY->setValue(nextdatasp[4]);
    locYT->setValue(nextdatasp[5]);
    locXL->setValue(nextdatasp[6]);
    centerX->setValue(nextdatasp[7]);
    centerY->setValue(nextdatasp[8]);

    //sliders
    lightness->setValue(nextdatasp[9]);
    contrast->setValue(nextdatasp[10]);
    chroma->setValue(nextdatasp[11]);
    sensi->setValue(nextdatasp[12]);
    transit->setValue(nextdatasp[13]);

    //inverse
    if (nextdatasp[14] == 0) {
        invers->set_active(false);
    } else {
        invers->set_active(true);
    }

    //method cursor
    if (nextdatasp[15] == 0) {
        Smethod->set_active(0);
    } else if (nextdatasp[15] == 1) {
        Smethod->set_active(1);
    } else if (nextdatasp[15] == 2) {
        Smethod->set_active(2);
    } else if (nextdatasp[15] == 3) {
        Smethod->set_active(3);
    }

    //sliders blurr
    radius->setValue(nextdatasp[17]);
    strength->setValue(nextdatasp[18]);
    sensibn->setValue(nextdatasp[19]);

    //inverse
    if (nextdatasp[20] == 0) {
        inversrad->set_active(false);
    } else {
        inversrad->set_active(true);
    }

    //sliders retinex
    str->setValue(nextdatasp[21]);
    chrrt->setValue(nextdatasp[22]);
    neigh->setValue(nextdatasp[23]);
    vart->setValue(nextdatasp[24]);
    sensih->setValue(nextdatasp[25]);

    //inverse
    if (nextdatasp[26] == 0) {
        inversret->set_active(false);
    } else {
        inversret->set_active(true);
    }

    //method retinex
    if (nextdatasp[27] == 0) {
        retinexMethod->set_active(0);
    } else if (nextdatasp[27] == 1) {
        retinexMethod->set_active(1);
    } else if (nextdatasp[27] == 2) {
        retinexMethod->set_active(2);
    }

    //sharpening
    sharradius->setValue(nextdatasp[28]);
    sharamount->setValue(nextdatasp[29]);
    shardamping->setValue(nextdatasp[30]);
    shariter->setValue(nextdatasp[31]);
    sensisha->setValue(nextdatasp[32]);

    if (nextdatasp[33] == 0) {
        inverssha->set_active(false);
    } else {
        inverssha->set_active(true);
    }

    if (nextdatasp[34] == 0) {
        qualityMethod->set_active(0);
    } else if (nextdatasp[34] == 1) {
        qualityMethod->set_active(1);
    } else if (nextdatasp[34] == 2) {
        qualityMethod->set_active(2);
    }

    thres->setValue(nextdatasp[35]);
    proxi->setValue(nextdatasp[36]);

    //denoise
    noiselumf->setValue(nextdatasp[37]);
    noiselumc->setValue(nextdatasp[38]);
    noisechrof->setValue(nextdatasp[39]);
    noisechroc->setValue(nextdatasp[40]);

    //cbdl
    multiplier[0]->setValue(nextdatasp[41]);
    multiplier[1]->setValue(nextdatasp[42]);
    multiplier[2]->setValue(nextdatasp[43]);
    multiplier[3]->setValue(nextdatasp[44]);
    multiplier[4]->setValue(nextdatasp[45]);
    threshold->setValue(nextdatasp[46]);
    sensicb->setValue(nextdatasp[47]);

    //blur luma
    if (nextdatasp[48] == 0) {
        activlum->set_active(false);
    } else {
        activlum->set_active(true);
    }

//TM
    stren->setValue(nextdatasp[49]);
    gamma->setValue(nextdatasp[50]);
    estop->setValue(nextdatasp[51]);
    scaltm->setValue(nextdatasp[52]);
    rewei->setValue(nextdatasp[53]);
    sensitm->setValue(nextdatasp[54]);
    //  usleep(10000);

    //Reticurv
    retrab->setValue(nextdatasp[55]);

    //curvactiv
    if (nextdatasp[56] == 0) {
        curvactiv->set_active(false);
    } else {
        curvactiv->set_active(true);
    }

    if (nextdatasp[57] == 0) {
        qualitycurveMethod->set_active(0);
    } else if (nextdatasp[57] == 1) {
        qualitycurveMethod->set_active(1);
    } else if (nextdatasp[57] == 2) {
        qualitycurveMethod->set_active(2);
    }

    sensiv->setValue(nextdatasp[58]);
    pastels->setValue(nextdatasp[59]);
    saturated->setValue(nextdatasp[60]);

    //protectskin
    if (nextdatasp[61] == 0) {
        protectSkins->set_active(false);
    } else {
        protectSkins->set_active(true);
    }

    //avoidColorShift
    if (nextdatasp[62] == 0) {
        avoidColorShift->set_active(false);
    } else {
        avoidColorShift->set_active(true);
    }

    //pastSatTog
    if (nextdatasp[63] == 0) {
        pastSatTog->set_active(false);
        adjusterChanged(pastels, pastels->getValue());
        adjusterChanged(saturated, saturated->getValue());

    } else {
        pastSatTog->set_active(true);
        adjusterChanged(pastels, pastels->getValue());
        adjusterChanged(saturated, saturated->getValue());

    }

    expcomp->setValue(nextdatasp[64]);
    black->setValue(nextdatasp[65]);
    hlcompr->setValue(nextdatasp[66]);
    hlcomprthresh->setValue(nextdatasp[67]);
    shcompr->setValue(nextdatasp[68]);
    sensiex->setValue(nextdatasp[69]);

    centerXbuf->setValue(nextdatasp[70]);
    centerYbuf->setValue(nextdatasp[71]);
    adjblur->setValue(nextdatasp[72]);

    //protectskin
    if (nextdatasp[73] == 0) {
        cutpast->set_active(false);
    } else {
        cutpast->set_active(true);
    }

    chromacbdl->setValue(nextdatasp[74]);

    if (nextdatasp[75] == 0) {
        lastdust->set_active(false);
    } else {
        lastdust->set_active(true);
    }

    if (nextdatasp[76] == 0) {
        blurMethod->set_active(0);
    } else if (nextdatasp[76] == 1) {
        blurMethod->set_active(1);
    } else if (nextdatasp[76] == 2) {
        blurMethod->set_active(2);
    }

    if (nextdatasp[77] == 0) {
        dustMethod->set_active(0);
    } else if (nextdatasp[77] == 1) {
        dustMethod->set_active(1);
    } else if (nextdatasp[77] == 2) {
        dustMethod->set_active(2);
    }

    if (nextdatasp[78] == 0) {
        Exclumethod->set_active(0);
    } else if (nextdatasp[78] == 1) {
        Exclumethod->set_active(1);
    }

    sensiexclu->setValue(nextdatasp[79]);
    struc->setValue(nextdatasp[80]);

    double intermed = 0.01 * (double) nextdatasp[81];
    hueref->setValue(intermed);
    chromaref->setValue(nextdatasp[82]);
    lumaref->setValue(nextdatasp[83]);
    sobelref->setValue(nextdatasp[84]);

    int *s_datc;
    s_datc = new int[70];
    int siz;
    ImProcFunctions::strcurv_data(nextstr, s_datc, siz);


    std::vector<double>   creti;

    for (int j = 0; j < siz; j++) {
        creti.push_back((double)(s_datc[j]) / 1000.);
    }

    delete [] s_datc;

    cTgainshape->setCurve(creti);

    //LLcurv
    int *s_datcl;
    s_datcl = new int[70];
    int sizl;
    ImProcFunctions::strcurv_data(nextll_str, s_datcl, sizl);


    std::vector<double>   cll;

    for (int j = 0; j < sizl; j++) {
        cll.push_back((double)(s_datcl[j]) / 1000.);
    }

    delete [] s_datcl;
    llshape->setCurve(cll);

    //CCcurv
    int *s_datcc;
    s_datcc = new int[70];
    int sizc;
    ImProcFunctions::strcurv_data(nextcc_str, s_datcc, sizc);


    std::vector<double>   ccc;

    for (int j = 0; j < sizc; j++) {
        ccc.push_back((double)(s_datcc[j]) / 1000.);
    }

    delete [] s_datcc;
    ccshape->setCurve(ccc);


    //LHcurv
    int *s_datch;
    s_datch = new int[70];
    int sizh;
    ImProcFunctions::strcurv_data(nextlh_str, s_datch, sizh);


    std::vector<double>   clh;

    for (int j = 0; j < sizh; j++) {
        clh.push_back((double)(s_datch[j]) / 1000.);
    }

    delete [] s_datch;
    LHshape->setCurve(clh);

    //HHcurv
    int *s_datchh;
    s_datchh = new int[70];
    int sizhh;
    ImProcFunctions::strcurv_data(nexthh_str, s_datchh, sizhh);


    std::vector<double>   chh;

    for (int j = 0; j < sizhh; j++) {
        chh.push_back((double)(s_datchh[j]) / 1000.);
    }

    delete [] s_datchh;
    HHshape->setCurve(chh);

    //skinTonesCurve
    int *s_datcsk;
    s_datcsk = new int[70];
    int sizsk;
    ImProcFunctions::strcurv_data(nextsk_str, s_datcsk, sizsk);
    std::vector<double>   csk;

    for (int j = 0; j < sizsk; j++) {
        csk.push_back((double)(s_datcsk[j]) / 1000.);
    }

    delete [] s_datcsk;

    skinTonesCurve->setCurve(csk);

    //PSthreshold
    int sizps = 2;
    int s_datcps[sizps + 1];
    ImProcFunctions::strcurv_data(nextps_str, s_datcps, sizps);
    psThreshold->setValue(s_datcps[0], s_datcps[1]);


    //exCurve
    int *s_datcex;
    s_datcex = new int[70];
    int sizex;
    ImProcFunctions::strcurv_data(nextex_str, s_datcex, sizex);
    std::vector<double>   cex;

    for (int j = 0; j < sizex; j++) {
        cex.push_back((double)(s_datcex[j]) / 1000.);
    }

    delete [] s_datcex;
    shape->setCurve(cex);


    enableListener();

    //update all sliders by this strange process!
    if (anbspot->getValue() == 0) {
        anbspot->setValue(1);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged(anbspot, 1);

    } else if (anbspot->getValue() == 1) {
        anbspot->setValue(0);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged(anbspot, 0);

    }


    //update all curves
    std::vector<double>   cretirab;
    cretirab = cTgainshaperab->getCurve();

    if (cretirab.at(5) == 0.70) {
        cretirab.at(5) = 0.9;
        cTgainshaperab->setCurve(cretirab);

        curveChanged(cTgainshaperab);
    } else if (cretirab.at(5) == 0.90) {
        cretirab.at(5) = 0.7;
        cTgainshaperab->setCurve(cretirab);
        curveChanged(cTgainshaperab);

    }

    //

//   printf("G1 maj anbspot=%i  cretirab=%f\n", anbspot->getValue(), cretirab.at(5));


    //add events for each cases
    if (listener) { //for all sliders
        listener->panelChanged(Evlocallabanbspot, ""); //anbspot->getTextValue());
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabCTgainCurverab, M(""));
    }

    if (listener) {//for inverse color
        listener->panelChanged(Evlocallabinvers, M("GENERAL_ENABLED"));
    }

    if (listener) {//for cutpast
        //   listener->panelChanged (Evlocallabcutpast, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for cutpast
        //    listener->panelChanged (Evlocallablastdust, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for blur method
        listener->panelChanged(EvlocallabblurMethod, blurMethod->get_active_text());
    }

    if (listener) {//for dust method
        //    listener->panelChanged (EvlocallabdustMethod, dustMethod->get_active_text ());
    }

    if (listener) {//for Exclude method
        listener->panelChanged(Evlocallabexclumethod, Exclumethod->get_active_text());
    }

    if (listener) {//for curvactiv
        listener->panelChanged(Evlocallabcurvactiv, M("GENERAL_ENABLED"));
    }

    if (listener) {//for inverse blurr
        listener->panelChanged(Evlocallabinversrad, M("GENERAL_ENABLED"));
    }

    if (listener) {//for quality method
        listener->panelChanged(EvlocallabqualityMethod, qualityMethod->get_active_text());
    }

    if (listener) {//for quality method
        listener->panelChanged(EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text());
    }

    if (listener) {//for inverse retinex
        listener->panelChanged(Evlocallabinversret, M("GENERAL_ENABLED"));
    }

    if (listener) {//for inverse sharpen
        listener->panelChanged(Evlocallabinverssha, M("GENERAL_ENABLED"));
    }

    if (listener) {//for Smethod : position of mouse cursor
        listener->panelChanged(EvlocallabSmet, Smethod->get_active_text());
    }

    if (listener) {//for retinex method
        listener->panelChanged(EvlocallabretinexMethod, retinexMethod->get_active_text());
    }

    if (listener) {//for curve reti
        listener->panelChanged(EvlocallabCTgainCurve, M(""));
    }

    if (listener) {//for curve LL
        listener->panelChanged(Evlocallabllshape, M(""));
    }

    if (listener) {//for curve LH
        listener->panelChanged(EvlocallabLHshape, M(""));
    }

    if (listener) {//for curve LH
        listener->panelChanged(Evlocallabccshape, M(""));
    }

    if (listener) {//for curve LH
        listener->panelChanged(EvlocallabHHshape, M(""));
    }

    if (listener) {//for curve Skin
        listener->panelChanged(EvlocallabSkinTonesCurve, M(""));
    }

    if (listener) {//for PSthreshold
        listener->panelChanged(EvlocallabPastSatThreshold, M(""));
    }

    //for checkbox
    if (listener) {//f
        listener->panelChanged(EvlocallabProtectSkins, M(""));
    }

    if (listener) {//fo
        listener->panelChanged(EvlocallabAvoidColorShift, M(""));
    }

    if (listener) {//for
        listener->panelChanged(EvlocallabPastSatTog, M(""));
    }

    if (listener) {//for
        listener->panelChanged(Evlocallabshape, M(""));
    }

    return false;
}

void Locallab::localChanged(int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, std::string hh_str, std::string sk_str, std::string ps_str, std::string ex_str, int sp, int maxdat)
{
    for (int i = 2; i < 85; i++) {
        nextdatasp[i] = datasp[i][sp];
    }

    nextstr = datastr;
    nextll_str = ll_str;
    nextlh_str = lh_str;
    nextcc_str = cc_str;
    nexthh_str = hh_str;
    nextsk_str = sk_str;
    nextps_str = ps_str;
    nextex_str = ex_str;

    nextlength = maxdat;
    g_idle_add(localChangedUI, this);
}

void Locallab::localretChanged(int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, std::string hh_str, std::string sk_str, std::string ps_str, std::string ex_str, int sp, int maxdat)
{
    nextlength = maxdat;
    nextstr2 = datastr;
    nextll_str2 = ll_str;
    nextlh_str2 = lh_str;
    nextcc_str2 = cc_str;
    nexthh_str2 = hh_str;
    nextsk_str2 = sk_str;
    nextps_str2 = ps_str;
    nextex_str2 = ex_str;

    g_idle_add(localretChangedUI, this);
}


void Locallab::read(const ProcParams* pp, const ParamsEdited* pedited)
{
    anbspot->hide();
    hueref->hide();
    chromaref->hide();
    lumaref->hide();
    sobelref->hide();
    centerXbuf->hide();
    centerYbuf->hide();

    disableListener();
    enablecolorConn.block(true);
    enablevibranceConn.block(true);
    enableblurConn.block(true);
    enabletonemapConn.block(true);
    enableretiConn.block(true);
    enablesharpConn.block(true);
    enablecbdlConn.block(true);
    enabledenoiConn.block(true);


    if (pedited) {
        degree->setEditedState(pedited->locallab.degree ? Edited : UnEdited);
        locY->setEditedState(pedited->locallab.locY ? Edited : UnEdited);
        locX->setEditedState(pedited->locallab.locX ? Edited : UnEdited);
        locYT->setEditedState(pedited->locallab.locYT ? Edited : UnEdited);
        locXL->setEditedState(pedited->locallab.locXL ? Edited : UnEdited);
        centerX->setEditedState(pedited->locallab.centerX ? Edited : UnEdited);
        centerY->setEditedState(pedited->locallab.centerY ? Edited : UnEdited);
        circrad->setEditedState(pedited->locallab.circrad ? Edited : UnEdited);
        centerXbuf->setEditedState(pedited->locallab.centerXbuf ? Edited : UnEdited);
        centerYbuf->setEditedState(pedited->locallab.centerYbuf ? Edited : UnEdited);
        adjblur->setEditedState(pedited->locallab.adjblur ? Edited : UnEdited);
        thres->setEditedState(pedited->locallab.thres ? Edited : UnEdited);
        proxi->setEditedState(pedited->locallab.proxi ? Edited : UnEdited);
        lightness->setEditedState(pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setEditedState(pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setEditedState(pedited->locallab.chroma ? Edited : UnEdited);
        expcomp->setEditedState(pedited->locallab.expcomp ? Edited : UnEdited);
        hlcompr->setEditedState(pedited->locallab.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setEditedState(pedited->locallab.hlcomprthresh ? Edited : UnEdited);
        black->setEditedState(pedited->locallab.black ? Edited : UnEdited);
        shcompr->setEditedState(pedited->locallab.shcompr ? Edited : UnEdited);
        sensiex->setEditedState(pedited->locallab.sensiex ? Edited : UnEdited);

        sharradius->setEditedState(pedited->locallab.sharradius ? Edited : UnEdited);
        sharamount->setEditedState(pedited->locallab.sharamount ? Edited : UnEdited);
        shardamping->setEditedState(pedited->locallab.shardamping ? Edited : UnEdited);
        shariter->setEditedState(pedited->locallab.shariter ? Edited : UnEdited);
        sensisha->setEditedState(pedited->locallab.sensisha ? Edited : UnEdited);
        noiselumf->setEditedState(pedited->locallab.noiselumf ? Edited : UnEdited);
        noiselumc->setEditedState(pedited->locallab.noiselumc ? Edited : UnEdited);
        noisechrof->setEditedState(pedited->locallab.noisechrof ? Edited : UnEdited);
        noisechroc->setEditedState(pedited->locallab.noisechroc ? Edited : UnEdited);

        pastels->setEditedState(pedited->locallab.pastels ? Edited : UnEdited);
        saturated->setEditedState(pedited->locallab.saturated ? Edited : UnEdited);
        psThreshold->setEditedState(pedited->locallab.psthreshold ? Edited : UnEdited);
        protectSkins->set_inconsistent(!pedited->locallab.protectskins);
        avoidColorShift->set_inconsistent(!pedited->locallab.avoidcolorshift);
        pastSatTog->set_inconsistent(!pedited->locallab.pastsattog);
        skinTonesCurve->setUnChanged(!pedited->locallab.skintonescurve);
        sensiv->setEditedState(pedited->locallab.sensiv ? Edited : UnEdited);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setEditedState(pedited->locallab.mult[i] ? Edited : UnEdited);
        }

        threshold->setEditedState(pedited->locallab.threshold ? Edited : UnEdited);
        chromacbdl->setEditedState(pedited->locallab.chromacbdl ? Edited : UnEdited);
        sensiexclu->setEditedState(pedited->locallab.sensiexclu ? Edited : UnEdited);
        struc->setEditedState(pedited->locallab.struc ? Edited : UnEdited);

        sensi->setEditedState(pedited->locallab.sensi ? Edited : UnEdited);
        sensih->setEditedState(pedited->locallab.sensih ? Edited : UnEdited);
        retrab->setEditedState(pedited->locallab.retrab ? Edited : UnEdited);
        sensicb->setEditedState(pedited->locallab.sensicb ? Edited : UnEdited);
        sensibn->setEditedState(pedited->locallab.sensibn ? Edited : UnEdited);
        sensitm->setEditedState(pedited->locallab.sensitm ? Edited : UnEdited);
        radius->setEditedState(pedited->locallab.radius ? Edited : UnEdited);
        strength->setEditedState(pedited->locallab.strength ? Edited : UnEdited);
        stren->setEditedState(pedited->locallab.stren ? Edited : UnEdited);
        gamma->setEditedState(pedited->locallab.gamma ? Edited : UnEdited);
        estop->setEditedState(pedited->locallab.estop ? Edited : UnEdited);
        scaltm->setEditedState(pedited->locallab.scaltm ? Edited : UnEdited);
        rewei->setEditedState(pedited->locallab.rewei ? Edited : UnEdited);

        nbspot->setEditedState(pedited->locallab.nbspot ? Edited : UnEdited);
        anbspot->setEditedState(pedited->locallab.anbspot ? Edited : UnEdited);
        hueref->setEditedState(pedited->locallab.hueref ? Edited : UnEdited);
        chromaref->setEditedState(pedited->locallab.chromaref ? Edited : UnEdited);
        lumaref->setEditedState(pedited->locallab.lumaref ? Edited : UnEdited);
        sobelref->setEditedState(pedited->locallab.sobelref ? Edited : UnEdited);
        transit->setEditedState(pedited->locallab.transit ? Edited : UnEdited);
        str->setEditedState(pedited->locallab.str ? Edited : UnEdited);
        neigh->setEditedState(pedited->locallab.neigh ? Edited : UnEdited);
        vart->setEditedState(pedited->locallab.vart ? Edited : UnEdited);
        chrrt->setEditedState(pedited->locallab.chrrt ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->locallab.enabled);
        avoid->set_inconsistent(multiImage && !pedited->locallab.avoid);
        activlum->set_inconsistent(multiImage && !pedited->locallab.activlum);
        invers->set_inconsistent(multiImage && !pedited->locallab.invers);
        cutpast->set_inconsistent(multiImage && !pedited->locallab.cutpast);
        lastdust->set_inconsistent(multiImage && !pedited->locallab.lastdust);
        curvactiv->set_inconsistent(multiImage && !pedited->locallab.curvactiv);
        inversrad->set_inconsistent(multiImage && !pedited->locallab.inversrad);
        inverssha->set_inconsistent(multiImage && !pedited->locallab.inverssha);
        cTgainshape->setUnChanged(!pedited->locallab.localTgaincurve);
        llshape->setUnChanged(!pedited->locallab.llcurve);
        ccshape->setUnChanged(!pedited->locallab.cccurve);
        LHshape->setUnChanged(!pedited->locallab.LHcurve);
        HHshape->setUnChanged(!pedited->locallab.HHcurve);
        shape->setUnChanged(!pedited->locallab.excurve);
        inversret->set_inconsistent(multiImage && !pedited->locallab.inversret);
        cTgainshaperab->setUnChanged(!pedited->locallab.localTgaincurverab);
        expcolor->set_inconsistent(!pedited->locallab.expcolor);
        expexpose->set_inconsistent(!pedited->locallab.expexpose);
        expvibrance->set_inconsistent(!pedited->locallab.expvibrance);
        expblur->set_inconsistent(!pedited->locallab.expblur);
        exptonemap->set_inconsistent(!pedited->locallab.exptonemap);
        expreti->set_inconsistent(!pedited->locallab.expreti);
        expsharp->set_inconsistent(!pedited->locallab.expsharp);
        expcbdl->set_inconsistent(!pedited->locallab.expcbdl);
        expdenoi->set_inconsistent(!pedited->locallab.expdenoi);

        if (!pedited->locallab.Smethod) {
            Smethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.Exclumethod) {
            Exclumethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.retinexMethod) {
            retinexMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.blurMethod) {
            blurMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.dustMethod) {
            dustMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.qualityMethod) {
            qualityMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.qualitycurveMethod) {
            qualitycurveMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

    }

    setEnabled(pp->locallab.enabled);

    Smethodconn.block(true);
    Exclumethodconn.block(true);
    retinexMethodConn.block(true);
    qualityMethodConn.block(true);
    qualitycurveMethodConn.block(true);
    blurMethodConn.block(true);
    dustMethodConn.block(true);

    avoidConn.block(true);
    avoid->set_active(pp->locallab.avoid);
    avoidConn.block(false);
    activlumConn.block(true);
    activlum->set_active(pp->locallab.activlum);
    activlumConn.block(false);
    inversConn.block(true);
    invers->set_active(pp->locallab.invers);
    inversConn.block(false);
    cutpastConn.block(true);
    cutpast->set_active(pp->locallab.cutpast);
    cutpastConn.block(false);
    lastdustConn.block(true);
    lastdust->set_active(pp->locallab.lastdust);
    lastdustConn.block(false);

    curvactivConn.block(true);
    curvactiv->set_active(pp->locallab.curvactiv);
    curvactivConn.block(false);
    inversradConn.block(true);
    inversrad->set_active(pp->locallab.inversrad);
    inversradConn.block(false);
    inversretConn.block(true);
    inversret->set_active(pp->locallab.inversret);
    inversretConn.block(false);
    inversshaConn.block(true);
    inverssha->set_active(pp->locallab.inverssha);
    inversshaConn.block(false);

    degree->setValue(pp->locallab.degree);
    locY->setValue(pp->locallab.locY);
    locX->setValue(pp->locallab.locX);
    locYT->setValue(pp->locallab.locYT);
    locXL->setValue(pp->locallab.locXL);
    centerX->setValue(pp->locallab.centerX);
    centerY->setValue(pp->locallab.centerY);
    circrad->setValue(pp->locallab.circrad);
    centerXbuf->setValue(pp->locallab.centerXbuf);
    centerYbuf->setValue(pp->locallab.centerYbuf);
    adjblur->setValue(pp->locallab.adjblur);
    thres->setValue(pp->locallab.thres);
    proxi->setValue(pp->locallab.proxi);
    lightness->setValue(pp->locallab.lightness);
    contrast->setValue(pp->locallab.contrast);
    chroma->setValue(pp->locallab.chroma);
    expcomp->setValue(pp->locallab.expcomp);
    hlcompr->setValue(pp->locallab.hlcompr);
    hlcomprthresh->setValue(pp->locallab.hlcomprthresh);
    black->setValue(pp->locallab.black);
    shcompr->setValue(pp->locallab.shcompr);
    shcompr->set_sensitive(!((int)black->getValue() == 0));     //at black=0 shcompr value has no effect
    sensiexclu->setValue(pp->locallab.sensiexclu);
    struc->setValue(pp->locallab.struc);

    sharradius->setValue(pp->locallab.sharradius);
    sharamount->setValue(pp->locallab.sharamount);
    shardamping->setValue(pp->locallab.shardamping);
    shariter->setValue(pp->locallab.shariter);
    sensisha->setValue(pp->locallab.sensisha);
    sensi->setValue(pp->locallab.sensi);
    sensiex->setValue(pp->locallab.sensiex);
    sensih->setValue(pp->locallab.sensih);
    retrab->setValue(pp->locallab.retrab);
    sensicb->setValue(pp->locallab.sensicb);
    sensibn->setValue(pp->locallab.sensibn);
    sensitm->setValue(pp->locallab.sensitm);
    transit->setValue(pp->locallab.transit);
    radius->setValue(pp->locallab.radius);
    strength->setValue(pp->locallab.strength);
    stren->setValue(pp->locallab.stren);
    gamma->setValue(pp->locallab.gamma);
    estop->setValue(pp->locallab.estop);
    scaltm->setValue(pp->locallab.scaltm);
    rewei->setValue(pp->locallab.rewei);
    str->setValue(pp->locallab.str);
    neigh->setValue(pp->locallab.neigh);
    nbspot->setValue(pp->locallab.nbspot);
    anbspot->setValue(pp->locallab.anbspot);
    hueref->setValue(pp->locallab.hueref);
    chromaref->setValue(pp->locallab.chromaref);
    lumaref->setValue(pp->locallab.lumaref);
    sobelref->setValue(pp->locallab.sobelref);
    vart->setValue(pp->locallab.vart);
    chrrt->setValue(pp->locallab.chrrt);
    cTgainshape->setCurve(pp->locallab.localTgaincurve);
    cTgainshaperab->setCurve(pp->locallab.localTgaincurverab);
    llshape->setCurve(pp->locallab.llcurve);
    ccshape->setCurve(pp->locallab.cccurve);
    LHshape->setCurve(pp->locallab.LHcurve);
    shape->setCurve(pp->locallab.excurve);
    HHshape->setCurve(pp->locallab.HHcurve);
    lastactivlum = pp->locallab.activlum;
    lastanbspot = pp->locallab.anbspot;
    noiselumf->setValue(pp->locallab.noiselumf);
    noiselumc->setValue(pp->locallab.noiselumc);
    noisechrof->setValue(pp->locallab.noisechrof);
    noisechroc->setValue(pp->locallab.noisechroc);
    expcolor->setEnabled(pp->locallab.expcolor);
    expexpose->setEnabled(pp->locallab.expexpose);
    expvibrance->setEnabled(pp->locallab.expvibrance);
    sensiv->setValue(pp->locallab.sensiv);

    pskinsconn.block(true);
    protectSkins->set_active(pp->locallab.protectskins);
    pskinsconn.block(false);
    lastProtectSkins = pp->locallab.protectskins;

    ashiftconn.block(true);
    avoidColorShift->set_active(pp->locallab.avoidcolorshift);
    ashiftconn.block(false);
    lastAvoidColorShift = pp->locallab.avoidcolorshift;

    pastsattogconn.block(true);
    pastSatTog->set_active(pp->locallab.pastsattog);
    pastsattogconn.block(false);
    lastPastSatTog = pp->locallab.pastsattog;

    pastels->setValue(pp->locallab.pastels);
    psThreshold->setValue<int> (pp->locallab.psthreshold);

    if (lastPastSatTog) {
        // Link both slider, so we set saturated and psThresholds unsensitive
        psThreshold->set_sensitive(false);
        saturated->set_sensitive(false);
        saturated->setValue(pp->locallab.pastels);     // Pastels and Saturated are linked
    } else {
        // Separate sliders, so we set saturated and psThresholds sensitive again
        psThreshold->set_sensitive(true);
        saturated->set_sensitive(true);
        saturated->setValue(pp->locallab.saturated);   // Pastels and Saturated are separate
    }

    skinTonesCurve->setCurve(pp->locallab.skintonescurve);



    expblur->setEnabled(pp->locallab.expblur);
    exptonemap->setEnabled(pp->locallab.exptonemap);
    expreti->setEnabled(pp->locallab.expreti);
    expsharp->setEnabled(pp->locallab.expsharp);
    expcbdl->setEnabled(pp->locallab.expcbdl);
    expdenoi->setEnabled(pp->locallab.expdenoi);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue(pp->locallab.mult[i]);
    }

    threshold->setValue(pp->locallab.threshold);
    chromacbdl->setValue(pp->locallab.chromacbdl);

    lastavoid = pp->locallab.avoid;
    lastinvers = pp->locallab.invers;
    lastcutpast = pp->locallab.cutpast;
    lastlastdust = pp->locallab.lastdust;
    lastcurvactiv = pp->locallab.curvactiv;
    lastinversrad = pp->locallab.inversrad;
    lastinversret = pp->locallab.inversret;
    lastinverssha = pp->locallab.inverssha;
    activlumChanged();
    inversChanged();
    cutpastChanged();
    lastdustChanged();
    curvactivChanged();
    inversradChanged();
    inversretChanged();
    inversshaChanged();

    updateGeometry(pp->locallab.centerX, pp->locallab.centerY, pp->locallab.circrad, pp->locallab.locY, pp->locallab.degree,  pp->locallab.locX, pp->locallab.locYT, pp->locallab.locXL);

    if (pp->locallab.Smethod == "IND") {
        Smethod->set_active(0);
    } else if (pp->locallab.Smethod == "SYM") {
        Smethod->set_active(1);
    } else if (pp->locallab.Smethod == "INDSL") {
        Smethod->set_active(2);
    } else if (pp->locallab.Smethod == "SYMSL") {
        Smethod->set_active(3);
    }

    SmethodChanged();
    Smethodconn.block(false);

    if (pp->locallab.Exclumethod == "norm") {
        Exclumethod->set_active(0);
    } else if (pp->locallab.Exclumethod == "exc") {
        Exclumethod->set_active(1);
    }

    ExclumethodChanged();
    Exclumethodconn.block(false);

    if (pp->locallab.retinexMethod == "low") {
        retinexMethod->set_active(0);
    } else if (pp->locallab.retinexMethod == "uni") {
        retinexMethod->set_active(1);
    } else if (pp->locallab.retinexMethod == "high") {
        retinexMethod->set_active(2);
    }

    retinexMethodChanged();
    retinexMethodConn.block(false);

    if (pp->locallab.blurMethod == "norm") {
        blurMethod->set_active(0);
    } else if (pp->locallab.blurMethod == "inv") {
        blurMethod->set_active(1);
    } else if (pp->locallab.blurMethod == "sym") {
        blurMethod->set_active(2);
    }

    blurMethodChanged();
    blurMethodConn.block(false);

    if (pp->locallab.dustMethod == "cop") {
        dustMethod->set_active(0);
    } else if (pp->locallab.dustMethod == "mov") {
        dustMethod->set_active(1);
    } else if (pp->locallab.dustMethod == "pas") {
        dustMethod->set_active(2);
    }

    dustMethodChanged();
    dustMethodConn.block(false);

    if (pp->locallab.qualityMethod == "std") {
        qualityMethod->set_active(0);
    } else if (pp->locallab.qualityMethod == "enh") {
        qualityMethod->set_active(1);
    } else if (pp->locallab.qualityMethod == "enhden") {
        qualityMethod->set_active(2);
    }

    qualityMethodChanged();
    qualityMethodConn.block(false);

    if (pp->locallab.qualitycurveMethod == "none") {
        qualitycurveMethod->set_active(0);
    } else if (pp->locallab.qualitycurveMethod == "std") {
        qualitycurveMethod->set_active(1);
    } else if (pp->locallab.qualitycurveMethod == "enh") {
        qualitycurveMethod->set_active(2);
    }

    qualitycurveMethodChanged();
    qualitycurveMethodConn.block(false);

    anbspot->hide();
    hueref->hide();
    chromaref->hide();
    lumaref->hide();
    sobelref->hide();
    centerXbuf->hide();
    centerYbuf->hide();

    if (pp->locallab.Smethod == "SYM" || pp->locallab.Smethod == "SYMSL") {
        locXL->setValue(locX->getValue());
        locYT->setValue(locY->getValue());
    } else if (pp->locallab.Smethod == "LOC") {
        locXL->setValue(locX->getValue());
        locYT->setValue(locX->getValue());
        locY->setValue(locX->getValue());
    } else if (pp->locallab.Smethod == "INDSL" || pp->locallab.Smethod == "IND") {
        locX->setValue(pp->locallab.locX);
        locY->setValue(pp->locallab.locY);
        locXL->setValue(pp->locallab.locXL);
        locYT->setValue(pp->locallab.locYT);

    }

    enablecolorConn.block(false);
    enablevibranceConn.block(false);
    enableblurConn.block(false);
    enabletonemapConn.block(false);
    enableretiConn.block(false);
    enablesharpConn.block(false);
    enablecbdlConn.block(false);
    enabledenoiConn.block(false);

    enableListener();
}

void Locallab::updateGeometry(const int centerX_, const int centerY_, const int circrad_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth, const int fullHeight)
{
    EditDataProvider* dataProvider = getEditProvider();


    if (!dataProvider) {
        return;
    }

    int imW = 0;
    int imH = 0;

    if (fullWidth != -1 && fullHeight != -1) {
        imW = fullWidth;
        imH = fullHeight;
    } else {
        dataProvider->getImageSize(imW, imH);

        if (!imW || !imH) {
            return;
        }
    }

    PolarCoord polCoord1, polCoord2, polCoord0;
    // dataProvider->getImageSize(imW, imH);
    double decayY = (locY_) * double (imH) / 2000.;
    double decayYT = (locYT_) * double (imH) / 2000.;
    double decayX = (locX_) * (double (imW)) / 2000.;
    double decayXL = (locXL_) * (double (imW)) / 2000.;
    rtengine::Coord origin(imW / 2 + centerX_ * imW / 2000.f, imH / 2 + centerY_ * imH / 2000.f);
//   printf("deX=%f dexL=%f deY=%f deyT=%f locX=%i locY=%i\n", decayX, decayXL, decayY, decayYT, locX_, locY_);

    if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
        decayYT = decayY;
        decayXL = decayX;
    }

//    Line *currLine;
//    Circle *currCircle;
    //  Arcellipse *currArcellipse;
//    Beziers *currBeziers;
    double decay;
    /*
    const auto updateLine = [&] (Geometry * geometry, const float radius, const float begin, const float end) {
        const auto line = static_cast<Line*> (geometry);
        line->begin = PolarCoord (radius, -degree_ + begin);
        line->begin += origin;
        line->end = PolarCoord (radius, -degree_ + end);
        line->end += origin;
    };
    */
    const auto updateLineWithDecay = [&](Geometry * geometry, const float radius, const float decal, const float offSetAngle) {
        const auto line = static_cast<Line*>(geometry);  //180
        line->begin = PolarCoord(radius, -degree_ + decal) + PolarCoord(decay, -degree_ + offSetAngle);
        line->begin += origin;//0
        line->end = PolarCoord(radius, -degree_ + (decal - 180)) + PolarCoord(decay, -degree_ + offSetAngle);
        line->end += origin;
    };

    const auto updateCircle = [&](Geometry * geometry) {
        const auto circle = static_cast<Circle*>(geometry);
        circle->center = origin;
        circle->radius = circrad_;
    };

    const auto updateBeziers = [&](Geometry * geometry, const double dX_, const double dI_, const double dY_,  const float begi, const float inte, const float en) {
        const auto beziers = static_cast<Beziers*>(geometry);
        beziers->begin = PolarCoord(dX_, begi);
        beziers->begin += origin;//0
        beziers->inter = PolarCoord(dI_, inte);
        beziers->inter += origin;//0
        beziers->end = PolarCoord(dY_,  en);
        beziers->end += origin;
        //  printf("dX=%f dI=%f dY=%f begx=%i begy=%i intx=%i inty=%i endx=%i endy=%i\n", dX_, dI_, dY_, beziers->begin.x, beziers->begin.y, beziers->inter.x, beziers->inter.y, beziers->end.x, beziers->end.y);
    };

    /*
        const auto updateArcellipse = [&] (Geometry * geometry, const double dX_, const double dY_, const float kbegang, const float kendang) {
            const auto arcellipse = static_cast<Arcellipse*> (geometry);
            arcellipse->center = origin;
            arcellipse->radius = dY_;
            arcellipse->radius2 = dX_;
            arcellipse->translax = (double) imW /2.; //dX_ - dY_;
            arcellipse->translay = (double) imH /2.;
            arcellipse->scalx = dX_ / dY_; // double(locX_) / double (locY_); //arcellipse->radius2 / arcellipse->radius ; // dX_ / dY_;
            arcellipse->scaly = 1.; //dX_ / dY_; //locY_/locX_;
            arcellipse->begang = kbegang * M_PI;
            arcellipse->endang = kendang * M_PI;


        };
    */
    double dimline = 100.;

    if (options.showdelimspot) {
        dimline = 500.;
    }


    decay = decayX;
    updateLineWithDecay(visibleGeometry.at(0), dimline, 90., 0.);
    updateLineWithDecay(mouseOverGeometry.at(0), dimline, 90., 0.);

    decay = decayXL;

    updateLineWithDecay(visibleGeometry.at(1), dimline, 90., 180.);
    updateLineWithDecay(mouseOverGeometry.at(1), dimline, 90., 180.);

    decay = decayYT;
    updateLineWithDecay(visibleGeometry.at(2), dimline, 180., 270.);
    updateLineWithDecay(mouseOverGeometry.at(2), dimline, 180., 270.);

    decay = decayY;

    updateLineWithDecay(visibleGeometry.at(3), dimline, 180, 90.);
    updateLineWithDecay(mouseOverGeometry.at(3), dimline, 180., 90.);


    updateCircle(visibleGeometry.at(4));
    updateCircle(mouseOverGeometry.at(4));

    if (options.showdelimspot) {
        //this decayww evaluate approximation of a point in the ellipse for an angle alpha
        //this decayww evaluate approximation of a point in the ellipse for an angle alpha
        double decay5 = 1.003819 * ((decayX * decayY) / sqrt(0.00765 * SQR(decayX) + SQR(decayY)));    //0.07179 = SQR(sin(15)/cos(15))  1.0038 = 1 / cos(5)
        double decay15 = 1.03527 * ((decayX * decayY) / sqrt(0.07179 * SQR(decayX) + SQR(decayY)));    //0.07179 = SQR(sin(15)/cos(15))  1.03527 = 1 / cos(15)
        double decay30 = 1.15473 * ((decayX * decayY) / sqrt(0.33335 * SQR(decayX) + SQR(decayY)));
        double decay60 = 2. * ((decayX * decayY) / sqrt(3.0 * SQR(decayX) + SQR(decayY)));
        double decay75 = 3.86398 * ((decayX * decayY) / sqrt(13.929 * SQR(decayX) + SQR(decayY)));
        double decay85 = 11.473 * ((decayX * decayY) / sqrt(130.64 * SQR(decayX) + SQR(decayY)));

        double decay5L = 1.003819 * ((decayXL * decayY) / sqrt(0.00765 * SQR(decayXL) + SQR(decayY)));    //0.07179 = SQR(sin(15)/cos(15))  1.0038 = 1 / cos(5)
        double decay15L = 1.03527 * ((decayXL * decayY) / sqrt(0.07179 * SQR(decayXL) + SQR(decayY)));
        double decay30L = 1.15473 * ((decayXL * decayY) / sqrt(0.33335 * SQR(decayXL) + SQR(decayY)));
        double decay60L = 2. * ((decayXL * decayY) / sqrt(3.0 * SQR(decayXL) + SQR(decayY)));
        double decay75L = 3.86398 * ((decayXL * decayY) / sqrt(13.929 * SQR(decayXL) + SQR(decayY)));
        double decay85L = 11.473 * ((decayXL * decayY) / sqrt(130.64 * SQR(decayXL) + SQR(decayY)));

        double decay5LT = 1.003819 * ((decayXL * decayYT) / sqrt(0.00765 * SQR(decayXL) + SQR(decayYT)));    //0.07179 = SQR(sin(15)/cos(15))  1.0038 = 1 / cos(5)
        double decay15LT = 1.03527 * ((decayXL * decayYT) / sqrt(0.07179 * SQR(decayXL) + SQR(decayYT)));
        double decay30LT = 1.15473 * ((decayXL * decayYT) / sqrt(0.33335 * SQR(decayXL) + SQR(decayYT)));
        double decay60LT = 2. * ((decayXL * decayYT) / sqrt(3.0 * SQR(decayXL) + SQR(decayYT)));
        double decay75LT = 3.86398 * ((decayXL * decayYT) / sqrt(13.929 * SQR(decayXL) + SQR(decayYT)));
        double decay85LT = 11.473 * ((decayXL * decayYT) / sqrt(130.64 * SQR(decayXL) + SQR(decayYT)));

        double decay5T = 1.003819 * ((decayX * decayYT) / sqrt(0.00765 * SQR(decayX) + SQR(decayYT)));    //0.07179 = SQR(sin(15)/cos(15))  1.0038 = 1 / cos(5)
        double decay15T = 1.03527 * ((decayX * decayYT) / sqrt(0.07179 * SQR(decayX) + SQR(decayYT)));
        double decay30T = 1.15473 * ((decayX * decayYT) / sqrt(0.33335 * SQR(decayX) + SQR(decayYT)));
        double decay60T = 2. * ((decayX * decayYT) / sqrt(3.0 * SQR(decayX) + SQR(decayYT)));
        double decay75T = 3.86398 * ((decayX * decayYT) / sqrt(13.929 * SQR(decayX) + SQR(decayYT)));
        double decay85T = 11.473 * ((decayX * decayYT) / sqrt(130.64 * SQR(decayX) + SQR(decayYT)));

        double decay45 = (1.414 * decayX * decayY) / sqrt(SQR(decayX) + SQR(decayY));
        double decay45L = (1.414 * decayXL * decayY) / sqrt(SQR(decayXL) + SQR(decayY));
        double decay45LT = (1.414 * decayXL * decayYT) / sqrt(SQR(decayXL) + SQR(decayYT));
        double decay45T = (1.414 * decayX * decayYT) / sqrt(SQR(decayX) + SQR(decayYT));

        //printf("decayX=%f decayY=%f decay10=%f decay45=%f oriX=%i origY=%i\n", decayX, decayY, decay10, decay45, origin.x, origin.y);
        updateBeziers(visibleGeometry.at(5), decayX, decay5, decay15, 0., 5., 15.);
        updateBeziers(mouseOverGeometry.at(5), decayX, decay5, decay15, 0., 5., 15.);

        updateBeziers(visibleGeometry.at(6), decay15, decay30, decay45, 15., 30., 45.);
        updateBeziers(mouseOverGeometry.at(6), decay15, decay30, decay45, 15., 30., 45.);

        updateBeziers(visibleGeometry.at(7), decay45, decay60, decay75, 45., 60., 75.);
        updateBeziers(mouseOverGeometry.at(7), decay45, decay60, decay75, 45., 60., 75.);

        updateBeziers(visibleGeometry.at(8), decay75, decay85, decayY, 75., 85., 90.);
        updateBeziers(mouseOverGeometry.at(8), decay75, decay85, decayY, 75., 85., 90.);

        updateBeziers(visibleGeometry.at(9), decayY, decay85L, decay75L, 90., 95., 105.);
        updateBeziers(mouseOverGeometry.at(9), decayY, decay85L, decay75L, 90., 95., 105.);

        updateBeziers(visibleGeometry.at(10), decay75L, decay60L, decay45L, 105., 120., 135.);
        updateBeziers(mouseOverGeometry.at(10), decay75L, decay60L, decay45L, 105., 120., 135.);

        updateBeziers(visibleGeometry.at(11), decay45L, decay30L, decay15L, 135., 150., 165.);
        updateBeziers(mouseOverGeometry.at(11), decay45L, decay30L, decay15L, 135., 150., 165.);

        updateBeziers(visibleGeometry.at(12), decay15L, decay5L, decayXL, 165., 175., 180.);
        updateBeziers(mouseOverGeometry.at(12), decay15L, decay5L, decayXL, 165., 175., 180.);


        updateBeziers(visibleGeometry.at(13), decayXL, decay5LT, decay15LT, 180., 185., 195.);
        updateBeziers(mouseOverGeometry.at(13), decayXL, decay5LT, decay15LT, 180., 185., 195.);

        updateBeziers(visibleGeometry.at(14), decay15LT, decay30LT, decay45LT, 195., 210., 225.);
        updateBeziers(mouseOverGeometry.at(14), decay15LT, decay30LT, decay45LT, 195., 210., 225.);

        updateBeziers(visibleGeometry.at(15), decay45LT, decay60LT, decay75LT, 225., 240., 255.);
        updateBeziers(mouseOverGeometry.at(15), decay45LT, decay60LT, decay75LT, 225., 240., 255.);

        updateBeziers(visibleGeometry.at(16), decay75LT, decay85LT, decayYT, 255., 265., 270.);
        updateBeziers(mouseOverGeometry.at(16), decay75LT, decay85LT, decayYT, 255., 265., 270.);

        updateBeziers(visibleGeometry.at(17), decayYT, decay85T, decay75T, 270., 275., 285.);
        updateBeziers(mouseOverGeometry.at(17), decayYT, decay85T, decay75T, 270., 275., 285.);

        updateBeziers(visibleGeometry.at(18), decay75T, decay60T, decay45T, 285., 300., 315.);
        updateBeziers(mouseOverGeometry.at(18), decay75T, decay60T, decay45T, 285., 300., 315.);

        updateBeziers(visibleGeometry.at(19), decay45T, decay30T, decay15T, 315., 330., 345.);
        updateBeziers(mouseOverGeometry.at(19), decay45T, decay30T, decay15T, 315., 330., 345.);

        updateBeziers(visibleGeometry.at(20), decay15T, decay5T, decayX, 345., 355., 360.);
        updateBeziers(mouseOverGeometry.at(20), decay15T, decay5T, decayX, 345., 355., 360.);

    }

    //  updateArcellipse (visibleGeometry.at (5), decayX, decayY, 0., 0.5);
    //  updateArcellipse (mouseOverGeometry.at (5), decayX, decayY, 0., 0.5);

}

void Locallab::write(ProcParams* pp, ParamsEdited* pedited)
{
    pp->locallab.degree = degree->getValue();
    pp->locallab.locY = locY->getIntValue();
    pp->locallab.locX = locX->getValue();
    pp->locallab.locYT = locYT->getIntValue();
    pp->locallab.locXL = locXL->getValue();
    pp->locallab.centerX = centerX->getIntValue();
    pp->locallab.centerY = centerY->getIntValue();
    pp->locallab.circrad = circrad->getIntValue();
    pp->locallab.centerXbuf = centerXbuf->getIntValue();
    pp->locallab.centerYbuf = centerYbuf->getIntValue();
    pp->locallab.adjblur = adjblur->getIntValue();
    pp->locallab.proxi = proxi->getIntValue();
    pp->locallab.thres = thres->getIntValue();
    pp->locallab.lightness = lightness->getIntValue();
    pp->locallab.contrast = contrast->getIntValue();
    pp->locallab.chroma = chroma->getIntValue();
    pp->locallab.expcomp = expcomp->getValue();
    pp->locallab.black = (int)black->getValue();
    pp->locallab.hlcompr = (int)hlcompr->getValue();
    pp->locallab.hlcomprthresh = (int)hlcomprthresh->getValue();
    pp->locallab.shcompr = (int)shcompr->getValue();
    pp->locallab.noiselumc = noiselumc->getIntValue();
    pp->locallab.noiselumf = noiselumf->getIntValue();
    pp->locallab.noisechrof = noisechrof->getIntValue();
    pp->locallab.noisechroc = noisechroc->getIntValue();
    pp->locallab.sharradius = sharradius->getIntValue();
    pp->locallab.sharamount = sharamount->getIntValue();
    pp->locallab.shardamping = shardamping->getIntValue();
    pp->locallab.shariter = shariter->getIntValue();
    pp->locallab.sensisha = sensisha->getIntValue();
    pp->locallab.sensi = sensi->getIntValue();
    pp->locallab.sensiex = sensiex->getIntValue();
    pp->locallab.sensih = sensih->getIntValue();
    pp->locallab.retrab = retrab->getIntValue();
    pp->locallab.sensicb = sensicb->getIntValue();
    pp->locallab.sensiexclu = sensiexclu->getIntValue();
    pp->locallab.struc = struc->getIntValue();
    pp->locallab.sensibn = sensibn->getIntValue();
    pp->locallab.sensitm = sensitm->getIntValue();
    pp->locallab.radius = radius->getIntValue();
    pp->locallab.strength = strength->getIntValue();
    pp->locallab.stren = stren->getIntValue();
    pp->locallab.gamma = gamma->getIntValue();
    pp->locallab.estop = estop->getIntValue();
    pp->locallab.scaltm = scaltm->getIntValue();
    pp->locallab.rewei = rewei->getIntValue();
    pp->locallab.enabled = getEnabled();
    pp->locallab.transit = transit->getIntValue();
    pp->locallab.avoid = avoid->get_active();
    pp->locallab.activlum = activlum->get_active();
    pp->locallab.invers = invers->get_active();
    pp->locallab.cutpast = cutpast->get_active();
    pp->locallab.lastdust = lastdust->get_active();
    pp->locallab.curvactiv = curvactiv->get_active();
    pp->locallab.inversrad = inversrad->get_active();
    pp->locallab.inversret = inversret->get_active();
    pp->locallab.inverssha = inverssha->get_active();
    pp->locallab.str = str->getIntValue();
    pp->locallab.neigh = neigh->getIntValue();
    pp->locallab.nbspot = nbspot->getIntValue();
    pp->locallab.anbspot = anbspot->getIntValue();
    pp->locallab.hueref = hueref->getValue();
    pp->locallab.chromaref = chromaref->getValue();
    pp->locallab.lumaref = lumaref->getValue();
    pp->locallab.sobelref = sobelref->getValue();
    pp->locallab.vart = vart->getIntValue();
    pp->locallab.chrrt = chrrt->getIntValue();
    pp->locallab.localTgaincurve       = cTgainshape->getCurve();
    pp->locallab.localTgaincurverab       = cTgainshaperab->getCurve();
    pp->locallab.llcurve       = llshape->getCurve();
    pp->locallab.cccurve       = ccshape->getCurve();
    pp->locallab.LHcurve       = LHshape->getCurve();
    pp->locallab.HHcurve       = HHshape->getCurve();
    pp->locallab.excurve       = shape->getCurve();
    pp->locallab.expcolor      = expcolor->getEnabled();
    pp->locallab.expexpose      = expexpose->getEnabled();
    pp->locallab.expvibrance      = expvibrance->getEnabled();
    pp->locallab.expblur      = expblur->getEnabled();
    pp->locallab.exptonemap      = exptonemap->getEnabled();
    pp->locallab.expreti      = expreti->getEnabled();
    pp->locallab.expsharp      = expsharp->getEnabled();
    pp->locallab.expcbdl      = expcbdl->getEnabled();
    pp->locallab.expdenoi      = expdenoi->getEnabled();

    for (int i = 0; i < 5; i++) {
        pp->locallab.mult[i] = multiplier[i]->getIntValue();
    }

    pp->locallab.threshold = threshold->getIntValue();
    pp->locallab.chromacbdl = chromacbdl->getIntValue();

    pp->locallab.pastels         = pastels->getIntValue();
    pp->locallab.saturated       = pastSatTog->get_active() ? pp->locallab.pastels : saturated->getIntValue();
    pp->locallab.psthreshold     = psThreshold->getValue<int> ();
    pp->locallab.protectskins    = protectSkins->get_active();
    pp->locallab.avoidcolorshift = avoidColorShift->get_active();
    pp->locallab.pastsattog      = pastSatTog->get_active();
    pp->locallab.skintonescurve  = skinTonesCurve->getCurve();
    pp->locallab.sensiv = sensiv->getIntValue();


    if (pedited) {
        pedited->locallab.degree = degree->getEditedState();
        pedited->locallab.Smethod  = Smethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.Exclumethod  = Exclumethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.retinexMethod    = retinexMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.qualityMethod    = qualityMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.qualitycurveMethod    = qualitycurveMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.blurMethod    = blurMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.dustMethod    = dustMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.locY = locY->getEditedState();
        pedited->locallab.locX = locX->getEditedState();
        pedited->locallab.locYT = locYT->getEditedState();
        pedited->locallab.locXL = locXL->getEditedState();
        pedited->locallab.centerX = centerX->getEditedState();
        pedited->locallab.centerY = centerY->getEditedState();
        pedited->locallab.circrad = circrad->getEditedState();
        pedited->locallab.centerXbuf = centerXbuf->getEditedState();
        pedited->locallab.centerYbuf = centerYbuf->getEditedState();
        pedited->locallab.adjblur = adjblur->getEditedState();
        pedited->locallab.proxi = proxi->getEditedState();
        pedited->locallab.thres = thres->getEditedState();
        pedited->locallab.lightness = lightness->getEditedState();
        pedited->locallab.contrast = contrast->getEditedState();
        pedited->locallab.chroma = chroma->getEditedState();
        pedited->locallab.expcomp    = expcomp->getEditedState();
        pedited->locallab.black      = black->getEditedState();
        pedited->locallab.hlcompr    = hlcompr->getEditedState();
        pedited->locallab.hlcomprthresh = hlcomprthresh->getEditedState();
        pedited->locallab.shcompr    = shcompr->getEditedState();

        pedited->locallab.noiselumf = noiselumf->getEditedState();
        pedited->locallab.noiselumc = noiselumc->getEditedState();
        pedited->locallab.noisechrof = noisechrof->getEditedState();
        pedited->locallab.noisechroc = noisechroc->getEditedState();
        pedited->locallab.sharradius = sharradius->getEditedState();
        pedited->locallab.sharamount = sharamount->getEditedState();
        pedited->locallab.shardamping = shardamping->getEditedState();
        pedited->locallab.shariter = shariter->getEditedState();
        pedited->locallab.sensisha = sensisha->getEditedState();
        pedited->locallab.sensi = sensi->getEditedState();
        pedited->locallab.sensiex = sensiex->getEditedState();
        pedited->locallab.sensih = sensih->getEditedState();
        pedited->locallab.retrab = retrab->getEditedState();
        pedited->locallab.sensicb = sensicb->getEditedState();
        pedited->locallab.sensiexclu = sensiexclu->getEditedState();
        pedited->locallab.struc = struc->getEditedState();
        pedited->locallab.sensibn = sensibn->getEditedState();
        pedited->locallab.sensitm = sensitm->getEditedState();
        pedited->locallab.radius = radius->getEditedState();
        pedited->locallab.strength = strength->getEditedState();
        pedited->locallab.stren = stren->getEditedState();
        pedited->locallab.gamma = gamma->getEditedState();
        pedited->locallab.estop = estop->getEditedState();
        pedited->locallab.scaltm = scaltm->getEditedState();
        pedited->locallab.rewei = rewei->getEditedState();
        pedited->locallab.transit = transit->getEditedState();
        pedited->locallab.enabled = !get_inconsistent();
        pedited->locallab.avoid = !avoid->get_inconsistent();
        pedited->locallab.invers = !invers->get_inconsistent();
        pedited->locallab.cutpast = !cutpast->get_inconsistent();
        pedited->locallab.lastdust = !lastdust->get_inconsistent();
        pedited->locallab.curvactiv = !curvactiv->get_inconsistent();
        pedited->locallab.activlum = !activlum->get_inconsistent();
        pedited->locallab.inversret = !inversret->get_inconsistent();
        pedited->locallab.inversrad = !inversrad->get_inconsistent();
        pedited->locallab.inverssha = !inverssha->get_inconsistent();
        pedited->locallab.str = str->getEditedState();
        pedited->locallab.neigh = neigh->getEditedState();
        pedited->locallab.nbspot = nbspot->getEditedState();
        pedited->locallab.anbspot = anbspot->getEditedState();
        pedited->locallab.hueref = hueref->getEditedState();
        pedited->locallab.chromaref = chromaref->getEditedState();
        pedited->locallab.lumaref = lumaref->getEditedState();
        pedited->locallab.sobelref = sobelref->getEditedState();

        pedited->locallab.vart = vart->getEditedState();
        pedited->locallab.chrrt = chrrt->getEditedState();
        pedited->locallab.localTgaincurve        = !cTgainshape->isUnChanged();
        pedited->locallab.localTgaincurverab        = !cTgainshaperab->isUnChanged();
        pedited->locallab.llcurve        = !llshape->isUnChanged();
        pedited->locallab.cccurve        = !ccshape->isUnChanged();
        pedited->locallab.LHcurve        = !LHshape->isUnChanged();
        pedited->locallab.excurve        = !shape->isUnChanged();
        pedited->locallab.HHcurve        = !HHshape->isUnChanged();
        pedited->locallab.expcolor     = !expcolor->get_inconsistent();
        pedited->locallab.expexpose     = !expexpose->get_inconsistent();
        pedited->locallab.expvibrance     = !expvibrance->get_inconsistent();
        pedited->locallab.expblur     = !expblur->get_inconsistent();
        pedited->locallab.exptonemap     = !exptonemap->get_inconsistent();
        pedited->locallab.expreti     = !expreti->get_inconsistent();
        pedited->locallab.expsharp     = !expsharp->get_inconsistent();
        pedited->locallab.expcbdl     = !expcbdl->get_inconsistent();
        pedited->locallab.expdenoi     = !expdenoi->get_inconsistent();

        for (int i = 0; i < 5; i++) {
            pedited->locallab.mult[i] = multiplier[i]->getEditedState();
        }

        pedited->locallab.threshold = threshold->getEditedState();
        pedited->locallab.chromacbdl = chromacbdl->getEditedState();
        pedited->locallab.pastels         = pastels->getEditedState();
        pedited->locallab.saturated       = saturated->getEditedState();
        pedited->locallab.psthreshold     = psThreshold->getEditedState();
        pedited->locallab.protectskins    = !protectSkins->get_inconsistent();
        pedited->locallab.avoidcolorshift = !avoidColorShift->get_inconsistent();
        pedited->locallab.pastsattog      = !pastSatTog->get_inconsistent();
        pedited->locallab.skintonescurve  = !skinTonesCurve->isUnChanged();
        pedited->locallab.sensiv = sensiv->getEditedState();

    }

    if (retinexMethod->get_active_row_number() == 0) {
        pp->locallab.retinexMethod = "low";
    } else if (retinexMethod->get_active_row_number() == 1) {
        pp->locallab.retinexMethod = "uni";
    } else if (retinexMethod->get_active_row_number() == 2) {
        pp->locallab.retinexMethod = "high";
    }

    if (blurMethod->get_active_row_number() == 0) {
        pp->locallab.blurMethod = "norm";
    } else if (blurMethod->get_active_row_number() == 1) {
        pp->locallab.blurMethod = "inv";
    } else if (blurMethod->get_active_row_number() == 2) {
        pp->locallab.blurMethod = "sym";
    }

    if (dustMethod->get_active_row_number() == 0) {
        pp->locallab.dustMethod = "cop";
    } else if (dustMethod->get_active_row_number() == 1) {
        pp->locallab.dustMethod = "mov";
    } else if (dustMethod->get_active_row_number() == 2) {
        pp->locallab.dustMethod = "pas";
    }

    if (qualityMethod->get_active_row_number() == 0) {
        pp->locallab.qualityMethod = "std";
    } else if (qualityMethod->get_active_row_number() == 1) {
        pp->locallab.qualityMethod = "enh";
    } else if (qualityMethod->get_active_row_number() == 2) {
        pp->locallab.qualityMethod = "enhden";
    }

    if (qualitycurveMethod->get_active_row_number() == 0) {
        pp->locallab.qualitycurveMethod = "none";
    } else if (qualitycurveMethod->get_active_row_number() == 1) {
        pp->locallab.qualitycurveMethod = "std";
    } else if (qualitycurveMethod->get_active_row_number() == 2) {
        pp->locallab.qualitycurveMethod = "enh";
    }

    if (Exclumethod->get_active_row_number() == 0) {
        pp->locallab.Exclumethod = "norm";
    } else if (Exclumethod->get_active_row_number() == 1) {
        pp->locallab.Exclumethod = "exc";
    }

    if (Smethod->get_active_row_number() == 0) {
        pp->locallab.Smethod = "IND";
    } else if (Smethod->get_active_row_number() == 1) {
        pp->locallab.Smethod = "SYM";
    } else if (Smethod->get_active_row_number() == 2) {
        pp->locallab.Smethod = "INDSL";
    } else if (Smethod->get_active_row_number() == 3) {
        pp->locallab.Smethod = "SYMSL";
    }

    if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
//   if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {
        pp->locallab.locX = locX->getValue();
        pp->locallab.locY = locY->getValue();

        pp->locallab.locXL = pp->locallab.locX;
        pp->locallab.locYT = pp->locallab.locY;
    }
    /*  else if(Smethod->get_active_row_number()==2){
            pp->locallab.locXL=pp->locallab.locX;
            pp->locallab.locYT=pp->locallab.locX;
            pp->locallab.locY=pp->locallab.locX;
        }
        */
    else {
        pp->locallab.locXL = locXL->getValue();
        pp->locallab.locX = locX->getValue();
        pp->locallab.locY = locY->getValue();
        pp->locallab.locYT = locYT->getValue();
    }
}

void Locallab::protectskins_toggled()
{
    if (batchMode) {
        if (protectSkins->get_inconsistent()) {
            protectSkins->set_inconsistent(false);
            pskinsconn.block(true);
            protectSkins->set_active(false);
            pskinsconn.block(false);
        } else if (lastProtectSkins) {
            protectSkins->set_inconsistent(true);
        }

        lastProtectSkins = protectSkins->get_active();
    }

    if (listener && getEnabled()) {
        if (protectSkins->get_active()) {
            listener->panelChanged(EvlocallabProtectSkins, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvlocallabProtectSkins, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::avoidcolorshift_toggled()
{
    if (batchMode) {
        if (avoidColorShift->get_inconsistent()) {
            avoidColorShift->set_inconsistent(false);
            ashiftconn.block(true);
            avoidColorShift->set_active(false);
            ashiftconn.block(false);
        } else if (lastAvoidColorShift) {
            avoidColorShift->set_inconsistent(true);
        }

        lastAvoidColorShift = avoidColorShift->get_active();
    }

    if (listener && getEnabled()) {
        if (avoidColorShift->get_active()) {
            listener->panelChanged(EvlocallabAvoidColorShift, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvlocallabAvoidColorShift, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::pastsattog_toggled()
{
    if (batchMode) {
        if (pastSatTog->get_inconsistent()) {
            pastSatTog->set_inconsistent(false);
            pastsattogconn.block(true);
            pastSatTog->set_active(false);
            pastsattogconn.block(false);
        } else if (lastPastSatTog) {
            pastSatTog->set_inconsistent(true);
        }

        lastPastSatTog = pastSatTog->get_active();
    }

    if (pastSatTog->get_active()) {
        // Link both slider, so we set saturated and psThresholds unsensitive
        psThreshold->set_sensitive(false);
        saturated->set_sensitive(false);
        saturated->setValue(pastels->getValue());      // Pastels and Saturated are linked
    } else {
        // Separate sliders, so we set saturated and psThresholds sensitive again
        psThreshold->set_sensitive(true);
        saturated->set_sensitive(true);
    }

    if (listener && getEnabled()) {
        if (pastSatTog->get_active()) {
            listener->panelChanged(EvlocallabPastSatTog, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvlocallabPastSatTog, M("GENERAL_DISABLED"));
        }
    }
}




void Locallab::curveChanged(CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == cTgainshape) {
            listener->panelChanged(EvlocallabCTgainCurve, M("HISTORY_CUSTOMCURVE"));  //HISTORY_CUSTOMCURVE
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue(strval + 1);
            adjusterChanged(retrab, strval + 1);
            usleep(10000);  //to test
            retrab->setValue(strval);

            adjusterChanged(retrab, strval);
        }

        else if (ce == cTgainshaperab) {
            listener->panelChanged(EvlocallabCTgainCurverab, M(""));
        } else if (ce == LHshape) {
            listener->panelChanged(EvlocallabLHshape, M(""));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue(strval + 1);
            adjusterChanged(retrab, strval + 1);
            usleep(10000);  //to test
            retrab->setValue(strval);

            adjusterChanged(retrab, strval);
        } else if (ce == HHshape) {
            listener->panelChanged(EvlocallabHHshape, M("HISTORY_CUSTOMCURVE"));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue(strval + 1);
            adjusterChanged(retrab, strval + 1);
            usleep(10000);  //to test
            retrab->setValue(strval);
        } else if (ce == shape) {
            listener->panelChanged(Evlocallabshape, M("HISTORY_CUSTOMCURVE"));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue(strval + 1);
            adjusterChanged(retrab, strval + 1);
            usleep(10000);  //to test
            retrab->setValue(strval);

            adjusterChanged(retrab, strval);


        } else if (ce == llshape) {
            listener->panelChanged(Evlocallabllshape, M("HISTORY_CUSTOMCURVE"));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue(strval + 1);
            adjusterChanged(retrab, strval + 1);
            usleep(10000);  //to test
            retrab->setValue(strval);

            adjusterChanged(retrab, strval);

        } else if (ce == ccshape) {
            listener->panelChanged(Evlocallabccshape, M("HISTORY_CUSTOMCURVE"));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue(strval + 1);
            adjusterChanged(retrab, strval + 1);
            usleep(10000);  //to test
            retrab->setValue(strval);

            adjusterChanged(retrab, strval);


        } else if (ce == skinTonesCurve) {
            listener->panelChanged(EvlocallabSkinTonesCurve, M("HISTORY_CUSTOMCURVE"));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue(strval + 1);
            adjusterChanged(retrab, strval + 1);
            usleep(10000);  //to test
            retrab->setValue(strval);

            adjusterChanged(retrab, strval);



        }


    }
}

void Locallab::retinexMethodChanged()
{
    if (!batchMode) {

        retrab->hide();
        LocalcurveEditorgainTrab->hide();
    }

    if (listener) {
        listener->panelChanged(EvlocallabretinexMethod, retinexMethod->get_active_text());
    }
}

void Locallab::blurMethodChanged()
{
    if (!batchMode) {

    }

    if (blurMethod->get_active_row_number() == 0 || blurMethod->get_active_row_number() == 2) {
        sensibn->show();
    } else {
        sensibn->hide();
    }

    if (listener) {
        listener->panelChanged(EvlocallabblurMethod, blurMethod->get_active_text());
    }
}


void Locallab::dustMethodChanged()
{
    if (!batchMode) {

    }


    if (listener) {
        //   listener->panelChanged (EvlocallabdustMethod, dustMethod->get_active_text ());
    }
}

void Locallab::qualityMethodChanged()
{
    if (!batchMode) {
        /*
        if (qualityMethod->get_active_row_number() == 0) { //STD
            proxi->hide();
            thres->hide();
        } else {//enh
            proxi->show();
            thres->show();
        }
        */
    }

    if (listener) {
        listener->panelChanged(EvlocallabqualityMethod, qualityMethod->get_active_text());
    }
}

void Locallab::qualitycurveMethodChanged()
{
    if (!batchMode) {
        /*
        if (qualitycurveMethod->get_active_row_number() == 0  || qualitycurveMethod->get_active_row_number() == 1) { //None or STD
            artifFrame->hide();
        } else if (qualitycurveMethod->get_active_row_number() == 0 && qualityMethod->get_active_row_number() >= 1) {
            artifFrame->show();
        } else if (qualitycurveMethod->get_active_row_number() == 2){
            artifFrame->show();
        }
        */
    }

    if (listener) {
        listener->panelChanged(EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text());
    }
}

void Locallab::ExclumethodChanged()
{
    if (!batchMode) {
        if (Exclumethod->get_active_row_number() == 0) {
            excluFrame->hide();
        } else if (Exclumethod->get_active_row_number() == 1) {
            excluFrame->show();
        }
    }


    if (listener) {
        listener->panelChanged(Evlocallabexclumethod, Exclumethod->get_active_text());
    }

}


void Locallab::SmethodChanged()
{
    if (!batchMode) {
        if (Smethod->get_active_row_number() == 0) { //IND 0
            locX->hide();
            locXL->hide();
            locY->hide();
            locYT->hide();
            centerX->hide();
            centerY->hide();
        } else if (Smethod->get_active_row_number() == 1) {         // 1 SYM
            locX->hide();
            locXL->hide();
            locY->hide();
            locYT->hide();
            centerX->hide();
            centerY->hide();

        } else if (Smethod->get_active_row_number() == 2) {         //2 SYM
            locX->show();
            locXL->show();
            locY->show();
            locYT->show();
            centerX->show();
            centerY->show();

        } else if (Smethod->get_active_row_number() == 3) {         // 3 SYM
            locX->show();
            locXL->hide();
            locY->show();
            locYT->hide();
            centerX->show();
            centerY->show();

        }

        /*      else if(Smethod->get_active_row_number()==2) {              // LOC
                    locX->show();
                    locXL->hide();
                    locY->hide();
                    locYT->hide();
                }   */
    }

    if (listener && getEnabled()) {
        if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
            listener->panelChanged(EvlocallabSmet, Smethod->get_active_text());
            locXL->setValue(locX->getValue());
            locYT->setValue(locY->getValue());
        }
        //   else if(Smethod->get_active_row_number()==2) {
        //          listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
        //           locXL->setValue (locX->getValue());
        //           locYT->setValue (locX->getValue());
        //          locY->setValue (locX->getValue());
        //     }
        else

        {
            listener->panelChanged(EvlocallabSmet, Smethod->get_active_text());

        }
    }
}
void Locallab::inversChanged()
{

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

    if (invers->get_active()) {
        //sensi->hide();
        sensi->show();
        llCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        artifFrame->show();
        labqualcurv->hide();

    } else {
        sensi->show();
        llCurveEditorG->show();
        curvactiv->show();
        qualitycurveMethod->show();
        artifFrame->show();
        labqualcurv->show();

    }

    if (listener) {
        if (getEnabled()) {
            listener->panelChanged(Evlocallabinvers, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabinvers, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::cutpastChanged()
{

    if (batchMode) {
        if (cutpast->get_inconsistent()) {
            cutpast->set_inconsistent(false);
            cutpastConn.block(true);
            cutpast->set_active(false);
            cutpastConn.block(false);
        } else if (lastcutpast) {
            cutpast->set_inconsistent(true);
        }

        lastcutpast = cutpast->get_active();
    }


    if (listener) {
        if (getEnabled()) {
            //     listener->panelChanged (Evlocallabcutpast, M ("GENERAL_ENABLED"));
        } else {
            //     listener->panelChanged (Evlocallabcutpast, M ("GENERAL_DISABLED"));
        }
    }
}

void Locallab::lastdustChanged()
{

    if (batchMode) {
        if (lastdust->get_inconsistent()) {
            lastdust->set_inconsistent(false);
            lastdustConn.block(true);
            lastdust->set_active(false);
            lastdustConn.block(false);
        } else if (lastlastdust) {
            lastdust->set_inconsistent(true);
        }

        lastlastdust = lastdust->get_active();
    }


    if (listener) {
        if (getEnabled()) {
            //    listener->panelChanged (Evlocallablastdust, M ("GENERAL_ENABLED"));
        } else {
            //    listener->panelChanged (Evlocallablastdust, M ("GENERAL_DISABLED"));
        }
    }
}


void Locallab::curvactivChanged()
{

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


    if (listener) {
        if (getEnabled()) {
            listener->panelChanged(Evlocallabcurvactiv, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabcurvactiv, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::activlumChanged()
{

    if (batchMode) {
        if (activlum->get_inconsistent()) {
            activlum->set_inconsistent(false);
            activlumConn.block(true);
            activlum->set_active(false);
            activlumConn.block(false);
        } else if (lastactivlum) {
            activlum->set_inconsistent(true);
        }

        lastactivlum = activlum->get_active();
    }


    if (listener) {

        if (getEnabled()) {
            listener->panelChanged(Evlocallabactivlum, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabactivlum, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::inversradChanged()
{

    if (batchMode) {
        if (inversrad->get_inconsistent()) {
            inversrad->set_inconsistent(false);
            inversradConn.block(true);
            inversrad->set_active(false);
            inversradConn.block(false);
        } else if (lastinversrad) {
            inversrad->set_inconsistent(true);
        }

        lastinversrad = inversrad->get_active();
    }

    if (inversrad->get_active()) {
        sensibn->hide();
    } else {
        sensibn->show();
    }


    if (listener) {
        if (getEnabled()) {
            listener->panelChanged(Evlocallabinversrad, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabinversrad, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::inversshaChanged()
{

    if (batchMode) {
        if (inverssha->get_inconsistent()) {
            inverssha->set_inconsistent(false);
            inversshaConn.block(true);
            inverssha->set_active(false);
            inversshaConn.block(false);
        } else if (lastinverssha) {
            inverssha->set_inconsistent(true);
        }

        lastinverssha = inverssha->get_active();
    }

    /*
        if(inverssha->get_active ()) {
            sensisha->hide();
        } else {
            sensisha->show();
        }
    */
    if (listener) {
        if (getEnabled()) {
            listener->panelChanged(Evlocallabinverssha, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabinverssha, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::inversretChanged()
{

    if (batchMode) {
        if (inversret->get_inconsistent()) {
            inversret->set_inconsistent(false);
            inversretConn.block(true);
            inversret->set_active(false);
            inversretConn.block(false);
        } else if (lastinversret) {
            inversret->set_inconsistent(true);
        }

        lastinversret = inversret->get_active();
    }

    if (inversret->get_active()) {
        sensih->hide();
    } else {
        sensih->show();
    }


    if (listener) {
        if (getEnabled()) {
            listener->panelChanged(Evlocallabinversret, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabinversret, M("GENERAL_DISABLED"));
        }
    }
}


void Locallab::setDefaults(const ProcParams * defParams, const ParamsEdited * pedited)
{
    degree->setDefault(defParams->locallab.degree);
    locY->setDefault(defParams->locallab.locY);
    locX->setDefault(defParams->locallab.locX);
    locYT->setDefault(defParams->locallab.locYT);
    locXL->setDefault(defParams->locallab.locXL);
    centerX->setDefault(defParams->locallab.centerX);
    centerY->setDefault(defParams->locallab.centerY);
    circrad->setDefault(defParams->locallab.circrad);
    centerXbuf->setDefault(defParams->locallab.centerXbuf);
    centerYbuf->setDefault(defParams->locallab.centerYbuf);
    adjblur->setDefault(defParams->locallab.adjblur);
    thres->setDefault(defParams->locallab.thres);
    proxi->setDefault(defParams->locallab.proxi);
    lightness->setDefault(defParams->locallab.lightness);
    contrast->setDefault(defParams->locallab.contrast);
    chroma->setDefault(defParams->locallab.chroma);
    expcomp->setDefault(defParams->locallab.expcomp);
    black->setDefault(defParams->locallab.black);
    hlcompr->setDefault(defParams->locallab.hlcompr);
    hlcomprthresh->setDefault(defParams->locallab.hlcomprthresh);
    shcompr->setDefault(defParams->locallab.shcompr);

    noiselumf->setDefault(defParams->locallab.noiselumf);
    noiselumc->setDefault(defParams->locallab.noiselumc);
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
    nbspot->setDefault(defParams->locallab.nbspot);
    anbspot->setDefault(defParams->locallab.anbspot);
    hueref->setDefault(defParams->locallab.hueref);
    chromaref->setDefault(defParams->locallab.chromaref);
    lumaref->setDefault(defParams->locallab.lumaref);
    sobelref->setDefault(defParams->locallab.sobelref);

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
        centerXbuf->setDefaultEditedState(pedited->locallab.centerXbuf ? Edited : UnEdited);
        centerYbuf->setDefaultEditedState(pedited->locallab.centerYbuf ? Edited : UnEdited);
        adjblur->setDefaultEditedState(pedited->locallab.adjblur ? Edited : UnEdited);
        thres->setDefaultEditedState(pedited->locallab.thres ? Edited : UnEdited);
        proxi->setDefaultEditedState(pedited->locallab.proxi ? Edited : UnEdited);
        lightness->setDefaultEditedState(pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setDefaultEditedState(pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setDefaultEditedState(pedited->locallab.chroma ? Edited : UnEdited);
        expcomp->setDefaultEditedState(pedited->locallab.expcomp ? Edited : UnEdited);
        black->setDefaultEditedState(pedited->locallab.black ? Edited : UnEdited);
        hlcompr->setDefaultEditedState(pedited->locallab.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setDefaultEditedState(pedited->locallab.hlcomprthresh ? Edited : UnEdited);
        shcompr->setDefaultEditedState(pedited->locallab.shcompr ? Edited : UnEdited);

        noiselumf->setDefaultEditedState(pedited->locallab.noiselumf ? Edited : UnEdited);
        noiselumc->setDefaultEditedState(pedited->locallab.noiselumc ? Edited : UnEdited);
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
        nbspot->setDefaultEditedState(pedited->locallab.nbspot ? Edited : UnEdited);
        anbspot->setDefaultEditedState(pedited->locallab.anbspot ? Edited : UnEdited);
        hueref->setDefaultEditedState(pedited->locallab.hueref ? Edited : UnEdited);
        chromaref->setDefaultEditedState(pedited->locallab.chromaref ? Edited : UnEdited);
        lumaref->setDefaultEditedState(pedited->locallab.lumaref ? Edited : UnEdited);
        sobelref->setDefaultEditedState(pedited->locallab.sobelref ? Edited : UnEdited);
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
        centerXbuf->setDefaultEditedState(Irrelevant);
        centerYbuf->setDefaultEditedState(Irrelevant);
        adjblur->setDefaultEditedState(Irrelevant);
        thres->setDefaultEditedState(Irrelevant);
        proxi->setDefaultEditedState(Irrelevant);
        lightness->setDefaultEditedState(Irrelevant);
        contrast->setDefaultEditedState(Irrelevant);
        chroma->setDefaultEditedState(Irrelevant);
        expcomp->setDefaultEditedState(Irrelevant);
        black->setDefaultEditedState(Irrelevant);
        hlcompr->setDefaultEditedState(Irrelevant);
        hlcomprthresh->setDefaultEditedState(Irrelevant);
        shcompr->setDefaultEditedState(Irrelevant);

        noiselumf->setDefaultEditedState(Irrelevant);
        noiselumc->setDefaultEditedState(Irrelevant);
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
        nbspot->setDefaultEditedState(Irrelevant);
        anbspot->setDefaultEditedState(Irrelevant);
        hueref->setDefaultEditedState(Irrelevant);
        chromaref->setDefaultEditedState(Irrelevant);
        lumaref->setDefaultEditedState(Irrelevant);
        sobelref->setDefaultEditedState(Irrelevant);
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
}

void Locallab::adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop)
{
    if (listener && getEnabled()) {
        listener->panelChanged(EvlocallabPastSatThreshold, psThreshold->getHistoryString());
    }
}


void Locallab::adjusterChanged(Adjuster * a, double newval)
{

    updateGeometry(int (centerX->getValue()), int (centerY->getValue()), int (circrad->getValue()), (int)locY->getValue(), degree->getValue(), (int)locX->getValue(), (int)locYT->getValue(), (int)locXL->getValue());
    anbspot->hide();
    hueref->hide();
    chromaref->hide();
    lumaref->hide();
    sobelref->hide();
    centerXbuf->hide();
    centerYbuf->hide();

    if (a == pastels && pastSatTog->get_active()) {
        saturated->setValue(newval);
    }

    if (listener && getEnabled()) {
        if (a == degree) {
            listener->panelChanged(EvlocallabDegree, degree->getTextValue());
        } else if (a == locY) {
            if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) { // 0 2
                listener->panelChanged(EvlocallablocY, locY->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocY, locY->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged(EvlocallablocY, locY->getTextValue());
                locYT->setValue(locY->getValue());
            }
        } else if (a == locX) {
            //listener->panelChanged (EvlocallablocX, locX->getTextValue());
            if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
                listener->panelChanged(EvlocallablocX, locX->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged(EvlocallablocX, locX->getTextValue());
                locXL->setValue(locX->getValue());
            }
        } else if (a == locYT) {
            if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
                listener->panelChanged(EvlocallablocYT, locYT->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocYT, locYT->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            else if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
                listener->panelChanged(EvlocallablocYT, locYT->getTextValue());
                locYT->setValue(locY->getValue());
            }
        } else if (a == locXL) {
            if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
                listener->panelChanged(EvlocallablocXL, locXL->getTextValue());
                listener->panelChanged(EvlocallablocXL, locXL->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocXL, locXL->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged(EvlocallablocXL, locXL->getTextValue());
                locXL->setValue(locX->getValue());
            }
        } else if (a == lightness) {
            listener->panelChanged(Evlocallablightness, lightness->getTextValue());
        } else if (a == contrast) {
            listener->panelChanged(Evlocallabcontrast, contrast->getTextValue());
        } else if (a == chroma) {
            listener->panelChanged(Evlocallabchroma, chroma->getTextValue());
        } else if (a == expcomp) {
            listener->panelChanged(Evlocallabexpcomp, expcomp->getTextValue());
        } else if (a == hlcompr) {
            listener->panelChanged(Evlocallabhlcompr, hlcompr->getTextValue());
        } else if (a == hlcomprthresh) {
            listener->panelChanged(Evlocallabhlcomprthresh, hlcomprthresh->getTextValue());
        } else if (a == black) {
            listener->panelChanged(Evlocallabblack, black->getTextValue());
            shcompr->set_sensitive(!((int)black->getValue() == 0));     //at black=0 shcompr value has no effect
        } else if (a == shcompr) {
            listener->panelChanged(Evlocallabshcompr, shcompr->getTextValue());
        } else if (a == sensiex) {
            listener->panelChanged(Evlocallabsensiex, sensiex->getTextValue());


        } else if (a == pastels) {
            listener->panelChanged(EvlocallabPastels, pastels->getTextValue());
        } else if (a == saturated && !pastSatTog->get_active()) {
            listener->panelChanged(EvlocallabSaturated, saturated->getTextValue());
        } else if (a == sensiv) {
            listener->panelChanged(Evlocallabsensiv, sensiv->getTextValue());





        } else if (a == noiselumf) {
            listener->panelChanged(Evlocallabnoiselumf, noiselumf->getTextValue());
        } else if (a == noiselumc) {
            listener->panelChanged(Evlocallabnoiselumc, noiselumc->getTextValue());
        } else if (a == noisechrof) {
            listener->panelChanged(Evlocallabnoisechrof, noisechrof->getTextValue());
        } else if (a == noisechroc) {
            listener->panelChanged(Evlocallabnoisechroc, noisechroc->getTextValue());
        } else if (a == sharradius) {
            listener->panelChanged(Evlocallabsharradius, sharradius->getTextValue());
        } else if (a == sharamount) {
            listener->panelChanged(Evlocallabsharamount, sharamount->getTextValue());
        } else if (a == shardamping) {
            listener->panelChanged(Evlocallabshardamping, shardamping->getTextValue());
        } else if (a == shariter) {
            listener->panelChanged(Evlocallabshariter, shariter->getTextValue());
        } else if (a == sensisha) {
            listener->panelChanged(Evlocallabsensis, sensisha->getTextValue());
        } else if (a == sensi) {
            listener->panelChanged(Evlocallabsensi, sensi->getTextValue());
        } else if (a == sensih) {
            listener->panelChanged(Evlocallabsensih, sensih->getTextValue());
        } else if (a == retrab) {
            listener->panelChanged(Evlocallabretrab, ""); //retrab->getTextValue());
        } else if (a == radius) {
            listener->panelChanged(Evlocallabradius, radius->getTextValue());
        } else if (a == strength) {
            listener->panelChanged(Evlocallabstrength, strength->getTextValue());
        } else if (a == stren) {
            listener->panelChanged(Evlocallabstren, stren->getTextValue());
        } else if (a == gamma) {
            listener->panelChanged(Evlocallabgamma, gamma->getTextValue());
        } else if (a == estop) {
            listener->panelChanged(Evlocallabestop, estop->getTextValue());
        } else if (a == scaltm) {
            listener->panelChanged(Evlocallabscaltm, scaltm->getTextValue());
        } else if (a == rewei) {
            listener->panelChanged(Evlocallabrewei, rewei->getTextValue());
        } else if (a == sensitm) {
            listener->panelChanged(Evlocallabsensitm, sensitm->getTextValue());
        } else if (a == transit) {
            listener->panelChanged(Evlocallabtransit, transit->getTextValue());
        } else if (a == str) {
            listener->panelChanged(Evlocallabstr, str->getTextValue());
        } else if (a == neigh) {
            listener->panelChanged(Evlocallabneigh, neigh->getTextValue());
        } else if (a == nbspot) {
            listener->panelChanged(Evlocallabnbspot, nbspot->getTextValue());
        } else if (a == anbspot) {
            listener->panelChanged(Evlocallabanbspot, ""); //anbspot->getTextValue());
        } else if (a == hueref) {
            listener->panelChanged(Evlocallabhueref, ""); //anbspot->getTextValue());
        } else if (a == chromaref) {
            listener->panelChanged(Evlocallabchromaref, ""); //anbspot->getTextValue());
        } else if (a == lumaref) {
            listener->panelChanged(Evlocallablumaref, ""); //anbspot->getTextValue());
        } else if (a == sobelref) {
            listener->panelChanged(Evlocallabsobelref, ""); //anbspot->getTextValue());
        } else if (a == vart) {
            listener->panelChanged(Evlocallabvart, vart->getTextValue());
        } else if (a == chrrt) {
            listener->panelChanged(Evlocallabchrrt, chrrt->getTextValue());
        } else if (a == circrad) {
            listener->panelChanged(Evlocallabcircrad, circrad->getTextValue());
        } else if (a == thres) {
            listener->panelChanged(Evlocallabthres, thres->getTextValue());
        } else if (a == threshold) {
            listener->panelChanged(EvlocallabThresho, threshold->getTextValue());
        } else if (a == chromacbdl) {
            listener->panelChanged(Evlocallabchromacbdl, chromacbdl->getTextValue());
        } else if (a == sensicb) {
            listener->panelChanged(Evlocallabsensicb, sensicb->getTextValue());
        } else if (a == sensiexclu) {
            listener->panelChanged(Evlocallabsensiexclu, sensiexclu->getTextValue());
        } else if (a == struc) {
            listener->panelChanged(Evlocallabstruc, struc->getTextValue());
        } else if (a == sensibn) {
            listener->panelChanged(Evlocallabsensibn, sensibn->getTextValue());
        } else if (a == proxi) {
            listener->panelChanged(Evlocallabproxi, proxi->getTextValue());
        } else if (a == adjblur) {
            //    listener->panelChanged (Evlocallabadjblur, adjblur->getTextValue());
        } else if (a == centerX || a == centerY) {
            listener->panelChanged(EvlocallabCenter, Glib::ustring::compose("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
        } else if (a == centerXbuf || a == centerYbuf) {
            listener->panelChanged(EvlocallabCenterbuf, Glib::ustring::compose("X=%1\nY=%2", centerXbuf->getTextValue(), centerYbuf->getTextValue()));
        } else {
            listener->panelChanged(EvlocallabEqualizer,
                                   Glib::ustring::compose("%1, %2, %3, %4, %5",
                                           Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[0]->getValue()),
                                           Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[1]->getValue()),
                                           Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[2]->getValue()),
                                           Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[3]->getValue()),
                                           Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[4]->getValue()))
                                  );

        }
    }
}

void Locallab::enabledChanged()
{
    anbspot->hide();
    hueref->hide();
    chromaref->hide();
    lumaref->hide();
    sobelref->hide();

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

    if (batchMode) {
        if (avoid->get_inconsistent()) {
            avoid->set_inconsistent(false);
            avoidConn.block(true);
            avoid->set_active(false);
            avoidConn.block(false);
        } else if (lastavoid) {
            avoid->set_inconsistent(true);
        }

        lastavoid = avoid->get_active();
    }

    if (listener) {
        if (getEnabled()) {
            listener->panelChanged(Evlocallabavoid, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabavoid, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::setAdjusterBehavior(bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd,  bool strengthadd)
{
    degree->setAddMode(degreeadd);
    locY->setAddMode(locYadd);
    locX->setAddMode(locXadd);
    locYT->setAddMode(locYTadd);
    locXL->setAddMode(locXLadd);
    centerX->setAddMode(centeradd);
    centerY->setAddMode(centeradd);
    lightness->setAddMode(lightnessadd);
    contrast->setAddMode(contrastadd);
    chroma->setAddMode(chromaadd);
    sensi->setAddMode(sensiadd);
    transit->setAddMode(transitadd);
    radius->setAddMode(radiusadd);
    strength->setAddMode(strengthadd);

}

void Locallab::trimValues(rtengine::procparams::ProcParams * pp)
{
    degree->trimValue(pp->locallab.degree);
    locY->trimValue(pp->locallab.locY);
    locX->trimValue(pp->locallab.locX);
    locYT->trimValue(pp->locallab.locYT);
    locXL->trimValue(pp->locallab.locXL);
    centerX->trimValue(pp->locallab.centerX);
    centerY->trimValue(pp->locallab.centerY);
    circrad->trimValue(pp->locallab.circrad);
    centerXbuf->trimValue(pp->locallab.centerXbuf);
    centerYbuf->trimValue(pp->locallab.centerYbuf);
    adjblur->trimValue(pp->locallab.adjblur);
    thres->trimValue(pp->locallab.thres);
    proxi->trimValue(pp->locallab.proxi);
    lightness->trimValue(pp->locallab.lightness);
    contrast->trimValue(pp->locallab.contrast);
    chroma->trimValue(pp->locallab.chroma);
    expcomp->trimValue(pp->locallab.expcomp);
    hlcompr->trimValue(pp->locallab.hlcompr);
    hlcomprthresh->trimValue(pp->locallab.hlcomprthresh);
    black->trimValue(pp->locallab.black);
    shcompr->trimValue(pp->locallab.shcompr);

    noiselumf->trimValue(pp->locallab.noiselumf);
    noiselumc->trimValue(pp->locallab.noiselumc);
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
    retrab->trimValue(pp->locallab.retrab);
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
    nbspot->trimValue(pp->locallab.nbspot);
    anbspot->trimValue(pp->locallab.anbspot);
    hueref->trimValue(pp->locallab.hueref);
    chromaref->trimValue(pp->locallab.chromaref);
    lumaref->trimValue(pp->locallab.lumaref);
    sobelref->trimValue(pp->locallab.sobelref);

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

}

void Locallab::setBatchMode(bool batchMode)
{
    removeIfThere(this, edit, false);
    ToolPanel::setBatchMode(batchMode);

    hueref->hide();
    chromaref->hide();
    lumaref->hide();
    sobelref->hide();

    degree->showEditedCB();
    locY->showEditedCB();
    locX->showEditedCB();
    locYT->showEditedCB();
    locXL->showEditedCB();
    centerX->showEditedCB();
    centerY->showEditedCB();
    circrad->showEditedCB();
    centerXbuf->showEditedCB();
    centerYbuf->showEditedCB();
    adjblur->showEditedCB();
    thres->showEditedCB();
    proxi->showEditedCB();
    lightness->showEditedCB();
    contrast->showEditedCB();
    chroma->showEditedCB();
    expcomp->showEditedCB();
    black->showEditedCB();
    hlcompr->showEditedCB();
    hlcomprthresh->showEditedCB();
    shcompr->showEditedCB();

    noiselumf->showEditedCB();
    noiselumc->showEditedCB();
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
    retrab->showEditedCB();
    sensiexclu->showEditedCB();
    struc->showEditedCB();
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
    transit->showEditedCB();
    Smethod->append(M("GENERAL_UNCHANGED"));
    str->showEditedCB();
    neigh->showEditedCB();
    nbspot->showEditedCB();
    anbspot->showEditedCB();
    hueref->showEditedCB();
    chromaref->showEditedCB();
    lumaref->showEditedCB();
    sobelref->showEditedCB();
    vart->showEditedCB();
    LocalcurveEditorgainT->setBatchMode(batchMode);
    LocalcurveEditorgainTrab->setBatchMode(batchMode);
    llCurveEditorG->setBatchMode(batchMode);
//    llCurveEditorG2->setBatchMode (batchMode);
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
    EditSubscriber::setEditProvider(provider);
    cTgainshape->setEditProvider(provider);
    cTgainshaperab->setEditProvider(provider);

}

void Locallab::editToggled()
{
    if (edit->get_active()) {
        subscribe();
    } else {
        unsubscribe();
    }
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



CursorShape Locallab::getCursor(int objectID)
{
    switch (objectID) {
        case (2): {
            int angle = degree->getIntValue();

            if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
                return CSMove1DV;
            }

            return CSMove1DH;
        }

        case (3): {
            int angle = degree->getIntValue();

            if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
                return CSMove1DV;
            }

            return CSMove1DH;
        }

        case (0): {
            int angle = degree->getIntValue();

            if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
                return CSMove1DH;
            }

            return CSMove1DV;
        }

        case (1): {
            int angle = degree->getIntValue();

            if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
                return CSMove1DH;
            }

            return CSMove1DV;
        }

        case (4):
            return CSMove2D;

        default:
            return CSOpenHand;
    }
}

bool Locallab::mouseOver(int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (editProvider && editProvider->object != lastObject) {
        if (lastObject > -1) {
            if (lastObject == 2 || lastObject == 3) {
                EditSubscriber::visibleGeometry.at(2)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at(3)->state = Geometry::NORMAL;

            } else if (lastObject == 0 || lastObject == 1) {
                EditSubscriber::visibleGeometry.at(0)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at(1)->state = Geometry::NORMAL;

            }

            else {
                EditSubscriber::visibleGeometry.at(4)->state = Geometry::NORMAL;
//               EditSubscriber::visibleGeometry.at (lastObject)->state = Geometry::NORMAL;
            }
        }

        if (editProvider->object > -1) {
            if (editProvider->object == 2 || editProvider->object == 3) {
                EditSubscriber::visibleGeometry.at(2)->state = Geometry::PRELIGHT;
                EditSubscriber::visibleGeometry.at(3)->state = Geometry::PRELIGHT;

            } else if (editProvider->object == 0 || editProvider->object == 1) {
                EditSubscriber::visibleGeometry.at(0)->state = Geometry::PRELIGHT;
                EditSubscriber::visibleGeometry.at(1)->state = Geometry::PRELIGHT;

            }

            else {
                EditSubscriber::visibleGeometry.at(4)->state = Geometry::PRELIGHT;
                //              EditSubscriber::visibleGeometry.at (editProvider->object)->state = Geometry::PRELIGHT;
            }
        }

        lastObject = editProvider->object;
        return true;
    }

    return false;
}

bool Locallab::button1Pressed(int modifierKey)
{
    if (lastObject < 0) {
        return false;
    }

    EditDataProvider *provider = getEditProvider();

    if (!(modifierKey & GDK_CONTROL_MASK)) {
        // button press is valid (no modifier key)
        PolarCoord pCoord;
        //  EditDataProvider *provider = getEditProvider();
        int imW, imH;
        provider->getImageSize(imW, imH);
        double halfSizeW = imW / 2.;
        double halfSizeH = imH / 2.;
        draggedCenter.set(int (halfSizeW + halfSizeW * (centerX->getValue() / 1000.)), int (halfSizeH + halfSizeH * (centerY->getValue() / 1000.)));

        // trick to get the correct angle (clockwise/counter-clockwise)
        rtengine::Coord p1 = draggedCenter;
        rtengine::Coord p2 = provider->posImage;
        int p = p1.y;
        p1.y = p2.y;
        p2.y = p;
        pCoord = p2 - p1;
        draggedPointOldAngle = pCoord.angle;
        draggedPointAdjusterAngle = degree->getValue();

        if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
            if (lastObject == 2) {
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double verti = double (imH);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

                if (lastObject == 2) {
                    //draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locYT->getValue() / 2000. * verti);

                }
            } else if (lastObject == 3) {
                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;

                double verti = double (imH);

                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                // draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

                if (lastObject == 3) {
                    draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locY->getValue() / 2000. * verti);

                }

            }

        } else if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
            if (lastObject == 2 || lastObject == 3) {
                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double verti = double (imH);
                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //    draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

                if (lastObject == 3) {
                    draggedlocYOffset = -draggedlocYOffset;
                }

                draggedlocYOffset -= (locY->getValue() / 2000. * verti);
            }
        }

        if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
            if (lastObject == 0) {
                // Dragging a line to change the angle

                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;

                double horiz = double (imW);

                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //     draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                //printf ("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);
                //  if (lastObject==1)
                //      draggedlocXOffset = -draggedlocXOffset;//-
                draggedlocXOffset -= (locX->getValue() / 2000. * horiz);
            } else if (lastObject == 1) {

                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double horiz = double (imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //     draggedPoint.setFromCartesian(centerPos, currPos);
                // printf ("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

                if (lastObject == 1) {
                    draggedlocXOffset = -draggedlocXOffset;    //-
                }

                draggedlocXOffset -= (locXL->getValue() / 2000. * horiz);
            }

        } else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {

            if (lastObject == 0 || lastObject == 1) {
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double horiz = double (imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //    draggedPoint.setFromCartesian(centerPos, currPos);
                //printf ("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

                if (lastObject == 1) {
                    draggedlocXOffset = -draggedlocXOffset;    //-
                }

                draggedlocXOffset -= (locX->getValue() / 2000. * horiz);
            }
        }

        /*  else if(Smethod->get_active_row_number()==2) {
                if (lastObject==0 || lastObject==1 || lastObject==2 || lastObject==3) {
                if (lastObject==2 || lastObject==3) {
                    // Dragging a line to change the angle
                    PolarCoord draggedPoint;
                    Coord currPos;
                    currPos = provider->posImage;
                    Coord centerPos = draggedCenter;
                    double verti = double(imH);
                    // trick to get the correct angle (clockwise/counter-clockwise)
                    int p = centerPos.y;
                    centerPos.y = currPos.y;
                    currPos.y = p;

                    draggedPoint.setFromCartesian(centerPos, currPos);
                    // compute the projected value of the dragged point
                    draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue())/180.*M_PI);
                    if (lastObject==3)
                        draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locY->getValue() / 200. * verti);
                }


                if (lastObject==0 || lastObject==1) {
                    PolarCoord draggedPoint;
                    Coord currPos;
                    currPos = provider->posImage;
                    Coord centerPos = draggedCenter;
                    double horiz = double(imW);
                    int p = centerPos.y;
                    centerPos.y = currPos.y;
                    currPos.y = p;
                    draggedPoint.setFromCartesian(centerPos, currPos);
                    printf("rad=%f ang=%f\n",draggedPoint.radius,draggedPoint.angle-degree->getValue());
                    draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue()+90.)/180.*M_PI);
                    if (lastObject==1)
                        draggedlocXOffset = -draggedlocXOffset;//-
                    draggedlocXOffset -= (locX->getValue() / 200. * horiz);
                }

                }
            }
            */
        //    EditSubscriber::dragging = true;
        EditSubscriber::action = ES_ACTION_DRAGGING;
        return false;
    } else {
        // this will let this class ignore further drag events
        if (lastObject > -1) { // should theoretically always be true
            if (lastObject == 2 || lastObject == 3) {
                EditSubscriber::visibleGeometry.at(2)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at(3)->state = Geometry::NORMAL;
            }

            if (lastObject == 0 || lastObject == 1) {
                EditSubscriber::visibleGeometry.at(0)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at(1)->state = Geometry::NORMAL;

            } else {
                EditSubscriber::visibleGeometry.at(4)->state = Geometry::NORMAL;
//               EditSubscriber::visibleGeometry.at (lastObject)->state = Geometry::NORMAL;
            }
        }

        lastObject = -1;
        return true;
    }
}

bool Locallab::button1Released()
{
    draggedPointOldAngle = -1000.;
    EditSubscriber::action = ES_ACTION_NONE;

    return true;
}

bool Locallab::drag1(int modifierKey)
{
    // compute the polar coordinate of the mouse position
    EditDataProvider *provider = getEditProvider();
    int imW, imH;
    provider->getImageSize(imW, imH);
    double halfSizeW = imW / 2.;
    double halfSizeH = imH / 2.;

    if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
        if (lastObject == 2) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double (imH);
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //  draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

            if (lastObject == 2) {
                currDraggedlocYOffset -= draggedlocYOffset;
            }

            //else if (lastObject==3)
            // Dragging the lower locY bar
            //  currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            currDraggedlocYOffset = currDraggedlocYOffset * 2000. / verti;

            if (int (currDraggedlocYOffset) != locYT->getIntValue()) {
                locYT->setValue((int (currDraggedlocYOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry(centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged(EvlocallablocY, locYT->getTextValue());
                }

                return true;
            }
        } else if (lastObject == 3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double (imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //  draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

            //  if (lastObject==2)
            // Dragging the upper locY bar
            //      currDraggedlocYOffset -= draggedlocYOffset;
            //  else
            if (lastObject == 3)
                // Dragging the lower locY bar
            {
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            }

            currDraggedlocYOffset = currDraggedlocYOffset * 2000. / verti;

            if (int (currDraggedlocYOffset) != locY->getIntValue()) {

                locY->setValue((int (currDraggedlocYOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();

                updateGeometry(centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged(EvlocallablocY, locY->getTextValue());
                }

                return true;
            }
        }

    } else if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
        if (lastObject == 2 || lastObject == 3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double (imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //   draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

            if (lastObject == 2)
                // Dragging the upper locY bar
            {
                currDraggedlocYOffset -= draggedlocYOffset;
            } else if (lastObject == 3)
                // Dragging the lower locY bar
            {
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            }

            currDraggedlocYOffset = currDraggedlocYOffset * 2000. / verti;

            if (int (currDraggedlocYOffset) != locY->getIntValue()) {
                locY->setValue((int (currDraggedlocYOffset)));
                //Smethod->get_active_row_number()==2
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();

                updateGeometry(centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(),  locYT->getValue(), locXL->getValue());

                if (listener) {
                    if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
                        listener->panelChanged(EvlocallablocY, locY->getTextValue());
                    }

                    //  else listener->panelChanged (EvlocallablocY, locX->getTextValue());

                }

                return true;
            }
        }

    }

    if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
        //else if (lastObject==0) {
        if (lastObject == 0) {// >=4
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double (imW);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //    draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

            if (lastObject == 0) //>=4
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 2000. / horiz;

            if (int (currDraggedStrOffset) != locX->getIntValue()) {
                locX->setValue((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry(centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged(EvlocallablocX, locX->getTextValue());
                }

                return true;
            }
        } else if (lastObject == 1) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double (imW);
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

            if (lastObject == 0)
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 2000. / horiz;

            if (int (currDraggedStrOffset) != locXL->getIntValue()) {
                locXL->setValue((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry(centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged(EvlocallablocX, locX->getTextValue());
                }

                return true;
            }
        }

    } else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
        if (lastObject == 0 || lastObject == 1) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double (imW);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            // draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

            if (lastObject == 0)
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 2000. / horiz;

            if (int (currDraggedStrOffset) != locX->getIntValue()) {
                locX->setValue((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry(centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged(EvlocallablocX, locX->getTextValue());
                }

                return true;
            }
        }
    }

    /*  else if(Smethod->get_active_row_number()==2) {
            if (lastObject==0 || lastObject==1 || lastObject==2 || lastObject==3) {
        if (lastObject==2 || lastObject==3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            Coord currPos;
            currPos = provider->posImage+provider->deltaImage;
            Coord centerPos = draggedCenter;
            double verti = double(imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue())/180.*M_PI);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue() +90.)/180.*M_PI);

            if (lastObject==2)
                currDraggedlocYOffset -= draggedlocYOffset;
            else if (lastObject==3)
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            currDraggedlocYOffset = currDraggedlocYOffset * 200. / verti;
        //  if (int(currDraggedlocYOffset) != locY->getIntValue()) {
        //      locY->setValue((int(currDraggedlocYOffset)));
            if (int(currDraggedlocYOffset) != locX->getIntValue()) {//locX
        //  if (int(currDraggedStrOffset) != locX->getIntValue()) {//locX
                locX->setValue((int(currDraggedlocYOffset)));
                double centX,centY;
                centX=centerX->getValue();
                centY=centerY->getValue();

            //  updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(),  locYT->getValue(), locXL->getValue());
                updateGeometry (centX, centY, locX->getValue(), degree->getValue(), locX->getValue(),  locX->getValue(), locX->getValue());
                if (listener) {
                    if(Smethod->get_active_row_number()==1) listener->panelChanged (EvlocallablocY, locY->getTextValue());

                    }
                return true;
            }
        }
            if (lastObject==0 || lastObject==1) {
                // Dragging the upper or lower locY bar
                PolarCoord draggedPoint;
                Coord currPos;
                currPos = provider->posImage+provider->deltaImage;
                Coord centerPos = draggedCenter;
                double horiz = double(imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint.setFromCartesian(centerPos, currPos);
                double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue() +90.)/180.*M_PI);
                if (lastObject==0)
                    currDraggedStrOffset -= draggedlocXOffset;
                else if (lastObject==1)
                    currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;//-
                    currDraggedStrOffset = currDraggedStrOffset * 200. / horiz;

                if (int(currDraggedStrOffset) != locX->getIntValue()) {
                    locX->setValue((int(currDraggedStrOffset)));
                    double centX,centY;
                    centX=centerX->getValue();
                    centY=centerY->getValue();
                    updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(),locXL->getValue());
                    if (listener)
                        listener->panelChanged (EvlocallablocX, locX->getTextValue());
                    return true;
                }
            }


            }
        }
        */
    //else if (lastObject==4) {
    if (lastObject == 4) {

        // Dragging the circle to change the center
        rtengine::Coord currPos;
        draggedCenter += provider->deltaPrevImage;
        currPos = draggedCenter;
        currPos.clip(imW, imH);
        int newCenterX = int ((double (currPos.x) - halfSizeW) / halfSizeW * 1000.);
        int newCenterY = int ((double (currPos.y) - halfSizeH) / halfSizeH * 1000.);

        if (newCenterX != centerX->getIntValue() || newCenterY != centerY->getIntValue()) {
            centerX->setValue(newCenterX);
            centerY->setValue(newCenterY);
            updateGeometry(newCenterX, newCenterY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

            if (listener) {
                listener->panelChanged(EvlocallabCenter, Glib::ustring::compose("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
            }

            return true;
        }
    }

    return false;
}

void Locallab::switchOffEditMode()
{
    if (edit->get_active()) {
        // switching off the toggle button
        bool wasBlocked = editConn.block(true);
        edit->set_active(false);

        if (!wasBlocked) {
            editConn.block(false);
        }
    }

    EditSubscriber::switchOffEditMode();  // disconnect
}

