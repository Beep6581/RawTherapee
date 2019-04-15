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
 *  2017 Jacques Desmis <jdesmis@gmail.com>
 *  2018 Pierre Cabrera <pierre.cab@gmail.com>
 */

#include "locallab.h"
#include "rtimage.h"
#include <iomanip>
#include "../rtengine/rt_math.h"
#include "options.h"
#include <cmath>
#include "editcallbacks.h"
#include "guiutils.h"
#include <string>
#include <unistd.h>
#include "../rtengine/improcfun.h"
#include "labgrid.h"

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
    expshadhigh(new MyExpander(true, M("TP_LOCALLAB_SHADHIGH"))),
    expvibrance(new MyExpander(true, M("TP_LOCALLAB_VIBRANCE"))),
    expsoft(new MyExpander(true, M("TP_LOCALLAB_SOFT"))),
    expblur(new MyExpander(true, M("TP_LOCALLAB_BLUFR"))),
    exptonemap(new MyExpander(true, M("TP_LOCALLAB_TM"))),
    expreti(new MyExpander(true, new Gtk::HBox())),
    expsharp(new MyExpander(true, new Gtk::HBox())),
    expcontrast(new MyExpander(true, M("TP_LOCALLAB_LOC_CONTRAST"))),
    expcbdl(new MyExpander(true, new Gtk::HBox())),
    expdenoi(new MyExpander(true, new Gtk::HBox())),
    expmaskcol(new MyExpander(false, M("TP_LOCALLAB_SHOW"))),
    expmaskexp(new MyExpander(false, M("TP_LOCALLAB_SHOW"))),
    expmasksh(new MyExpander(false, M("TP_LOCALLAB_SHOW"))),
    expmaskcb(new MyExpander(false, M("TP_LOCALLAB_SHOW"))),

    // CurveEditorGroup widgets
    // Color & Light
    llCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_LUM"))),
    HCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_HLH"))),
    maskCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    // Exposure
    curveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVEEDITOR_TONES_LABEL"))),
    maskexpCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    //Shadows Highlight
    maskSHCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    // Vibrance
    curveEditorGG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_LABEL"))),
    // Retinex
    LocalcurveEditorgainT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAIN"))),
    //CBDL
    maskcbCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    
    // Adjuster widgets
    // Color & Light
    lightness(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LIGHTNESS"), -100, 100, 1, 0))),
    contrast(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0))),
    chroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMA"), -100, 150, 1, 0))),
    strengthgrid(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRGRID"), 0, 100, 1, 30))),
    sensi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    structcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUCCOL"), 0, 100, 1, 0))),
    blurcolde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    blendmaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 10.))),
    chromaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    softradiuscol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    // Exposure
    expcomp(Gtk::manage(new Adjuster(M("TP_EXPOSURE_EXPCOMP"), -2.0, 4.0, 0.05, 0.0))),
    hlcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 0))),
    hlcomprthresh(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 0))),
    black(Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0))),
    shcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50))),
    expchroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_EXPCHROMA"), -50, 100, 1, 30))),
    warm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WARM"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-orange-small.png"))))),
    sensiex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    structexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUCCOL"), 0, 100, 1, 0))),
    blurexpde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    blendmaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 10.))),
    chromaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    softradiusexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    //Shadow hightlights
    highlights(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), 0, 100, 1, 0))),
    h_tonalwidth(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_HLTONALW"), 10, 100, 1, 70))),
    shadows(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_SHADOWS"), 0, 100, 1, 0))),
    s_tonalwidth(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_SHTONALW"), 10, 100, 1, 30))),
    sh_radius(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_RADIUS"), 0, 100, 1, 40))),
    sensihs(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    blendmaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 10.))),
    blurSHde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    chromaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    // Vibrance
    saturated(Gtk::manage(new Adjuster(M("TP_VIBRANCE_SATURATED"), -100., 100., 1., 0.))),
    pastels(Gtk::manage(new Adjuster(M("TP_VIBRANCE_PASTELS"), -100., 100., 1., 0.))),
    sensiv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    //Soft Light
    streng(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENG"), 1, 100, 1, 1))),
    sensisf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 1, 100, 1, 15))),
    // Blur & Noise
    radius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADIUS"), 1.0, 100.0, 0.1, 1.0))),
    strength(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENGTH"), 0, 100, 1, 0))),
    sensibn(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIBN"), 0, 100, 1, 40))),
    // Tone Mapping
    stren(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STREN"), -50, 100, 1, 0))),
    gamma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAM"), 80, 150, 1, 100))),
    estop(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ESTOP"), 10, 400, 1, 140))),
    scaltm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALTM"), 1, 100, 1, 10))),
    rewei(Gtk::manage(new Adjuster(M("TP_LOCALLAB_REWEI"), 0, 9, 1, 0))),
    sensitm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    // Retinex
    str(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STR"), 0, 100, 1, 0))),
    chrrt(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHRRT"), 0, 100, 1, 0))),
    neigh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NEIGH"), 14, 200, 1, 50))),//14, 150, 1, 50))),
    vart(Gtk::manage(new Adjuster(M("TP_LOCALLAB_VART"), 50, 500, 1, 200))),
    dehaz(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DEHAZ"), 0, 100, 1, 0))),
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIH"), 0, 100, 1, 15))),
    softradiusret(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    // Sharpening
    sharcontrast(Gtk::manage(new Adjuster(M("TP_SHARPENING_CONTRAST"), 0, 200, 1, 20))),
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 0.4, 2.5, 0.01, 0.75))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 100))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 0))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sharblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARBLUR"), 0.2, 2.0, 0.05, 0.2))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    // Local Contrast
    lcradius(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_RADIUS"), 20, 200, 1, 80))),
    lcamount(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_AMOUNT"), 0, 1.0, 0.01, 0))),
    lcdarkness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_DARKNESS"), 0, 3.0, 0.01, 1.0))),
    lclightness(Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_LIGHTNESS"), 0, 3.0, 0.01, 1.0))),
    sensilc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    // Contrast by detail levels
    chromacbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMACBDL"), 0, 300, 1, 0))),
    threshold(Gtk::manage(new Adjuster(M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 1., 0.01, 0.2))),
    clarityml(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CLARITYML"), 0, 100, 1, 0))),
    contresid(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRESID"), -100, 100, 1, 0))),
    blurcbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURCBDL"), 0, 100, 1, 0))),
    sensicb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSICB"), 0, 100, 1, 15))),
    softradiuscb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    blendmaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 10.))),
    chromaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskcb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    // Denoise
    noiselumf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINE"), MINCHRO, MAXCHRO, 1, 0))),
    noiselumf0(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINEZERO"), MINCHRO, MAXCHRO, 1, 0))),
    noiselumf2(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINETWO"), MINCHRO, MAXCHRO, 1, 0))),
    noiselumc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMCOARSE"), MINCHRO, MAXCHROCC, 1, 0))),
    noiselumdetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMDETAIL"), 0, 100, 1, 0))),
    noiselequal(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELEQUAL"), -2, 10, 1, 7, Gtk::manage(new RTImage("circle-white-small.png")), Gtk::manage(new RTImage("circle-black-small.png"))))),
    noisechrof(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROFINE"), MINCHRO, MAXCHRO, 1, 0))),
    noisechroc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROCOARSE"), MINCHRO, MAXCHROCC, 1, 0))),
    noisechrodetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHRODETAIL"), 0, 100, 1, 0))),
    adjblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ADJ"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-red-small.png"))))),
    bilateral(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BILATERAL"), 0, 100, 1, 0))),
    sensiden(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIDEN"), 0, 100, 1, 20))),

    // ButtonCheck widgets
    // Color & Light
    curvactiv(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CURV")))),
    invers(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    enaColorMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    // Exposure
    enaExpMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    inversex(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    //Shadows Highlight
    enaSHMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    inverssh(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
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
    //CBDL
    enacbMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    
    // ComboBox widgets
    // Color & Light
    qualitycurveMethod(Gtk::manage(new MyComboBoxText())),
    gridMethod(Gtk::manage(new MyComboBoxText())),
    showmaskcolMethod(Gtk::manage(new MyComboBoxText())),
    //Exposure
    showmaskexpMethod(Gtk::manage(new MyComboBoxText())),
    //Shadows Highlight
    showmaskSHMethod(Gtk::manage(new MyComboBoxText())),
    // Blur & Noise
    blurMethod(Gtk::manage(new MyComboBoxText())),
    // Retinex
    retinexMethod(Gtk::manage(new MyComboBoxText())),
    //CBDL
    showmaskcbMethod(Gtk::manage(new MyComboBoxText())),

    // ThresholdAdjuster widgets
    // Vibrance
    psThreshold(Gtk::manage(new ThresholdAdjuster(M("TP_VIBRANCE_PSTHRESHOLD"), -100., 100., 0., M("TP_VIBRANCE_PSTHRESHOLD_WEIGTHING"), 0, 0., 100., 75., M("TP_VIBRANCE_PSTHRESHOLD_SATTHRESH"), 0, this, false))),

    // Other widgets
    labqualcurv(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_QUALCURV_METHOD") + ":"))),
    lumacontrastMinusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")))),
    lumaneutralButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")))),
    lumacontrastPlusButton(Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS")))),
    gridFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABGRID")))),
    residFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_RESID")))),

    // Others
    defparams(nullptr),
    defpedited(nullptr),
    pe(nullptr)
{
    ToolVBox* const panel = Gtk::manage(new ToolVBox());
    bool showtooltip = options.showtooltip;
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
    // expcolor->set_tooltip_text(M("TP_LOCALLAB_EXPCOLOR_TOOLTIP"));
    setExpandAlignProperties (expmaskcol, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    expmaskcol->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expmaskcol));
    expmaskcol->setLevel (2);

    curvactivConn = curvactiv->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::curvactivChanged));
    lightness->setAdjusterListener(this);
    if(showtooltip) lightness->set_tooltip_text(M("TP_LOCALLAB_EXPCOLOR_TOOLTIP"));
    contrast->setAdjusterListener(this);

    chroma->setAdjusterListener(this);

    if(showtooltip) sensi->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensi->setAdjusterListener(this);

    strengthgrid->setAdjusterListener(this);
    structcol->setAdjusterListener(this);
    blurcolde->setAdjusterListener(this);

    blendmaskcol->setAdjusterListener(this);
    radmaskcol->setAdjusterListener(this);
    chromaskcol->setAdjusterListener(this);
    gammaskcol->setAdjusterListener(this);
    slomaskcol->setAdjusterListener(this);
    softradiuscol->setAdjusterListener(this);

    qualitycurveMethod->append(M("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->set_active(0);
    if(showtooltip) qualitycurveMethod->set_tooltip_markup(M("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));
    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::qualitycurveMethodChanged));

    gridMethod->append(M("TP_LOCALLAB_GRIDONE"));
    gridMethod->append(M("TP_LOCALLAB_GRIDTWO"));
    gridMethod->set_active(0);
    gridMethodConn = gridMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::gridMethodChanged));

    llCurveEditorG->setCurveListener(this);
    llshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "L(L)"));
    llshape->setResetCurve(DiagonalCurveType(defSpot.llcurve.at(0)), defSpot.llcurve);
    if(showtooltip) llshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    std::vector<GradientMilestone> mllshape;
    mllshape.push_back(GradientMilestone(0., 0., 0., 0.));
    mllshape.push_back(GradientMilestone(1., 1., 1., 1.));
    llshape->setBottomBarBgGradient(mllshape);
    llshape->setLeftBarBgGradient(mllshape);

    ccshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "C(C)"));
    ccshape->setResetCurve(DiagonalCurveType(defSpot.cccurve.at(0)), defSpot.cccurve);
    if(showtooltip) ccshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    std::vector<GradientMilestone> mccshape;
    mccshape.push_back(GradientMilestone(0., 0., 0., 0.));
    mccshape.push_back(GradientMilestone(1., 1., 1., 1.));
    ccshape->setBottomBarBgGradient(mccshape);
    ccshape->setLeftBarBgGradient(mccshape);

    // llCurveEditorG->newLine();
    llCurveEditorG->curveListComplete();
    HCurveEditorG->setCurveListener(this);

    LHshape = static_cast<FlatCurveEditor*>(HCurveEditorG->addCurve(CT_Flat, "L(H)", nullptr, false, true));
    LHshape->setIdentityValue(0.);
    LHshape->setResetCurve(FlatCurveType(defSpot.LHcurve.at(0)), defSpot.LHcurve);
    if(showtooltip) LHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    LHshape->setCurveColorProvider(this, 1);
    std::vector<GradientMilestone> mLHshape;

    for (int i = 0; i < 7; i++) {
        float x = float (i) * (1.0f / 6.0);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        mLHshape.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    LHshape->setBottomBarBgGradient(mLHshape);

    HHshape = static_cast<FlatCurveEditor*>(HCurveEditorG->addCurve(CT_Flat, "H(H)", nullptr, false, true));
    HHshape->setIdentityValue(0.);
    HHshape->setResetCurve(FlatCurveType(defSpot.HHcurve.at(0)), defSpot.HHcurve);
    if(showtooltip) HHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    HHshape->setCurveColorProvider(this, 1);
    std::vector<GradientMilestone> mHHshape;

    for (int i = 0; i < 7; i++) {
        float x = float (i) * (1.0f / 6.0);

        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        mHHshape.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    HHshape->setBottomBarBgGradient(mHHshape);

    HCurveEditorG->curveListComplete();

    inversConn  = invers->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversChanged));

    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWSTRUC"));
    showmaskcolMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));

    showmaskcolMethod->set_active(0);
    if(showtooltip) showmaskcolMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskcolMethodConn  = showmaskcolMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::showmaskcolMethodChanged));

    enaColorMaskConn = enaColorMask->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::enaColorMaskChanged));

    maskCurveEditorG->setCurveListener(this);

    CCmaskshape = static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false));
    CCmaskshape->setIdentityValue(0.);
    CCmaskshape->setResetCurve(FlatCurveType(defSpot.CCmaskcurve.at(0)), defSpot.CCmaskcurve);
    if(showtooltip) CCmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    CCmaskshape->setBottomBarColorProvider(this, 7);

    LLmaskshape = static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false));
    LLmaskshape->setIdentityValue(0.);
    LLmaskshape->setResetCurve(FlatCurveType(defSpot.LLmaskcurve.at(0)), defSpot.LLmaskcurve);
    if(showtooltip) LLmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    LLmaskshape->setBottomBarBgGradient(mllshape);
    if(showtooltip) LLmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));

    HHmaskshape = static_cast<FlatCurveEditor *>(maskCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true));
    HHmaskshape->setIdentityValue(0.);
    HHmaskshape->setResetCurve(FlatCurveType(defSpot.HHmaskcurve.at(0)), defSpot.HHmaskcurve);
    if(showtooltip) HHmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    HHmaskshape->setCurveColorProvider(this, 6);
    HHmaskshape->setBottomBarColorProvider(this, 6);

    maskCurveEditorG->curveListComplete();
    labgrid = Gtk::manage(new LabGrid(EvLocallabLabGridValue, M("TP_LOCALLAB_LABGRID_VALUES")));

    ToolParamBlock* const colorBox = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const superFrame = Gtk::manage(new Gtk::Frame());
    superFrame->set_label_align(0.025, 0.5);
    superFrame->set_label_widget(*curvactiv);
    ToolParamBlock* const superBox = Gtk::manage(new ToolParamBlock());
    superBox->pack_start(*lightness);
    superBox->pack_start(*contrast);
    superBox->pack_start(*chroma);
    gridFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const gridBox = Gtk::manage(new ToolParamBlock());
    gridBox->pack_start(*labgrid);
    gridBox->pack_start(*gridMethod);
    gridBox->pack_start(*strengthgrid);
    gridFrame->add(*gridBox);
    superBox->pack_start(*gridFrame);

    superFrame->add(*superBox);
    colorBox->pack_start(*superFrame);
    colorBox->pack_start(*sensi);
    colorBox->pack_start(*structcol);
    colorBox->pack_start(*blurcolde);
    colorBox->pack_start(*softradiuscol);
    Gtk::HBox* const qualcurvbox = Gtk::manage(new Gtk::HBox());
    qualcurvbox->pack_start(*labqualcurv, Gtk::PACK_SHRINK, 4);
    qualcurvbox->pack_start(*qualitycurveMethod);
    colorBox->pack_start(*qualcurvbox);
    colorBox->pack_start(*llCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    colorBox->pack_start(*HCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    colorBox->pack_start(*invers);
    ToolParamBlock* const maskcolBox = Gtk::manage(new ToolParamBlock());
    maskcolBox->pack_start(*showmaskcolMethod, Gtk::PACK_SHRINK, 4);
    maskcolBox->pack_start(*enaColorMask, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*maskCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskcolBox->pack_start(*blendmaskcol, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*radmaskcol, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*chromaskcol, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*gammaskcol, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*slomaskcol, Gtk::PACK_SHRINK, 0);
    expmaskcol->add(*maskcolBox);
    colorBox->pack_start(*expmaskcol);

    expcolor->add(*colorBox);
    expcolor->setLevel(2);

    panel->pack_start(*expcolor, false, false);

    // Exposure
    expexpose->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expexpose));
    enableexposeConn = expexpose->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expexpose));
    if(showtooltip) expexpose->set_tooltip_text(M("TP_LOCALLAB_EXPOSURE_TOOLTIP"));

    setExpandAlignProperties (expmaskexp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    expmaskexp->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expmaskexp));
    expmaskexp->setLevel (2);

    expcomp->setAdjusterListener(this);

    hlcompr->setAdjusterListener(this);

    hlcomprthresh->setAdjusterListener(this);

    black->setAdjusterListener(this);

    shcompr->setAdjusterListener(this);
    expchroma->setAdjusterListener(this);

    if(showtooltip) warm->set_tooltip_text(M("TP_LOCALLAB_WARM_TOOLTIP"));
    warm->setAdjusterListener(this);

    if(showtooltip) sensiex->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensiex->setAdjusterListener(this);
    inversexConn  = inversex->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversexChanged));

    structexp->setAdjusterListener(this);

    blurexpde->setAdjusterListener(this);

    blendmaskexp->setAdjusterListener(this);
    radmaskexp->setAdjusterListener(this);
    chromaskexp->setAdjusterListener(this);
    gammaskexp->setAdjusterListener(this);
    slomaskexp->setAdjusterListener(this);
    softradiusexp->setAdjusterListener(this);

    curveEditorG->setCurveListener(this);

    shapeexpos = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, ""));
    shapeexpos->setResetCurve(DiagonalCurveType(defSpot.excurve.at(0)), defSpot.excurve);
    if(showtooltip) shapeexpos->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_TONES_TOOLTIP"));
    std::vector<GradientMilestone> mshapeexpos;
    mshapeexpos.push_back(GradientMilestone(0., 0., 0., 0.));
    mshapeexpos.push_back(GradientMilestone(1., 1., 1., 1.));
    shapeexpos->setBottomBarBgGradient(mshapeexpos);
    shapeexpos->setLeftBarBgGradient(mshapeexpos);

    curveEditorG->curveListComplete();

    enaExpMaskConn = enaExpMask->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::enaExpMaskChanged));

    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWSTRUC"));
    showmaskexpMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));

    showmaskexpMethod->set_active(0);
    if(showtooltip) showmaskexpMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskexpMethodConn  = showmaskexpMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::showmaskexpMethodChanged));

    maskexpCurveEditorG->setCurveListener(this);

    CCmaskexpshape = static_cast<FlatCurveEditor*>(maskexpCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false));
    CCmaskexpshape->setIdentityValue(0.);
    CCmaskexpshape->setResetCurve(FlatCurveType(defSpot.CCmaskexpcurve.at(0)), defSpot.CCmaskexpcurve);
    if(showtooltip) CCmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    CCmaskexpshape->setBottomBarColorProvider(this, 7);

    LLmaskexpshape = static_cast<FlatCurveEditor*>(maskexpCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false));
    LLmaskexpshape->setIdentityValue(0.);
    LLmaskexpshape->setResetCurve(FlatCurveType(defSpot.LLmaskexpcurve.at(0)), defSpot.LLmaskexpcurve);
    if(showtooltip) LLmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    LLmaskexpshape->setBottomBarBgGradient(mllshape);

    HHmaskexpshape = static_cast<FlatCurveEditor *>(maskexpCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true));
    HHmaskexpshape->setIdentityValue(0.);
    HHmaskexpshape->setResetCurve(FlatCurveType(defSpot.HHmaskexpcurve.at(0)), defSpot.HHmaskexpcurve);
    if(showtooltip) HHmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    HHmaskexpshape->setCurveColorProvider(this, 6);
    HHmaskexpshape->setBottomBarColorProvider(this, 6);

    maskexpCurveEditorG->curveListComplete();

    ToolParamBlock* const exposeBox = Gtk::manage(new ToolParamBlock());
    exposeBox->pack_start(*expcomp);
    exposeBox->pack_start(*hlcompr);
    exposeBox->pack_start(*hlcomprthresh);
    exposeBox->pack_start(*black);
    exposeBox->pack_start(*shcompr);
    exposeBox->pack_start(*expchroma);
    exposeBox->pack_start(*warm);
    exposeBox->pack_start(*sensiex);
    exposeBox->pack_start(*structexp);
    exposeBox->pack_start(*blurexpde);
    exposeBox->pack_start(*softradiusexp);
    exposeBox->pack_start(*curveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    exposeBox->pack_start(*inversex);
    ToolParamBlock* const maskexpBox = Gtk::manage(new ToolParamBlock());
    maskexpBox->pack_start(*showmaskexpMethod, Gtk::PACK_SHRINK, 4);
    maskexpBox->pack_start(*enaExpMask, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*maskexpCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskexpBox->pack_start(*blendmaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*radmaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*chromaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*gammaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*slomaskexp, Gtk::PACK_SHRINK, 0);
    expmaskexp->add(*maskexpBox);
    exposeBox->pack_start(*expmaskexp);

    expexpose->add(*exposeBox);
    expexpose->setLevel(2);

    panel->pack_start(*expexpose, false, false);



//shadow highlight
    expshadhigh->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expshadhigh));
    enableshadhighConn = expshadhigh->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expshadhigh));
    if(showtooltip) expshadhigh->set_tooltip_text(M("TP_LOCALLAB_SHADOWHIGHLIGHT_TOOLTIP"));

    setExpandAlignProperties (expmasksh, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    expmasksh->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expmasksh));
    expmasksh->setLevel (2);

    highlights->setAdjusterListener(this);
    h_tonalwidth->setAdjusterListener(this);
    shadows->setAdjusterListener(this);
    s_tonalwidth->setAdjusterListener(this);
    sh_radius->setAdjusterListener(this);
    sensihs->setAdjusterListener(this);
    blendmaskSH->setAdjusterListener(this);
    radmaskSH->setAdjusterListener(this);
    blurSHde->setAdjusterListener(this);
    chromaskSH->setAdjusterListener(this);
    gammaskSH->setAdjusterListener(this);
    slomaskSH->setAdjusterListener(this);

    enaSHMaskConn = enaSHMask->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::enaSHMaskChanged));
    inversshConn  = inverssh->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversshChanged));

    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskSHMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    
//    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWSTRUC"));

    showmaskSHMethod->set_active(0);
    if(showtooltip) showmaskSHMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskSHMethodConn  = showmaskSHMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::showmaskSHMethodChanged));

    maskSHCurveEditorG->setCurveListener(this);

    CCmaskSHshape = static_cast<FlatCurveEditor*>(maskSHCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false));
    CCmaskSHshape->setIdentityValue(0.);
    CCmaskSHshape->setResetCurve(FlatCurveType(defSpot.CCmaskSHcurve.at(0)), defSpot.CCmaskSHcurve);
    if(showtooltip) CCmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    CCmaskSHshape->setBottomBarColorProvider(this, 7);

    LLmaskSHshape = static_cast<FlatCurveEditor*>(maskSHCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false));
    LLmaskSHshape->setIdentityValue(0.);
    LLmaskSHshape->setResetCurve(FlatCurveType(defSpot.LLmaskSHcurve.at(0)), defSpot.LLmaskSHcurve);
    if(showtooltip) LLmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    LLmaskSHshape->setBottomBarBgGradient(mllshape);

    HHmaskSHshape = static_cast<FlatCurveEditor *>(maskSHCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true));
    HHmaskSHshape->setIdentityValue(0.);
    HHmaskSHshape->setResetCurve(FlatCurveType(defSpot.HHmaskSHcurve.at(0)), defSpot.HHmaskSHcurve);
    if(showtooltip) HHmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    HHmaskSHshape->setCurveColorProvider(this, 6);
    HHmaskSHshape->setBottomBarColorProvider(this, 6);

    maskSHCurveEditorG->curveListComplete();

    ToolParamBlock* const shadhighBox = Gtk::manage(new ToolParamBlock());
    shadhighBox->pack_start(*highlights);
    shadhighBox->pack_start(*h_tonalwidth);
    shadhighBox->pack_start(*shadows);
    shadhighBox->pack_start(*s_tonalwidth);
    shadhighBox->pack_start(*sh_radius);
    shadhighBox->pack_start(*sensihs);
    shadhighBox->pack_start(*blurSHde);
    shadhighBox->pack_start(*inverssh);


    ToolParamBlock* const maskSHBox = Gtk::manage(new ToolParamBlock());
    maskSHBox->pack_start(*showmaskSHMethod, Gtk::PACK_SHRINK, 4);
    maskSHBox->pack_start(*enaSHMask, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*maskSHCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskSHBox->pack_start(*blendmaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*radmaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*chromaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*gammaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*slomaskSH, Gtk::PACK_SHRINK, 0);
    expmasksh->add(*maskSHBox);
    shadhighBox->pack_start(*expmasksh);

    expshadhigh->add(*shadhighBox);
    expshadhigh->setLevel(2);

    panel->pack_start(*expshadhigh, false, false);

    // Vibrance
    expvibrance->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expvibrance));
    enablevibranceConn = expvibrance->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expvibrance));

    saturated->setAdjusterListener(this);

    pastels->setAdjusterListener(this);

    if(showtooltip) psThreshold->set_tooltip_markup(M("TP_VIBRANCE_PSTHRESHOLD_TOOLTIP"));
    psThreshold->setAdjusterListener(this);

    pskinsconn = protectSkins->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::protectskins_toggled));

    ashiftconn = avoidColorShift->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::avoidcolorshift_toggled));

    pastsattogconn = pastSatTog->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::pastsattog_toggled));

    sensiv->setAdjusterListener(this);

    curveEditorGG->setCurveListener(this);

    skinTonesCurve = static_cast<DiagonalCurveEditor*>(curveEditorGG->addCurve(CT_Diagonal, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES")));
    if(showtooltip) skinTonesCurve->setTooltip(M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_TOOLTIP"));
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


    // Blur & Noise
    expblur->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expblur));
    enableblurConn = expblur->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expblur));

    radius->setAdjusterListener(this);

    strength->setAdjusterListener(this);

    if(showtooltip) sensibn->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensibn->setAdjusterListener(this);

    blurMethod->append(M("TP_LOCALLAB_BLNORM"));
    blurMethod->append(M("TP_LOCALLAB_BLINV"));
    blurMethod->set_active(0);
    if(showtooltip) blurMethod->set_tooltip_markup(M("TP_LOCALLAB_BLMETHOD_TOOLTIP"));
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

    if(showtooltip) sensitm->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
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
    Gtk::HBox* const retiTitleHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const retiLabel = Gtk::manage(new Gtk::Label());
    retiLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_RETI")) + Glib::ustring("</b>"));
    retiLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    retiTitleHBox->pack_start(*retiLabel, Gtk::PACK_EXPAND_WIDGET, 0);
    RTImage *retiImage = Gtk::manage(new RTImage("one-to-one-small.png"));
    if(showtooltip) retiImage->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
    retiTitleHBox->pack_end(*retiImage, Gtk::PACK_SHRINK, 0);
    expreti->setLabel(retiTitleHBox);
    expreti->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expreti));
    enableretiConn = expreti->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expreti));

    retinexMethod->append(M("TP_RETINEX_LOW"));
    retinexMethod->append(M("TP_RETINEX_UNIFORM"));
    retinexMethod->append(M("TP_RETINEX_HIGH"));
    retinexMethod->set_active(0);
    if(showtooltip) retinexMethod->set_tooltip_markup(M("TP_LOCRETI_METHOD_TOOLTIP"));
    retinexMethodConn = retinexMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::retinexMethodChanged));

    str->setAdjusterListener(this);

    neigh->setAdjusterListener(this);

    vart->setAdjusterListener(this);

    dehaz->setAdjusterListener(this);

    chrrt->setAdjusterListener(this);

    if(showtooltip) sensih->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensih->setAdjusterListener(this);
    softradiusret->setAdjusterListener(this);

    LocalcurveEditorgainT->setCurveListener(this);

    cTgainshape = static_cast<FlatCurveEditor*>(LocalcurveEditorgainT->addCurve(CT_Flat, "", nullptr, false, false));
    cTgainshape->setIdentityValue(0.);
    cTgainshape->setResetCurve(FlatCurveType(defSpot.localTgaincurve.at(0)), defSpot.localTgaincurve);
    if(showtooltip) cTgainshape->setTooltip(M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));

    LocalcurveEditorgainT->curveListComplete();

    inversretConn  = inversret->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversretChanged));

    ToolParamBlock* const retiBox = Gtk::manage(new ToolParamBlock());
    retiBox->pack_start(*retinexMethod);
    retiBox->pack_start(*str);
    retiBox->pack_start(*chrrt);
    retiBox->pack_start(*neigh);
    retiBox->pack_start(*vart);
    retiBox->pack_start(*dehaz);
    retiBox->pack_start(*softradiusret);
    retiBox->pack_start(*sensih);
    retiBox->pack_start(*LocalcurveEditorgainT, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    retiBox->pack_start(*inversret);
    expreti->add(*retiBox);
    expreti->setLevel(2);

    panel->pack_start(*expreti, false, false);

    // Sharpening
    Gtk::HBox* const sharpTitleHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const sharpLabel = Gtk::manage(new Gtk::Label());
    sharpLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_SHARP")) + Glib::ustring("</b>"));
    sharpLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    sharpTitleHBox->pack_start(*sharpLabel, Gtk::PACK_EXPAND_WIDGET, 0);
    RTImage *sharpImage = Gtk::manage(new RTImage("one-to-one-small.png"));
    if(showtooltip) sharpImage->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
    sharpTitleHBox->pack_end(*sharpImage, Gtk::PACK_SHRINK, 0);
    expsharp->setLabel(sharpTitleHBox);
    expsharp->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsharp));
    enablesharpConn = expsharp->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expsharp));

    sharcontrast->setAdjusterListener(this);

    sharradius->setAdjusterListener(this);

    sharamount->setAdjusterListener(this);

    shardamping->setAdjusterListener(this);

    shariter->setAdjusterListener(this);

    sharblur->setAdjusterListener(this);

    if(showtooltip) sensisha->set_tooltip_text(M("TP_LOCALLAB_SENSIS_TOOLTIP"));
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
    Gtk::HBox* const cbdlTitleHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const cbdlLabel = Gtk::manage(new Gtk::Label());
    cbdlLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_CBDL")) + Glib::ustring("</b>"));
    cbdlLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    cbdlTitleHBox->pack_start(*cbdlLabel, Gtk::PACK_EXPAND_WIDGET, 0);
    RTImage *cbdlImage = Gtk::manage(new RTImage("one-to-one-small.png"));
    if(showtooltip) cbdlImage->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
    cbdlTitleHBox->pack_end(*cbdlImage, Gtk::PACK_SHRINK, 0);
    expcbdl->setLabel(cbdlTitleHBox);
    expcbdl->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcbdl));
    enablecbdlConn = expcbdl->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcbdl));
    if(showtooltip) expcbdl->set_tooltip_text(M("TP_LOCALLAB_EXPCBDL_TOOLTIP"));

    setExpandAlignProperties (expmaskcb, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    expmaskcb->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expmaskcb));
    expmaskcb->setLevel (2);

    for (int i = 0; i < 5; i++) {
        Glib::ustring ss;
        ss = Glib::ustring::format(i);

        if (i == 0) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMAFINEST"));
        } else if (i == 4) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        }

        multiplier[i] = Gtk::manage(new Adjuster(ss, 0.0, 4.0, 0.01, 1.0));
        multiplier[i]->setAdjusterListener(this);
    }

    if(showtooltip) chromacbdl->set_tooltip_text(M("TP_LOCALLAB_CHROMACB_TOOLTIP"));
    chromacbdl->setAdjusterListener(this);

    threshold->setAdjusterListener(this);

    if(showtooltip) sensicb->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensicb->setAdjusterListener(this);
    softradiuscb->setAdjusterListener(this);
    clarityml->setAdjusterListener(this);
    contresid->setAdjusterListener(this);
    blurcbdl->setAdjusterListener(this);
    blendmaskcb->setAdjusterListener(this);
    radmaskcb->setAdjusterListener(this);
    chromaskcb->setAdjusterListener(this);
    gammaskcb->setAdjusterListener(this);
    slomaskcb->setAdjusterListener(this);

    enacbMaskConn = enacbMask->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::enacbMaskChanged));

    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskcbMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcbMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    

    showmaskcbMethod->set_active(0);
    if(showtooltip) showmaskcbMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskcbMethodConn  = showmaskcbMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::showmaskcbMethodChanged));

    maskcbCurveEditorG->setCurveListener(this);

    CCmaskcbshape = static_cast<FlatCurveEditor*>(maskcbCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false));
    CCmaskcbshape->setIdentityValue(0.);
    CCmaskcbshape->setResetCurve(FlatCurveType(defSpot.CCmaskcbcurve.at(0)), defSpot.CCmaskcbcurve);
    if(showtooltip) CCmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    CCmaskcbshape->setBottomBarColorProvider(this, 7);

    LLmaskcbshape = static_cast<FlatCurveEditor*>(maskcbCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false));
    LLmaskcbshape->setIdentityValue(0.);
    LLmaskcbshape->setResetCurve(FlatCurveType(defSpot.LLmaskcbcurve.at(0)), defSpot.LLmaskcbcurve);
    if(showtooltip) LLmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    LLmaskcbshape->setBottomBarBgGradient(mllshape);

    HHmaskcbshape = static_cast<FlatCurveEditor *>(maskcbCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true));
    HHmaskcbshape->setIdentityValue(0.);
    HHmaskcbshape->setResetCurve(FlatCurveType(defSpot.HHmaskcbcurve.at(0)), defSpot.HHmaskcbcurve);
    if(showtooltip) HHmaskcbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    HHmaskcbshape->setCurveColorProvider(this, 6);
    HHmaskcbshape->setBottomBarColorProvider(this, 6);

    maskcbCurveEditorG->curveListComplete();



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

    ToolParamBlock* const maskcbBox = Gtk::manage(new ToolParamBlock());
    maskcbBox->pack_start(*showmaskcbMethod, Gtk::PACK_SHRINK, 4);
    maskcbBox->pack_start(*enacbMask, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*maskcbCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskcbBox->pack_start(*blendmaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*radmaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*chromaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*gammaskcb, Gtk::PACK_SHRINK, 0);
    maskcbBox->pack_start(*slomaskcb, Gtk::PACK_SHRINK, 0);
    expmaskcb->add(*maskcbBox);

    Gtk::HSeparator *separator = Gtk::manage(new  Gtk::HSeparator());
    cbdlBox->pack_start(*separator, Gtk::PACK_SHRINK, 2);
    cbdlBox->pack_start(*chromacbdl);
    cbdlBox->pack_start(*threshold);
    cbdlBox->pack_start(*blurcbdl);
    residFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const residBox = Gtk::manage(new ToolParamBlock());    
    residBox->pack_start(*clarityml);
    residBox->pack_start(*contresid);
    residFrame->add(*residBox);
    cbdlBox->pack_start(*residFrame);
    cbdlBox->pack_start(*softradiuscb);
    cbdlBox->pack_start(*sensicb);
    cbdlBox->pack_start(*expmaskcb);
    expcbdl->add(*cbdlBox);
    expcbdl->setLevel(2);

    panel->pack_start(*expcbdl, false, false);

    // Denoise
    Gtk::HBox* const denoiTitleHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const denoiLabel = Gtk::manage(new Gtk::Label());
    denoiLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(M("TP_LOCALLAB_DENOIS")) + Glib::ustring("</b>"));
    denoiLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    denoiTitleHBox->pack_start(*denoiLabel, Gtk::PACK_EXPAND_WIDGET, 0);
    RTImage *denoiImage = Gtk::manage(new RTImage("one-to-one-small.png"));
    if(showtooltip) denoiImage->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
    denoiTitleHBox->pack_end(*denoiImage, Gtk::PACK_SHRINK, 0);
    expdenoi->setLabel(denoiTitleHBox);
    expdenoi->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expdenoi));
    enabledenoiConn = expdenoi->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expdenoi));

    noiselumf->setAdjusterListener(this);
    noiselumf0->setAdjusterListener(this);
    noiselumf2->setAdjusterListener(this);

    if(showtooltip) noiselumc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    noiselumc->setAdjusterListener(this);

    noiselumdetail->setAdjusterListener(this);

    noiselequal->setAdjusterListener(this);

    noisechrof->setAdjusterListener(this);

    if(showtooltip) noisechroc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    noisechroc->setAdjusterListener(this);

    noisechrodetail->setAdjusterListener(this);

    adjblur->setAdjusterListener(this);

    bilateral->setAdjusterListener(this);

    sensiden->setAdjusterListener(this);

    ToolParamBlock* const denoisBox = Gtk::manage(new ToolParamBlock());
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
    denoisBox->pack_start(*wavFrame);
    denoisBox->pack_start(*bilateral);
    denoisBox->pack_start(*sensiden);
    expdenoi->add(*denoisBox);
    expdenoi->setLevel(2);

    panel->pack_start(*expdenoi, false, false);

    pack_start(*panel);

    setParamEditable(false);

    show_all();
}

Locallab::~Locallab()
{
    idle_register.destroy();

    delete llCurveEditorG;
    delete HCurveEditorG;
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
        expshadhigh->set_expanded(expshadhigh == expander);
        expvibrance->set_expanded(expvibrance == expander);
        expsoft->set_expanded(expsoft == expander);
        expblur->set_expanded(expblur == expander);
        exptonemap->set_expanded(exptonemap == expander);
        expreti->set_expanded(expreti == expander);
        expsharp->set_expanded(expsharp == expander);
        expcontrast->set_expanded(expcontrast == expander);
        expcbdl->set_expanded(expcbdl == expander);
        expdenoi->set_expanded(expdenoi == expander);
        expmaskcol->set_expanded(expmaskcol == expander);
        expmaskexp->set_expanded(expmaskexp == expander);
        expmasksh->set_expanded(expmasksh == expander);
        expmaskcb->set_expanded(expmaskcb == expander);

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
    } else if (expander == expshadhigh) {
        event = EvLocenashadhigh;
        expConn = &enableshadhighConn;
    } else if (expander == expvibrance) {
        event = EvLocenavibrance;
        expConn = &enablevibranceConn;
    } else if (expander == expsoft) {
        event = EvLocenasoft;
        expConn = &enablesoftConn;
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
    tpOpen.push_back(expshadhigh->get_expanded());
    tpOpen.push_back(expvibrance->get_expanded());
    tpOpen.push_back(expsoft->get_expanded());
    tpOpen.push_back(expblur->get_expanded());
    tpOpen.push_back(exptonemap->get_expanded());
    tpOpen.push_back(expreti->get_expanded());
    tpOpen.push_back(expsharp->get_expanded());
    tpOpen.push_back(expcontrast->get_expanded());
    tpOpen.push_back(expcbdl->get_expanded());
    tpOpen.push_back(expdenoi->get_expanded());
    tpOpen.push_back(expmaskcol->get_expanded());
    tpOpen.push_back(expmaskexp->get_expanded());
    tpOpen.push_back(expmasksh->get_expanded());
    tpOpen.push_back(expmaskcb->get_expanded());

}

void Locallab::refChanged(double huer, double lumar, double chromar)
{
    if (!batchMode) {
        // Hue reference normalization (between 0 and 1)
        double normHuer = huer;
        float h = Color::huelab_to_huehsv2(normHuer);
        h += 1.f / 6.f;

        if (h > 1.f) {
            h -= 1.f;
        }

        normHuer = h;

        // Luma reference normalization (between 0 and 1)
        double normLumar = lumar / 100.f;

        // Chroma reference normalization (between 0 and 1)
        double normChromar = chromar / 137.4f;

        // printf("nh=%f nl=%f nc=%f\n", normHuer, normLumar, normChromar);

        idle_register.add(
        [this, normHuer, normLumar, normChromar]() -> bool {
            GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

            // Update Color & Light mask background
            CCmaskshape->updateLocallabBackground(normChromar);
            LLmaskshape->updateLocallabBackground(normLumar);
            HHmaskshape->updateLocallabBackground(normHuer);

            // Update Exposure mask background
            CCmaskexpshape->updateLocallabBackground(normChromar);
            LLmaskexpshape->updateLocallabBackground(normLumar);
            HHmaskexpshape->updateLocallabBackground(normHuer);

            // Update Shadow Highlight mask background
            CCmaskSHshape->updateLocallabBackground(normChromar);
            LLmaskSHshape->updateLocallabBackground(normLumar);
            HHmaskSHshape->updateLocallabBackground(normHuer);

            // Update CBDL mask background
            CCmaskcbshape->updateLocallabBackground(normChromar);
            LLmaskcbshape->updateLocallabBackground(normLumar);
            HHmaskcbshape->updateLocallabBackground(normHuer);

            return false;
        }
        );
    }
}

void Locallab::updateToolState(std::vector<int> &tpOpen)
{
    if (tpOpen.size() >= 18) {
        expsettings->setExpanded(tpOpen.at(0));
        expcolor->set_expanded(tpOpen.at(1));
        expexpose->set_expanded(tpOpen.at(2));
        expshadhigh->set_expanded(tpOpen.at(3));
        expvibrance->set_expanded(tpOpen.at(4));
        expsoft->set_expanded(tpOpen.at(5));
        expblur->set_expanded(tpOpen.at(6));
        exptonemap->set_expanded(tpOpen.at(7));
        expreti->set_expanded(tpOpen.at(8));
        expsharp->set_expanded(tpOpen.at(9));
        expcontrast->set_expanded(tpOpen.at(10));
        expcbdl->set_expanded(tpOpen.at(11));
        expdenoi->set_expanded(tpOpen.at(12));
        expmaskcol->set_expanded(tpOpen.at(13));
        expmaskexp->set_expanded(tpOpen.at(14));
        expmasksh->set_expanded(tpOpen.at(15));
        expmaskcb->set_expanded(tpOpen.at(16));
    }
}

void Locallab::lumaneutralPressed()
{
    // printf("lumaneutralPressed\n");

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue(1.0);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void Locallab::lumacontrastPlusPressed()
{
    // printf("lumacontrastPlusPressed\n");

    for (int i = 0; i < 5; i++) {
        float inc = (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + 0.01f * inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

void Locallab::lumacontrastMinusPressed()
{
    // printf("lumacontrastMinusPressed\n");

    for (int i = 0; i < 5; i++) {
        float inc = - (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + 0.01f * inc);
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
        r->structexclu = pp->locallab.spots.at(i).structexclu;
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

        if (pp->locallab.spots.at(i).qualityMethod == "enh") {
            r->qualityMethod = 0;
        } else {
            r->qualityMethod = 1;
        }

        r->transit = pp->locallab.spots.at(i).transit;
        r->thresh = pp->locallab.spots.at(i).thresh;
        r->iter = pp->locallab.spots.at(i).iter;
        r->balan = pp->locallab.spots.at(i).balan;
        r->transitweak = pp->locallab.spots.at(i).transitweak;
        r->avoid = pp->locallab.spots.at(i).avoid;

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
            r->structexclu = newSpot->structexclu;
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

            if (newSpot->qualityMethod == "enh") {
                r->qualityMethod = 0;
            } else {
                r->qualityMethod = 1;
            }

            r->transit = newSpot->transit;
            r->thresh = newSpot->thresh;
            r->iter = newSpot->iter;
            r->balan = newSpot->balan;
            r->transitweak = newSpot->transitweak;
            r->avoid = newSpot->avoid;
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

                    // Update Locallab tools GUI with selected spot
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
            r->structexclu = newSpot->structexclu;
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

            if (newSpot->qualityMethod == "enh") {
                r->qualityMethod = 0;
            } else {
                r->qualityMethod = 1;
            }

            r->transit = newSpot->transit;
            r->thresh = newSpot->thresh;
            r->iter = newSpot->iter;
            r->balan = newSpot->balan;
            r->transitweak = newSpot->transitweak;
            r->avoid = newSpot->avoid;
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
                    pp->locallab.spots.at(pp->locallab.selspot).structexclu = r->structexclu;
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
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "enh";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "enhden";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).transit = r->transit;
                    pp->locallab.spots.at(pp->locallab.selspot).thresh = r->thresh;
                    pp->locallab.spots.at(pp->locallab.selspot).iter = r->iter;
                    pp->locallab.spots.at(pp->locallab.selspot).balan = r->balan;
                    pp->locallab.spots.at(pp->locallab.selspot).transitweak = r->transitweak;
                    pp->locallab.spots.at(pp->locallab.selspot).avoid = r->avoid;
                    // Color & Light
                    pp->locallab.spots.at(pp->locallab.selspot).expcolor = expcolor->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).curvactiv = curvactiv->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).lightness = lightness->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).contrast = contrast->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).chroma = chroma->getIntValue();
                    labgrid->getParams(pp->locallab.spots.at(pp->locallab.selspot).labgridALow, pp->locallab.spots.at(pp->locallab.selspot).labgridBLow, pp->locallab.spots.at(pp->locallab.selspot).labgridAHigh, pp->locallab.spots.at(pp->locallab.selspot).labgridBHigh);
                    pp->locallab.spots.at(pp->locallab.selspot).strengthgrid = strengthgrid->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).labgridALow *= LocallabParams::LABGRIDL_CORR_MAX;
                    pp->locallab.spots.at(pp->locallab.selspot).labgridAHigh *= LocallabParams::LABGRIDL_CORR_MAX;
                    pp->locallab.spots.at(pp->locallab.selspot).labgridBLow *= LocallabParams::LABGRIDL_CORR_MAX;
                    pp->locallab.spots.at(pp->locallab.selspot).labgridBHigh *= LocallabParams::LABGRIDL_CORR_MAX;
                    pp->locallab.spots.at(pp->locallab.selspot).sensi = sensi->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).structcol = structcol->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).blurcolde = blurcolde->getIntValue();

                    if (qualitycurveMethod->get_active_row_number() == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = "none";
                    } else if (qualitycurveMethod->get_active_row_number() == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = "std";
                    }

                    if (gridMethod->get_active_row_number() == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).gridMethod = "one";
                    } else if (gridMethod->get_active_row_number() == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).gridMethod = "two";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).llcurve = llshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).cccurve = ccshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).LHcurve = LHshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).HHcurve = HHshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).invers = invers->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).enaColorMask = enaColorMask->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).CCmaskcurve = CCmaskshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).LLmaskcurve = LLmaskshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).HHmaskcurve = HHmaskshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).blendmaskcol = blendmaskcol->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).radmaskcol = radmaskcol->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).chromaskcol = chromaskcol->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).gammaskcol = gammaskcol->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).slomaskcol = slomaskcol->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).softradiuscol = softradiuscol->getValue();
                    // Exposure
                    pp->locallab.spots.at(pp->locallab.selspot).expexpose = expexpose->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).expcomp = expcomp->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).hlcompr = hlcompr->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).hlcomprthresh = hlcomprthresh->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).black = black->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shcompr = shcompr->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).expchroma = expchroma->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).warm = warm->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensiex = sensiex->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).structexp = structexp->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).blurexpde = blurexpde->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).excurve = shapeexpos->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).inversex = inversex->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).enaExpMask = enaExpMask->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve = LLmaskexpshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve = CCmaskexpshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).HHmaskexpcurve = HHmaskexpshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).blendmaskexp = blendmaskexp->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).radmaskexp = radmaskexp->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).chromaskexp = chromaskexp->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).gammaskexp = gammaskexp->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).slomaskexp = slomaskexp->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).softradiusexp = softradiusexp->getValue();
                    // Shadow highlight
                    pp->locallab.spots.at(pp->locallab.selspot).expshadhigh = expshadhigh->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).highlights = highlights->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).h_tonalwidth = h_tonalwidth->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shadows = shadows->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).s_tonalwidth = s_tonalwidth->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sh_radius = sh_radius->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensihs = sensihs->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).enaSHMask = enaSHMask->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).LLmaskSHcurve = LLmaskSHshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).CCmaskSHcurve = CCmaskSHshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).HHmaskSHcurve = HHmaskSHshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).blendmaskSH = blendmaskSH->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).radmaskSH = radmaskSH->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).blurSHde = blurSHde->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).inverssh = inverssh->get_active();
                    pp->locallab.spots.at(pp->locallab.selspot).chromaskSH = chromaskSH->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).gammaskSH = gammaskSH->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).slomaskSH = slomaskSH->getValue();
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
                    // Blur & Noise
                    pp->locallab.spots.at(pp->locallab.selspot).expblur = expblur->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).radius = radius->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).strength = strength->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensibn = sensibn->getIntValue();

                    if (blurMethod->get_active_row_number() == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).blurMethod = "norm";
                    } else if (blurMethod->get_active_row_number() == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).blurMethod = "inv";
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
                    pp->locallab.spots.at(pp->locallab.selspot).softradiusret = softradiusret->getValue();
                    // Sharpening
                    pp->locallab.spots.at(pp->locallab.selspot).expsharp = expsharp->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).sharcontrast = sharcontrast->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sharradius = sharradius->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sharamount = sharamount->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shardamping = shardamping->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).shariter = shariter->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sharblur = sharblur->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensisha = sensisha->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).inverssha = inverssha->get_active();
                    // Local Contrast
                    pp->locallab.spots.at(pp->locallab.selspot).expcontrast = expcontrast->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).lcradius = lcradius->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).lcamount = lcamount->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).lcdarkness = lcdarkness->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).lclightness = lclightness->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensilc = sensilc->getIntValue();
                    // Contrast by detail levels
                    pp->locallab.spots.at(pp->locallab.selspot).expcbdl = expcbdl->getEnabled();

                    for (int i = 0; i < 5; i++) {
                        pp->locallab.spots.at(pp->locallab.selspot).mult[i] = multiplier[i]->getValue();
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).chromacbdl = chromacbdl->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).threshold = threshold->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensicb = sensicb->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).clarityml = clarityml->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).contresid = contresid->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).blurcbdl = blurcbdl->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).softradiuscb = softradiuscb->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).enacbMask = enacbMask->get_active();

                    pp->locallab.spots.at(pp->locallab.selspot).LLmaskcbcurve = LLmaskcbshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).CCmaskcbcurve = CCmaskcbshape->getCurve();
                    pp->locallab.spots.at(pp->locallab.selspot).HHmaskcbcurve = HHmaskcbshape->getCurve();

                    pp->locallab.spots.at(pp->locallab.selspot).blendmaskcb = blendmaskcb->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).radmaskcb = radmaskcb->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).chromaskcb = chromaskcb->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).gammaskcb = gammaskcb->getValue();
                    pp->locallab.spots.at(pp->locallab.selspot).slomaskcb = slomaskcb->getValue();
                    
                    // Denoise
                    pp->locallab.spots.at(pp->locallab.selspot).expdenoi = expdenoi->getEnabled();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselumf = noiselumf->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselumf0 = noiselumf0->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselumf2 = noiselumf2->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselumc = noiselumc->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselumdetail = noiselumdetail->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noiselequal = noiselequal->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noisechrof = noisechrof->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noisechroc = noisechroc->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).noisechrodetail = noisechrodetail->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).adjblur = adjblur->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).bilateral = bilateral->getIntValue();
                    pp->locallab.spots.at(pp->locallab.selspot).sensiden = sensiden->getIntValue();
                }

                ControlSpotPanel::SpotEdited* const se = expsettings->getEditedStates();

                if (pe) {
                    if (pp->locallab.selspot < (int)pe->locallab.spots.size()) {
                        pe->locallab.spots.at(pp->locallab.selspot).name = pe->locallab.spots.at(pp->locallab.selspot).name || se->name;
                        pe->locallab.spots.at(pp->locallab.selspot).isvisible = pe->locallab.spots.at(pp->locallab.selspot).isvisible || se->isvisible;
                        pe->locallab.spots.at(pp->locallab.selspot).shape = pe->locallab.spots.at(pp->locallab.selspot).shape || se->shape;
                        pe->locallab.spots.at(pp->locallab.selspot).spotMethod = pe->locallab.spots.at(pp->locallab.selspot).spotMethod || se->spotMethod;
                        pe->locallab.spots.at(pp->locallab.selspot).sensiexclu = pe->locallab.spots.at(pp->locallab.selspot).sensiexclu || se->sensiexclu;
                        pe->locallab.spots.at(pp->locallab.selspot).structexclu = pe->locallab.spots.at(pp->locallab.selspot).structexclu || se->structexclu;
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
                        pe->locallab.spots.at(pp->locallab.selspot).transitweak = pe->locallab.spots.at(pp->locallab.selspot).transitweak || se->transitweak;
                        pe->locallab.spots.at(pp->locallab.selspot).balan = pe->locallab.spots.at(pp->locallab.selspot).balan || se->balan;
                        pe->locallab.spots.at(pp->locallab.selspot).avoid = pe->locallab.spots.at(pp->locallab.selspot).avoid || se->avoid;
                        // Color & Light
                        pe->locallab.spots.at(pp->locallab.selspot).expcolor = pe->locallab.spots.at(pp->locallab.selspot).expcolor || !expcolor->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).curvactiv = pe->locallab.spots.at(pp->locallab.selspot).curvactiv || !curvactiv->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).lightness = pe->locallab.spots.at(pp->locallab.selspot).lightness || lightness->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).contrast = pe->locallab.spots.at(pp->locallab.selspot).contrast || contrast->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).labgridALow = pe->locallab.spots.at(pp->locallab.selspot).labgridBLow = pe->locallab.spots.at(pp->locallab.selspot).labgridAHigh = pe->locallab.spots.at(pp->locallab.selspot).labgridBHigh = labgrid->getEdited();
                        pe->locallab.spots.at(pp->locallab.selspot).strengthgrid = pe->locallab.spots.at(pp->locallab.selspot).strengthgrid || strengthgrid->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).chroma = pe->locallab.spots.at(pp->locallab.selspot).chroma || chroma->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensi = pe->locallab.spots.at(pp->locallab.selspot).sensi || sensi->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).structcol = pe->locallab.spots.at(pp->locallab.selspot).structcol || structcol->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = pe->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod || qualitycurveMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pe->locallab.spots.at(pp->locallab.selspot).gridMethod = pe->locallab.spots.at(pp->locallab.selspot).gridMethod || gridMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pe->locallab.spots.at(pp->locallab.selspot).llcurve = pe->locallab.spots.at(pp->locallab.selspot).llcurve || !llshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).cccurve = pe->locallab.spots.at(pp->locallab.selspot).cccurve || !ccshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).LHcurve = pe->locallab.spots.at(pp->locallab.selspot).LHcurve || !LHshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).HHcurve = pe->locallab.spots.at(pp->locallab.selspot).HHcurve || !HHshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).invers = pe->locallab.spots.at(pp->locallab.selspot).invers || !invers->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).enaColorMask = pe->locallab.spots.at(pp->locallab.selspot).enaColorMask || !enaColorMask->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).CCmaskcurve = pe->locallab.spots.at(pp->locallab.selspot).CCmaskcurve || !CCmaskshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).LLmaskcurve = pe->locallab.spots.at(pp->locallab.selspot).LLmaskcurve || !LLmaskshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).HHmaskcurve = pe->locallab.spots.at(pp->locallab.selspot).HHmaskcurve || !HHmaskshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).blurcolde = pe->locallab.spots.at(pp->locallab.selspot).blurcolde || blurcolde->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).blendmaskcol = pe->locallab.spots.at(pp->locallab.selspot).blendmaskcol || blendmaskcol->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).radmaskcol = pe->locallab.spots.at(pp->locallab.selspot).radmaskcol || radmaskcol->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).chromaskcol = pe->locallab.spots.at(pp->locallab.selspot).chromaskcol || chromaskcol->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).gammaskcol = pe->locallab.spots.at(pp->locallab.selspot).gammaskcol || gammaskcol->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).slomaskcol = pe->locallab.spots.at(pp->locallab.selspot).slomaskcol || slomaskcol->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).softradiuscol = pe->locallab.spots.at(pp->locallab.selspot).softradiuscol || softradiuscol->getEditedState();
                        // Exposure
                        pe->locallab.spots.at(pp->locallab.selspot).expexpose = pe->locallab.spots.at(pp->locallab.selspot).expexpose || !expexpose->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).expcomp = pe->locallab.spots.at(pp->locallab.selspot).expcomp || expcomp->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).hlcompr = pe->locallab.spots.at(pp->locallab.selspot).hlcompr || hlcompr->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).hlcomprthresh = pe->locallab.spots.at(pp->locallab.selspot).hlcomprthresh || hlcomprthresh->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).black = pe->locallab.spots.at(pp->locallab.selspot).black || black->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).shcompr = pe->locallab.spots.at(pp->locallab.selspot).shcompr || shcompr->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).expchroma = pe->locallab.spots.at(pp->locallab.selspot).expchroma || expchroma->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).warm = pe->locallab.spots.at(pp->locallab.selspot).warm || warm->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensiex = pe->locallab.spots.at(pp->locallab.selspot).sensiex || sensiex->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).structexp = pe->locallab.spots.at(pp->locallab.selspot).structexp || structexp->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).blurexpde = pe->locallab.spots.at(pp->locallab.selspot).blurexpde || blurexpde->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).excurve = pe->locallab.spots.at(pp->locallab.selspot).excurve || !shapeexpos->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).inversex = pe->locallab.spots.at(pp->locallab.selspot).inversex || !inversex->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).enaExpMask = pe->locallab.spots.at(pp->locallab.selspot).enaExpMask || !enaExpMask->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve = pe->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve || !CCmaskexpshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve = pe->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve || !LLmaskexpshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).HHmaskexpcurve = pe->locallab.spots.at(pp->locallab.selspot).HHmaskexpcurve || !HHmaskexpshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).blendmaskexp = pe->locallab.spots.at(pp->locallab.selspot).blendmaskexp || blendmaskexp->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).radmaskexp = pe->locallab.spots.at(pp->locallab.selspot).radmaskexp || radmaskexp->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).chromaskexp = pe->locallab.spots.at(pp->locallab.selspot).chromaskexp || chromaskexp->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).gammaskexp = pe->locallab.spots.at(pp->locallab.selspot).gammaskexp || gammaskexp->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).slomaskexp = pe->locallab.spots.at(pp->locallab.selspot).slomaskexp || slomaskexp->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).softradiusexp = pe->locallab.spots.at(pp->locallab.selspot).softradiusexp || softradiusexp->getEditedState();
                        // Shadow highlight
                        pe->locallab.spots.at(pp->locallab.selspot).expshadhigh = pe->locallab.spots.at(pp->locallab.selspot).expshadhigh || !expshadhigh->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).highlights = pe->locallab.spots.at(pp->locallab.selspot).highlights || highlights->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).h_tonalwidth = pe->locallab.spots.at(pp->locallab.selspot).h_tonalwidth || h_tonalwidth->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).shadows = pe->locallab.spots.at(pp->locallab.selspot).shadows || shadows->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).s_tonalwidth = pe->locallab.spots.at(pp->locallab.selspot).s_tonalwidth || s_tonalwidth->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sh_radius = pe->locallab.spots.at(pp->locallab.selspot).sh_radius || sh_radius->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensihs = pe->locallab.spots.at(pp->locallab.selspot).sensihs || sensihs->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).enaSHMask = pe->locallab.spots.at(pp->locallab.selspot).enaSHMask || !enaSHMask->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).CCmaskSHcurve = pe->locallab.spots.at(pp->locallab.selspot).CCmaskSHcurve || !CCmaskSHshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).LLmaskSHcurve = pe->locallab.spots.at(pp->locallab.selspot).LLmaskSHcurve || !LLmaskSHshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).HHmaskSHcurve = pe->locallab.spots.at(pp->locallab.selspot).HHmaskSHcurve || !HHmaskSHshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).blendmaskSH = pe->locallab.spots.at(pp->locallab.selspot).blendmaskSH || blendmaskSH->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).radmaskSH = pe->locallab.spots.at(pp->locallab.selspot).radmaskSH || radmaskSH->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).blurSHde = pe->locallab.spots.at(pp->locallab.selspot).blurSHde || blurSHde->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).inverssh = pe->locallab.spots.at(pp->locallab.selspot).inverssh || !inverssh->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).chromaskSH = pe->locallab.spots.at(pp->locallab.selspot).chromaskSH || chromaskSH->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).gammaskSH = pe->locallab.spots.at(pp->locallab.selspot).gammaskSH || gammaskSH->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).slomaskSH = pe->locallab.spots.at(pp->locallab.selspot).slomaskSH || slomaskSH->getEditedState();
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
                        pe->locallab.spots.at(pp->locallab.selspot).softradiusret = pe->locallab.spots.at(pp->locallab.selspot).softradiusret || softradiusret->getEditedState();
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
                        pe->locallab.spots.at(pp->locallab.selspot).clarityml = pe->locallab.spots.at(pp->locallab.selspot).clarityml || clarityml->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).contresid = pe->locallab.spots.at(pp->locallab.selspot).contresid || contresid->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).blurcbdl = pe->locallab.spots.at(pp->locallab.selspot).blurcbdl || blurcbdl->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).softradiuscb = pe->locallab.spots.at(pp->locallab.selspot).softradiuscb || softradiuscb->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).enacbMask = pe->locallab.spots.at(pp->locallab.selspot).enacbMask || !enacbMask->get_inconsistent();

                        pe->locallab.spots.at(pp->locallab.selspot).CCmaskcbcurve = pe->locallab.spots.at(pp->locallab.selspot).CCmaskcbcurve || !CCmaskcbshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).LLmaskcbcurve = pe->locallab.spots.at(pp->locallab.selspot).LLmaskcbcurve || !LLmaskcbshape->isUnChanged();
                        pe->locallab.spots.at(pp->locallab.selspot).HHmaskcbcurve = pe->locallab.spots.at(pp->locallab.selspot).HHmaskcbcurve || !HHmaskcbshape->isUnChanged();

                        pe->locallab.spots.at(pp->locallab.selspot).blendmaskcb = pe->locallab.spots.at(pp->locallab.selspot).blendmaskcb || blendmaskcb->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).radmaskcb = pe->locallab.spots.at(pp->locallab.selspot).radmaskcb || radmaskcb->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).chromaskcb = pe->locallab.spots.at(pp->locallab.selspot).chromaskcb || chromaskcb->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).gammaskcb = pe->locallab.spots.at(pp->locallab.selspot).gammaskcb || gammaskcb->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).slomaskcb = pe->locallab.spots.at(pp->locallab.selspot).slomaskcb || slomaskcb->getEditedState();

                        // Denoise
                        pe->locallab.spots.at(pp->locallab.selspot).expdenoi = pe->locallab.spots.at(pp->locallab.selspot).expdenoi || !expdenoi->get_inconsistent();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselumf = pe->locallab.spots.at(pp->locallab.selspot).noiselumf || noiselumf->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselumf0 = pe->locallab.spots.at(pp->locallab.selspot).noiselumf0 || noiselumf0->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselumf2 = pe->locallab.spots.at(pp->locallab.selspot).noiselumf2 || noiselumf2->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselumc = pe->locallab.spots.at(pp->locallab.selspot).noiselumc || noiselumc->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselumdetail = pe->locallab.spots.at(pp->locallab.selspot).noiselumdetail || noiselumdetail->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noiselequal = pe->locallab.spots.at(pp->locallab.selspot).noiselequal || noiselequal->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noisechrof = pe->locallab.spots.at(pp->locallab.selspot).noisechrof || noisechrof->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noisechroc = pe->locallab.spots.at(pp->locallab.selspot).noisechroc || noisechroc->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).noisechrodetail = pe->locallab.spots.at(pp->locallab.selspot).noisechrodetail || noisechrodetail->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).adjblur = pe->locallab.spots.at(pp->locallab.selspot).adjblur || adjblur->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).bilateral = pe->locallab.spots.at(pp->locallab.selspot).bilateral || bilateral->getEditedState();
                        pe->locallab.spots.at(pp->locallab.selspot).sensiden = pe->locallab.spots.at(pp->locallab.selspot).sensiden || sensiden->getEditedState();
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
                        pedited->locallab.spots.at(pp->locallab.selspot).structexclu = pedited->locallab.spots.at(pp->locallab.selspot).structexclu || se->structexclu;
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
                        pedited->locallab.spots.at(pp->locallab.selspot).balan = pedited->locallab.spots.at(pp->locallab.selspot).balan || se->balan;
                        pedited->locallab.spots.at(pp->locallab.selspot).transitweak = pedited->locallab.spots.at(pp->locallab.selspot).transitweak || se->transitweak;
                        pedited->locallab.spots.at(pp->locallab.selspot).avoid = pedited->locallab.spots.at(pp->locallab.selspot).avoid || se->avoid;
                        // Color & Light
                        pedited->locallab.spots.at(pp->locallab.selspot).expcolor = pedited->locallab.spots.at(pp->locallab.selspot).expcolor || !expcolor->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).curvactiv = pedited->locallab.spots.at(pp->locallab.selspot).curvactiv || !curvactiv->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).lightness = pedited->locallab.spots.at(pp->locallab.selspot).lightness || lightness->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).contrast = pedited->locallab.spots.at(pp->locallab.selspot).contrast || contrast->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).chroma = pedited->locallab.spots.at(pp->locallab.selspot).chroma || chroma->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).strengthgrid = pedited->locallab.spots.at(pp->locallab.selspot).strengthgrid || strengthgrid->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensi = pedited->locallab.spots.at(pp->locallab.selspot).sensi || sensi->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).structcol = pedited->locallab.spots.at(pp->locallab.selspot).structcol || structcol->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod = pedited->locallab.spots.at(pp->locallab.selspot).qualitycurveMethod || qualitycurveMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pedited->locallab.spots.at(pp->locallab.selspot).gridMethod = pedited->locallab.spots.at(pp->locallab.selspot).gridMethod || gridMethod->get_active_text() != M("GENERAL_UNCHANGED");
                        pedited->locallab.spots.at(pp->locallab.selspot).llcurve = pedited->locallab.spots.at(pp->locallab.selspot).llcurve || !llshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).cccurve = pedited->locallab.spots.at(pp->locallab.selspot).cccurve || !ccshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).LHcurve = pedited->locallab.spots.at(pp->locallab.selspot).LHcurve || !LHshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).HHcurve = pedited->locallab.spots.at(pp->locallab.selspot).HHcurve || !HHshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).invers = pedited->locallab.spots.at(pp->locallab.selspot).invers || !invers->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).enaColorMask = pedited->locallab.spots.at(pp->locallab.selspot).enaColorMask || !enaColorMask->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).CCmaskcurve = pedited->locallab.spots.at(pp->locallab.selspot).CCmaskcurve || !CCmaskshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).LLmaskcurve = pedited->locallab.spots.at(pp->locallab.selspot).LLmaskcurve || !LLmaskshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).HHmaskcurve = pedited->locallab.spots.at(pp->locallab.selspot).HHmaskcurve || !HHmaskshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).blurcolde = pedited->locallab.spots.at(pp->locallab.selspot).blurcolde || blurcolde->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).blendmaskcol = pedited->locallab.spots.at(pp->locallab.selspot).blendmaskcol || blendmaskcol->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).radmaskcol = pedited->locallab.spots.at(pp->locallab.selspot).radmaskcol || radmaskcol->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).chromaskcol = pedited->locallab.spots.at(pp->locallab.selspot).chromaskcol || chromaskcol->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).gammaskcol = pedited->locallab.spots.at(pp->locallab.selspot).gammaskcol || gammaskcol->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).slomaskcol = pedited->locallab.spots.at(pp->locallab.selspot).slomaskcol || slomaskcol->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).softradiuscol = pedited->locallab.spots.at(pp->locallab.selspot).softradiuscol || softradiuscol->getEditedState();
                        // Exposure
                        pedited->locallab.spots.at(pp->locallab.selspot).expexpose = pedited->locallab.spots.at(pp->locallab.selspot).expexpose || !expexpose->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).expcomp = pedited->locallab.spots.at(pp->locallab.selspot).expcomp || expcomp->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).hlcompr = pedited->locallab.spots.at(pp->locallab.selspot).hlcompr || hlcompr->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).hlcomprthresh = pedited->locallab.spots.at(pp->locallab.selspot).hlcomprthresh || hlcomprthresh->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).black = pedited->locallab.spots.at(pp->locallab.selspot).black || black->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).shcompr = pedited->locallab.spots.at(pp->locallab.selspot).shcompr || shcompr->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).expchroma = pedited->locallab.spots.at(pp->locallab.selspot).expchroma || expchroma->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).warm = pedited->locallab.spots.at(pp->locallab.selspot).warm || warm->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensiex = pedited->locallab.spots.at(pp->locallab.selspot).sensiex || sensiex->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).structexp = pedited->locallab.spots.at(pp->locallab.selspot).structexp || structexp->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).blurexpde = pedited->locallab.spots.at(pp->locallab.selspot).blurexpde || blurexpde->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).excurve = pedited->locallab.spots.at(pp->locallab.selspot).excurve || !shapeexpos->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).inversex = pedited->locallab.spots.at(pp->locallab.selspot).inversex || !inversex->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).enaExpMask = pedited->locallab.spots.at(pp->locallab.selspot).enaExpMask || !enaExpMask->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve = pedited->locallab.spots.at(pp->locallab.selspot).CCmaskexpcurve || !CCmaskexpshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve = pedited->locallab.spots.at(pp->locallab.selspot).LLmaskexpcurve || !LLmaskexpshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).HHmaskexpcurve = pedited->locallab.spots.at(pp->locallab.selspot).HHmaskexpcurve || !HHmaskexpshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).blendmaskexp = pedited->locallab.spots.at(pp->locallab.selspot).blendmaskexp || blendmaskexp->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).radmaskexp = pedited->locallab.spots.at(pp->locallab.selspot).radmaskexp || radmaskexp->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).chromaskexp = pedited->locallab.spots.at(pp->locallab.selspot).chromaskexp || chromaskexp->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).gammaskexp = pedited->locallab.spots.at(pp->locallab.selspot).gammaskexp || gammaskexp->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).slomaskexp = pedited->locallab.spots.at(pp->locallab.selspot).slomaskexp || slomaskexp->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).softradiusexp = pedited->locallab.spots.at(pp->locallab.selspot).softradiusexp || softradiusexp->getEditedState();
                        // Shadow highlight
                        pedited->locallab.spots.at(pp->locallab.selspot).expshadhigh = pedited->locallab.spots.at(pp->locallab.selspot).expshadhigh || !expshadhigh->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).highlights = pedited->locallab.spots.at(pp->locallab.selspot).highlights || highlights->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).h_tonalwidth = pedited->locallab.spots.at(pp->locallab.selspot).h_tonalwidth || h_tonalwidth->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).shadows = pedited->locallab.spots.at(pp->locallab.selspot).shadows || shadows->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).s_tonalwidth = pedited->locallab.spots.at(pp->locallab.selspot).s_tonalwidth || s_tonalwidth->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sh_radius = pedited->locallab.spots.at(pp->locallab.selspot).sh_radius || sh_radius->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensihs = pedited->locallab.spots.at(pp->locallab.selspot).sensihs || sensihs->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).enaSHMask = pedited->locallab.spots.at(pp->locallab.selspot).enaSHMask || !enaSHMask->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).CCmaskSHcurve = pedited->locallab.spots.at(pp->locallab.selspot).CCmaskSHcurve || !CCmaskSHshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).LLmaskSHcurve = pedited->locallab.spots.at(pp->locallab.selspot).LLmaskSHcurve || !LLmaskSHshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).HHmaskSHcurve = pedited->locallab.spots.at(pp->locallab.selspot).HHmaskSHcurve || !HHmaskSHshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).blendmaskSH = pedited->locallab.spots.at(pp->locallab.selspot).blendmaskSH || blendmaskSH->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).radmaskSH = pedited->locallab.spots.at(pp->locallab.selspot).radmaskSH || radmaskSH->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).blurSHde = pedited->locallab.spots.at(pp->locallab.selspot).blurSHde || blurSHde->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).inverssh = pedited->locallab.spots.at(pp->locallab.selspot).inverssh || !inverssh->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).chromaskSH = pedited->locallab.spots.at(pp->locallab.selspot).chromaskSH || chromaskSH->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).gammaskSH = pedited->locallab.spots.at(pp->locallab.selspot).gammaskSH || gammaskSH->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).slomaskSH = pedited->locallab.spots.at(pp->locallab.selspot).slomaskSH || slomaskSH->getEditedState();
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
                        pedited->locallab.spots.at(pp->locallab.selspot).softradiusret = pedited->locallab.spots.at(pp->locallab.selspot).softradiusret || softradiusret->getEditedState();
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
                        pedited->locallab.spots.at(pp->locallab.selspot).clarityml = pedited->locallab.spots.at(pp->locallab.selspot).clarityml || clarityml->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).contresid = pedited->locallab.spots.at(pp->locallab.selspot).contresid || contresid->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).blurcbdl = pedited->locallab.spots.at(pp->locallab.selspot).blurcbdl || blurcbdl->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).softradiuscb = pedited->locallab.spots.at(pp->locallab.selspot).softradiuscb || softradiuscb->getEditedState();

                        pedited->locallab.spots.at(pp->locallab.selspot).enacbMask = pedited->locallab.spots.at(pp->locallab.selspot).enacbMask || !enacbMask->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).CCmaskcbcurve = pedited->locallab.spots.at(pp->locallab.selspot).CCmaskcbcurve || !CCmaskcbshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).LLmaskcbcurve = pedited->locallab.spots.at(pp->locallab.selspot).LLmaskcbcurve || !LLmaskcbshape->isUnChanged();
                        pedited->locallab.spots.at(pp->locallab.selspot).HHmaskcbcurve = pedited->locallab.spots.at(pp->locallab.selspot).HHmaskcbcurve || !HHmaskcbshape->isUnChanged();

                        pedited->locallab.spots.at(pp->locallab.selspot).blendmaskcb = pedited->locallab.spots.at(pp->locallab.selspot).blendmaskcb || blendmaskcb->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).radmaskcb = pedited->locallab.spots.at(pp->locallab.selspot).radmaskcb || radmaskcb->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).chromaskcb = pedited->locallab.spots.at(pp->locallab.selspot).chromaskcb || chromaskcb->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).gammaskcb = pedited->locallab.spots.at(pp->locallab.selspot).gammaskcb || gammaskcb->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).slomaskcb = pedited->locallab.spots.at(pp->locallab.selspot).slomaskcb || slomaskcb->getEditedState();

                        // Denoise
                        pedited->locallab.spots.at(pp->locallab.selspot).expdenoi = pedited->locallab.spots.at(pp->locallab.selspot).expdenoi || !expdenoi->get_inconsistent();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselumf = pedited->locallab.spots.at(pp->locallab.selspot).noiselumf || noiselumf->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselumf0 = pedited->locallab.spots.at(pp->locallab.selspot).noiselumf0 || noiselumf0->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselumf2 = pedited->locallab.spots.at(pp->locallab.selspot).noiselumf2 || noiselumf2->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselumc = pedited->locallab.spots.at(pp->locallab.selspot).noiselumc || noiselumc->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselumdetail = pedited->locallab.spots.at(pp->locallab.selspot).noiselumdetail || noiselumdetail->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noiselequal = pedited->locallab.spots.at(pp->locallab.selspot).noiselequal || noiselequal->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noisechrof = pedited->locallab.spots.at(pp->locallab.selspot).noisechrof || noisechrof->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noisechroc = pedited->locallab.spots.at(pp->locallab.selspot).noisechroc || noisechroc->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).noisechrodetail = pedited->locallab.spots.at(pp->locallab.selspot).noisechrodetail || noisechrodetail->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).adjblur = pedited->locallab.spots.at(pp->locallab.selspot).adjblur || adjblur->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).bilateral = pedited->locallab.spots.at(pp->locallab.selspot).bilateral || bilateral->getEditedState();
                        pedited->locallab.spots.at(pp->locallab.selspot).sensiden = pedited->locallab.spots.at(pp->locallab.selspot).sensiden || sensiden->getEditedState();
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

        if (ce == HHmaskshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskshape, M("HISTORY_CUSTOMCURVE"));
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

        if (ce == HHmaskexpshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskexpshape, M("HISTORY_CUSTOMCURVE"));
            }
        }
    }

    //Shadows Highlight
    if (getEnabled() && expshadhigh->getEnabled()) {

        if (ce == CCmaskSHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskSHshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == LLmaskSHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskSHshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == HHmaskSHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskSHshape, M("HISTORY_CUSTOMCURVE"));
            }
        }
    }

    //CBDL
    if (getEnabled() && expcbdl->getEnabled()) {

        if (ce == CCmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskcbshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == LLmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskcbshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == HHmaskcbshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskcbshape, M("HISTORY_CUSTOMCURVE"));
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

void Locallab::gridMethodChanged()
{
    // printf("qualitycurveMethodChanged\n");

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvLocallabgridMethod, gridMethod->get_active_text());
        }
    }
}

void Locallab::showmaskcolMethodChanged()
{
    // printf("showmaskcolMethodChanged\n");

    // When one mask state is changed, other masks are deactivated
    disableListener();
        showmaskexpMethod->set_active(0);
        showmaskSHMethod->set_active(0);
        showmaskcbMethod->set_active(0);
    enableListener();

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskcolMethod, "");
    }
}

void Locallab::showmaskexpMethodChanged()
{
    // printf("showmaskexpMethodChanged\n");

    // When one mask state is changed, other masks are deactivated
    disableListener();
        showmaskcolMethod->set_active(0);
        showmaskcbMethod->set_active(0);
        showmaskSHMethod->set_active(0);
    enableListener();

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskexpMethod, "");
    }
}

void Locallab::showmaskSHMethodChanged()
{
    // printf("showmaskSHMethodChanged\n");

    // When one mask state is changed, other masks are deactivated
    disableListener();
        showmaskcolMethod->set_active(0);
        showmaskexpMethod->set_active(0);
        showmaskcbMethod->set_active(0);
    enableListener();

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskSHMethod, "");
    }
}

void Locallab::showmaskcbMethodChanged()
{
    // printf("showmaskSHMethodChanged\n");

    // When one mask state is changed, other masks are deactivated
    disableListener();
       showmaskcolMethod->set_active(0);
       showmaskSHMethod->set_active(0);
       showmaskexpMethod->set_active(0);
    enableListener();

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskcbMethod, "");
    }
}

void Locallab::resetMaskVisibility()
{
    // printf("resetMaskVisibility\n");

    disableListener();
    showmaskcolMethod->set_active(0);
    showmaskexpMethod->set_active(0);
    showmaskSHMethod->set_active(0);
    showmaskcbMethod->set_active(0);
    enableListener();
}

Locallab::llMaskVisibility* Locallab::getMaskVisibility()
{
    llMaskVisibility* maskStruct = new llMaskVisibility();
    maskStruct->colorMask = showmaskcolMethod->get_active_row_number();
    maskStruct->expMask = showmaskexpMethod->get_active_row_number();
    maskStruct->SHMask = showmaskSHMethod->get_active_row_number();
    maskStruct->cbMask = showmaskcbMethod->get_active_row_number();
 //   printf("SHmask=%i \n", maskStruct->SHMask);
 //   printf("cbmask=%i \n", maskStruct->cbMask);
    
    return maskStruct;
}

void Locallab::enaColorMaskChanged()
{
    // printf("enaColorMaskChanged\n");

    if (multiImage) {
        if (enaColorMask->get_inconsistent()) {
            enaColorMask->set_inconsistent(false);
            enaColorMaskConn.block(true);
            enaColorMask->set_active(false);
            enaColorMaskConn.block(false);
        }
    }

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            if (enaColorMask->get_active()) {
                listener->panelChanged(EvLocallabEnaColorMask, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvLocallabEnaColorMask, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::enaExpMaskChanged()
{
    // printf("enaExpMaskChanged\n");

    if (multiImage) {
        if (enaExpMask->get_inconsistent()) {
            enaExpMask->set_inconsistent(false);
            enaExpMaskConn.block(true);
            enaExpMask->set_active(false);
            enaExpMaskConn.block(false);
        }
    }

    if (getEnabled() && expexpose->getEnabled()) {
        if (listener) {
            if (enaExpMask->get_active()) {
                listener->panelChanged(EvLocallabEnaExpMask, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvLocallabEnaExpMask, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::enaSHMaskChanged()
{
    // printf("enaSHMaskChanged\n");

    if (multiImage) {
        if (enaSHMask->get_inconsistent()) {
            enaSHMask->set_inconsistent(false);
            enaSHMaskConn.block(true);
            enaSHMask->set_active(false);
            enaSHMaskConn.block(false);
        }
    }

    if (getEnabled() && expshadhigh->getEnabled()) {
        if (listener) {
            if (enaSHMask->get_active()) {
                listener->panelChanged(EvLocallabEnaSHMask, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvLocallabEnaSHMask, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::enacbMaskChanged()
{
    // printf("enacbMaskChanged\n");

    if (multiImage) {
        if (enacbMask->get_inconsistent()) {
            enacbMask->set_inconsistent(false);
            enacbMaskConn.block(true);
            enacbMask->set_active(false);
            enacbMaskConn.block(false);
        }
    }

    if (getEnabled() && expcbdl->getEnabled()) {
        if (listener) {
            if (enacbMask->get_active()) {
                listener->panelChanged(EvLocallabEnacbMask, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvLocallabEnacbMask, M("GENERAL_DISABLED"));
            }
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
        HCurveEditorG->show();
        curvactiv->hide();
        qualitycurveMethod->show();
        labqualcurv->show();
        expmaskcol->show();
        structcol->show();
        strengthgrid->hide();
        blurcolde->show();
        softradiuscol->show();
        showmaskcolMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
        gridFrame->hide();
    } else if (invers->get_active()) {
        sensi->show();
        llCurveEditorG->show();
        HCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        labqualcurv->hide();
        expmaskcol->hide();
        structcol->hide();
        blurcolde->show();
        gridFrame->hide();
        strengthgrid->hide();
        softradiuscol->hide();

    } else {
        sensi->show();
        llCurveEditorG->show();
        HCurveEditorG->show();
        curvactiv->hide();
        qualitycurveMethod->show();
        labqualcurv->show();
        expmaskcol->show();
        structcol->show();
        blurcolde->show();
        gridFrame->show();
        softradiuscol->show();

        if (batchMode) {
            showmaskcolMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
        }
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


void Locallab::inversexChanged()
{
    // printf("inversChanged\n");

    if (multiImage) {
        if (inversex->get_inconsistent()) {
            inversex->set_inconsistent(false);
            inversexConn.block(true);
            inversex->set_active(false);
            inversexConn.block(false);
        }
    }

    // Update Color & Light GUI according to invers button state (to be compliant with updateSpecificGUIState function)
    if (multiImage && inversex->get_inconsistent()) {
        sensiex->show();
        curveEditorG->show();
        expmaskexp->show();
        structexp->show();
        blurexpde->show();
        softradiusexp->show();
        showmaskexpMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
    } else if (inversex->get_active()) {
        sensiex->show();
        curveEditorG->show();
        expmaskexp->hide();
        structexp->hide();
        blurexpde->show();
        softradiusexp->hide();

    } else {
        sensiex->show();
        curveEditorG->show();
        expmaskexp->show();
        structexp->show();
        blurexpde->show();
        softradiusexp->show();

        if (batchMode) {
            showmaskexpMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
        }
    }

    if (getEnabled() && expexpose->getEnabled()) {
        if (listener) {
            if (inversex->get_active()) {
                listener->panelChanged(Evlocallabinversex, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinversex, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::inversshChanged()
{
    // printf("inversChanged\n");

    if (multiImage) {
        if (inverssh->get_inconsistent()) {
            inverssh->set_inconsistent(false);
            inversshConn.block(true);
            inverssh->set_active(false);
            inversshConn.block(false);
        }
    }

    // Update Color & Light GUI according to invers button state (to be compliant with updateSpecificGUIState function)
    if (multiImage && inverssh->get_inconsistent()) {
        
        sensihs->show();
        blurSHde->show();
        expmasksh->show();

        showmaskSHMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
    } else if (inverssh->get_active()) {
        sensihs->show();
        expmasksh->hide();
        blurSHde->show();

    } else {
        sensihs->show();
        expmasksh->show();
        blurSHde->show();

        if (batchMode) {
            showmaskSHMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
        }
    }

    if (getEnabled() && expshadhigh->getEnabled()) {
        if (listener) {
            if (inverssh->get_active()) {
                listener->panelChanged(Evlocallabinverssh, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinverssh, M("GENERAL_DISABLED"));
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
        sensisha->show();
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
        softradiusret->show();
        
    } else if (inversret->get_active()) {
        sensih->show();
        dehaz->show();
        softradiuscol->hide();
    } else {
        sensih->show();
        dehaz->show();
        softradiusret->show();
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
    // Shadow highlight
    expshadhigh->set_sensitive(cond);
    // Vibrance
    expvibrance->set_sensitive(cond);
    // Soft Light
    expsoft->set_sensitive(cond);
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
    labgrid->setDefault(defSpot->labgridALow / LocallabParams::LABGRIDL_CORR_MAX, defSpot->labgridBLow / LocallabParams::LABGRIDL_CORR_MAX, defSpot->labgridAHigh / LocallabParams::LABGRIDL_CORR_MAX, defSpot->labgridBHigh / LocallabParams::LABGRIDL_CORR_MAX);
    sensi->setDefault((double)defSpot->sensi);
    structcol->setDefault((double)defSpot->structcol);
    blurcolde->setDefault((double)defSpot->blurcolde);
    blendmaskcol->setDefault((double)defSpot->blendmaskcol);
    radmaskcol->setDefault(defSpot->radmaskcol);
    chromaskcol->setDefault(defSpot->chromaskcol);
    gammaskcol->setDefault(defSpot->gammaskcol);
    slomaskcol->setDefault(defSpot->slomaskcol);
    softradiuscol->setDefault(defSpot->softradiuscol);
    // Exposure
    expcomp->setDefault(defSpot->expcomp);
    hlcompr->setDefault((double)defSpot->hlcompr);
    hlcomprthresh->setDefault((double)defSpot->hlcomprthresh);
    black->setDefault((double)defSpot->black);
    shcompr->setDefault((double)defSpot->shcompr);
    expchroma->setDefault((double)defSpot->expchroma);
    warm->setDefault((double)defSpot->warm);
    sensiex->setDefault((double)defSpot->sensiex);
    structexp->setDefault((double)defSpot->structexp);
    blurexpde->setDefault((double)defSpot->blurexpde);
    blendmaskexp->setDefault((double)defSpot->blendmaskexp);
    radmaskexp->setDefault(defSpot->radmaskexp);
    chromaskexp->setDefault(defSpot->chromaskexp);
    gammaskexp->setDefault(defSpot->gammaskexp);
    slomaskexp->setDefault(defSpot->slomaskexp);
    softradiusexp->setDefault(defSpot->softradiusexp);
    // Shadow highlight
    highlights->setDefault((double)defSpot->highlights);
    h_tonalwidth->setDefault((double)defSpot->h_tonalwidth);
    shadows->setDefault((double)defSpot->shadows);
    s_tonalwidth->setDefault((double)defSpot->s_tonalwidth);
    sh_radius->setDefault((double)defSpot->sh_radius);
    sensihs->setDefault((double)defSpot->sensihs);
    blendmaskSH->setDefault((double)defSpot->blendmaskSH);
    radmaskSH->setDefault(defSpot->radmaskSH);
    blurSHde->setDefault((double)defSpot->blurSHde);
    chromaskSH->setDefault(defSpot->chromaskSH);
    gammaskSH->setDefault(defSpot->gammaskSH);
    slomaskSH->setDefault(defSpot->slomaskSH);
    // Vibrance
    saturated->setDefault((double)defSpot->saturated);
    pastels->setDefault((double)defSpot->pastels);
    psThreshold->setDefault<int>(defSpot->psthreshold);
    sensiv->setDefault((double)defSpot->sensiv);
    // Soft Light
    streng->setDefault((double)defSpot->streng);
    sensisf->setDefault((double)defSpot->sensisf);
    // Blur & Noise
    radius->setDefault(defSpot->radius);
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
    softradiusret->setDefault(defSpot->softradiusret);
    // Sharpening
    sharcontrast->setDefault((double)defSpot->sharcontrast);
    sharradius->setDefault(defSpot->sharradius);
    sharamount->setDefault((double)defSpot->sharamount);
    shardamping->setDefault((double)defSpot->shardamping);
    shariter->setDefault((double)defSpot->shariter);
    sharblur->setDefault(defSpot->sharblur);
    sensisha->setDefault((double)defSpot->sensisha);
    // Local Contrast
    lcradius->setDefault((double)defSpot->lcradius);
    lcamount->setDefault(defSpot->lcamount);
    lcdarkness->setDefault(defSpot->lcdarkness);
    lclightness->setDefault(defSpot->lclightness);
    sensilc->setDefault((double)defSpot->sensilc);
    // Contrast by detail levels
    for (int i = 0; i < 5; i++) {
        multiplier[i]->setDefault(defSpot->mult[i]);
    }

    chromacbdl->setDefault((double)defSpot->chromacbdl);
    threshold->setDefault(defSpot->threshold);
    sensicb->setDefault((double)defSpot->sensicb);
    clarityml->setDefault((double)defSpot->clarityml);
    contresid->setDefault((double)defSpot->contresid);
    blurcbdl->setDefault((double)defSpot->blurcbdl);
    softradiuscb->setDefault(defSpot->softradiuscb);
    blendmaskcb->setDefault((double)defSpot->blendmaskcb);
    radmaskcb->setDefault(defSpot->radmaskcb);
    chromaskcb->setDefault(defSpot->chromaskcb);
    gammaskcb->setDefault(defSpot->gammaskcb);
    slomaskcb->setDefault(defSpot->slomaskcb);
    
    // Denoise
    noiselumf->setDefault((double)defSpot->noiselumf);
    noiselumf0->setDefault((double)defSpot->noiselumf0);
    noiselumf2->setDefault((double)defSpot->noiselumf2);
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
        labgrid->setEdited(Edited);
        sensi->setDefaultEditedState(Irrelevant);
        structcol->setDefaultEditedState(Irrelevant);
        strengthgrid->setDefault((double)defSpot->strengthgrid);
        blurcolde->setDefaultEditedState(Irrelevant);
        blendmaskcol->setDefaultEditedState(Irrelevant);
        radmaskcol->setDefaultEditedState(Irrelevant);
        chromaskcol->setDefaultEditedState(Irrelevant);
        gammaskcol->setDefaultEditedState(Irrelevant);
        slomaskcol->setDefaultEditedState(Irrelevant);
        softradiuscol->setDefaultEditedState(Irrelevant);
        // Exposure
        expcomp->setDefaultEditedState(Irrelevant);
        hlcompr->setDefaultEditedState(Irrelevant);
        hlcomprthresh->setDefaultEditedState(Irrelevant);
        black->setDefaultEditedState(Irrelevant);
        shcompr->setDefaultEditedState(Irrelevant);
        expchroma->setDefaultEditedState(Irrelevant);
        warm->setDefaultEditedState(Irrelevant);
        sensiex->setDefaultEditedState(Irrelevant);
        structexp->setDefaultEditedState(Irrelevant);
        blurexpde->setDefaultEditedState(Irrelevant);
        blendmaskexp->setDefaultEditedState(Irrelevant);
        radmaskexp->setDefaultEditedState(Irrelevant);
        chromaskexp->setDefaultEditedState(Irrelevant);
        gammaskexp->setDefaultEditedState(Irrelevant);
        slomaskexp->setDefaultEditedState(Irrelevant);
        softradiusexp->setDefaultEditedState(Irrelevant);
        // Shadow highlight
        highlights->setDefaultEditedState(Irrelevant);
        h_tonalwidth->setDefaultEditedState(Irrelevant);
        shadows->setDefaultEditedState(Irrelevant);
        s_tonalwidth->setDefaultEditedState(Irrelevant);
        sh_radius->setDefaultEditedState(Irrelevant);
        sensihs->setDefaultEditedState(Irrelevant);
        blendmaskSH->setDefaultEditedState(Irrelevant);
        radmaskSH->setDefaultEditedState(Irrelevant);
        blurSHde->setDefaultEditedState(Irrelevant);
        chromaskSH->setDefaultEditedState(Irrelevant);
        gammaskSH->setDefaultEditedState(Irrelevant);
        slomaskSH->setDefaultEditedState(Irrelevant);
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
        softradiusret->setDefaultEditedState(Irrelevant);
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
        clarityml->setDefaultEditedState(Irrelevant);
        contresid->setDefaultEditedState(Irrelevant);
        blurcbdl->setDefaultEditedState(Irrelevant);
        softradiuscb->setDefaultEditedState(Irrelevant);
        blendmaskcb->setDefaultEditedState(Irrelevant);
        radmaskcb->setDefaultEditedState(Irrelevant);
        chromaskcb->setDefaultEditedState(Irrelevant);
        gammaskcb->setDefaultEditedState(Irrelevant);
        slomaskcb->setDefaultEditedState(Irrelevant);

        // Denoise
        noiselumf->setDefaultEditedState(Irrelevant);
        noiselumf0->setDefaultEditedState(Irrelevant);
        noiselumf2->setDefaultEditedState(Irrelevant);
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
        labgrid->setEdited((defSpotState->labgridALow || defSpotState->labgridBLow || defSpotState->labgridAHigh || defSpotState->labgridBHigh) ? Edited : UnEdited);
        sensi->setDefaultEditedState(defSpotState->sensi ? Edited : UnEdited);
        structcol->setDefaultEditedState(defSpotState->structcol ? Edited : UnEdited);
        strengthgrid->setDefaultEditedState(defSpotState->strengthgrid ? Edited : UnEdited);
        blurcolde->setDefaultEditedState(defSpotState->blurcolde ? Edited : UnEdited);
        blendmaskcol->setDefaultEditedState(defSpotState->blendmaskcol ? Edited : UnEdited);
        radmaskcol->setDefaultEditedState(defSpotState->radmaskcol ? Edited : UnEdited);
        chromaskcol->setDefaultEditedState(defSpotState->chromaskcol ? Edited : UnEdited);
        gammaskcol->setDefaultEditedState(defSpotState->gammaskcol ? Edited : UnEdited);
        slomaskcol->setDefaultEditedState(defSpotState->slomaskcol ? Edited : UnEdited);
        softradiuscol->setDefaultEditedState(defSpotState->softradiuscol ? Edited : UnEdited);
        // Exposure
        expcomp->setDefaultEditedState(defSpotState->expcomp ? Edited : UnEdited);
        hlcompr->setDefaultEditedState(defSpotState->hlcompr ? Edited : UnEdited);
        hlcomprthresh->setDefaultEditedState(defSpotState->hlcomprthresh ? Edited : UnEdited);
        black->setDefaultEditedState(defSpotState->black ? Edited : UnEdited);
        shcompr->setDefaultEditedState(defSpotState->shcompr ? Edited : UnEdited);
        expchroma->setDefaultEditedState(defSpotState->expchroma ? Edited : UnEdited);
        warm->setDefaultEditedState(defSpotState->warm ? Edited : UnEdited);
        sensiex->setDefaultEditedState(defSpotState->sensiex ? Edited : UnEdited);
        structexp->setDefaultEditedState(defSpotState->structexp ? Edited : UnEdited);
        blurexpde->setDefaultEditedState(defSpotState->blurexpde ? Edited : UnEdited);
        blendmaskexp->setDefaultEditedState(defSpotState->blendmaskexp ? Edited : UnEdited);
        radmaskexp->setDefaultEditedState(defSpotState->radmaskexp ? Edited : UnEdited);
        chromaskexp->setDefaultEditedState(defSpotState->chromaskexp ? Edited : UnEdited);
        gammaskexp->setDefaultEditedState(defSpotState->gammaskexp ? Edited : UnEdited);
        slomaskexp->setDefaultEditedState(defSpotState->slomaskexp ? Edited : UnEdited);
        softradiusexp->setDefaultEditedState(defSpotState->softradiusexp ? Edited : UnEdited);
        // Shadow highlight
        highlights->setDefaultEditedState(defSpotState->highlights ? Edited : UnEdited);
        h_tonalwidth->setDefaultEditedState(defSpotState->h_tonalwidth ? Edited : UnEdited);
        shadows->setDefaultEditedState(defSpotState->shadows ? Edited : UnEdited);
        s_tonalwidth->setDefaultEditedState(defSpotState->s_tonalwidth ? Edited : UnEdited);
        sh_radius->setDefaultEditedState(defSpotState->sh_radius ? Edited : UnEdited);
        sensihs->setDefaultEditedState(defSpotState->sensihs ? Edited : UnEdited);
        blendmaskSH->setDefaultEditedState(defSpotState->blendmaskSH ? Edited : UnEdited);
        radmaskSH->setDefaultEditedState(defSpotState->radmaskSH ? Edited : UnEdited);
        blurSHde->setDefaultEditedState(defSpotState->blurSHde ? Edited : UnEdited);
        chromaskSH->setDefaultEditedState(defSpotState->chromaskSH ? Edited : UnEdited);
        gammaskSH->setDefaultEditedState(defSpotState->gammaskSH ? Edited : UnEdited);
        slomaskSH->setDefaultEditedState(defSpotState->slomaskSH ? Edited : UnEdited);
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
        softradiusret->setDefaultEditedState(defSpotState->softradiusret ? Edited : UnEdited);
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
        clarityml->setDefaultEditedState(defSpotState->clarityml ? Edited : UnEdited);
        contresid->setDefaultEditedState(defSpotState->contresid ? Edited : UnEdited);
        blurcbdl->setDefaultEditedState(defSpotState->blurcbdl ? Edited : UnEdited);
        softradiuscb->setDefaultEditedState(defSpotState->softradiuscb ? Edited : UnEdited);

        blendmaskcb->setDefaultEditedState(defSpotState->blendmaskcb ? Edited : UnEdited);
        radmaskcb->setDefaultEditedState(defSpotState->radmaskcb ? Edited : UnEdited);
        chromaskcb->setDefaultEditedState(defSpotState->chromaskcb ? Edited : UnEdited);
        gammaskcb->setDefaultEditedState(defSpotState->gammaskcb ? Edited : UnEdited);
        slomaskcb->setDefaultEditedState(defSpotState->slomaskcb ? Edited : UnEdited);
        
        // Denoise
        noiselumf->setDefaultEditedState(defSpotState->noiselumf ? Edited : UnEdited);
        noiselumf0->setDefaultEditedState(defSpotState->noiselumf0 ? Edited : UnEdited);
        noiselumf2->setDefaultEditedState(defSpotState->noiselumf2 ? Edited : UnEdited);
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

        if (a == strengthgrid) {
            if (listener) {
                listener->panelChanged(EvLocallabLabstrengthgrid, strengthgrid->getTextValue());
            }
        }

        if (a == sensi) {
            if (listener) {
                listener->panelChanged(Evlocallabsensi, sensi->getTextValue());
            }
        }

        if (a == blurcolde) {
            if (listener) {
                listener->panelChanged(Evlocallabblurcolde, blurcolde->getTextValue());
            }
        }

        if (a == structcol) {
            if (listener) {
                listener->panelChanged(Evlocallabstructcol, structcol->getTextValue());
            }
        }

        if (a == blendmaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskcol, blendmaskcol->getTextValue());
            }
        }

        if (a == radmaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskcol, radmaskcol->getTextValue());
            }
        }

        if (a == chromaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskcol, chromaskcol->getTextValue());
            }
        }

        if (a == gammaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskcol, gammaskcol->getTextValue());
            }
        }

        if (a == slomaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskcol, slomaskcol->getTextValue());
            }
        }

        if (a == softradiuscol) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiuscol, softradiuscol->getTextValue());
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

        if (a == expchroma) {
            if (listener) {
                listener->panelChanged(Evlocallabexpchroma, expchroma->getTextValue());
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

        if (a == structexp) {
            if (listener) {
                listener->panelChanged(Evlocallabstructexp, structexp->getTextValue());
            }
        }

        if (a == blurexpde) {
            if (listener) {
                listener->panelChanged(Evlocallabblurexpde, blurexpde->getTextValue());
            }
        }

        if (a == blendmaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskexp, blendmaskexp->getTextValue());
            }
        }

        if (a == radmaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskexp, radmaskexp->getTextValue());
            }
        }

        if (a == chromaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskexp, chromaskexp->getTextValue());
            }
        }

        if (a == gammaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskexp, gammaskexp->getTextValue());
            }
        }

        if (a == slomaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskexp, slomaskexp->getTextValue());
            }
        }

        if (a == softradiusexp) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusexp, softradiusexp->getTextValue());
            }
        }

    }

    if (getEnabled() && expshadhigh->getEnabled()) {

        if (a == highlights) {
            if (listener) {
                listener->panelChanged(Evlocallabhighlights, highlights->getTextValue());
            }
        }

        if (a == h_tonalwidth) {
            if (listener) {
                listener->panelChanged(Evlocallabh_tonalwidth, h_tonalwidth->getTextValue());
            }
        }

        if (a == shadows) {
            if (listener) {
                listener->panelChanged(Evlocallabshadows, shadows->getTextValue());
            }
        }

        if (a == s_tonalwidth) {
            if (listener) {
                listener->panelChanged(Evlocallabs_tonalwidth, s_tonalwidth->getTextValue());
            }
        }

        if (a == sh_radius) {
            if (listener) {
                listener->panelChanged(Evlocallabsh_radius, sh_radius->getTextValue());
            }
        }

        if (a == sensihs) {
            if (listener) {
                listener->panelChanged(Evlocallabsensihs, sensihs->getTextValue());
            }
        }

        if (a == blendmaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabblendmaskSH, blendmaskSH->getTextValue());
            }
        }

        if (a == radmaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabradmaskSH, radmaskSH->getTextValue());
            }
        }

        if (a == blurSHde) {
            if (listener) {
                listener->panelChanged(EvlocallabblurSHde, blurSHde->getTextValue());
            }
        }

        if (a == chromaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabchromaskSH, chromaskSH->getTextValue());
            }
        }

        if (a == gammaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabgammaskSH, gammaskSH->getTextValue());
            }
        }

        if (a == slomaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabslomaskSH, slomaskSH->getTextValue());
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

        if (a == softradiusret) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusret, softradiusret->getTextValue());
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

        if (a == clarityml) {
            //contresid->setValue(0.);

            if (listener) {
                listener->panelChanged(EvLocallabclarityml, clarityml->getTextValue());
            }
        }

        if (a == contresid) {

            if (listener) {
                listener->panelChanged(EvLocallabcontresid, contresid->getTextValue());
            }
        }

        if (a == blurcbdl) {

            if (listener) {
                listener->panelChanged(EvLocallabblurcbdl, blurcbdl->getTextValue());
            }
        }

        if (a == softradiuscb) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiuscb, softradiuscb->getTextValue());
            }
        }

        if (a == blendmaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskcb, blendmaskcb->getTextValue());
            }
        }

        if (a == radmaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskcb, radmaskcb->getTextValue());
            }
        }

        if (a == chromaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskcb, chromaskcb->getTextValue());
            }
        }

        if (a == gammaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskcb, gammaskcb->getTextValue());
            }
        }

        if (a == slomaskcb) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskcb, slomaskcb->getTextValue());
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

        if (a == noiselumf0) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf0, noiselumf0->getTextValue());
            }
        }

        if (a == noiselumf2) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf2, noiselumf2->getTextValue());
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
    structcol->showEditedCB();
    strengthgrid->showEditedCB();
    blurcolde->showEditedCB();
    blendmaskcol->showEditedCB();
    radmaskcol->showEditedCB();
    chromaskcol->showEditedCB();
    gammaskcol->showEditedCB();
    slomaskcol->showEditedCB();
    softradiuscol->showEditedCB();
    // Exposure
    expcomp->showEditedCB();
    hlcompr->showEditedCB();
    hlcomprthresh->showEditedCB();
    black->showEditedCB();
    shcompr->showEditedCB();
    expchroma->showEditedCB();
    warm->showEditedCB();
    sensiex->showEditedCB();
    structexp->showEditedCB();
    blurexpde->showEditedCB();
    blendmaskexp->showEditedCB();
    radmaskexp->showEditedCB();
    chromaskexp->showEditedCB();
    gammaskexp->showEditedCB();
    slomaskexp->showEditedCB();
    softradiusexp->showEditedCB();
    //Shadow Highlight
    highlights->showEditedCB();
    h_tonalwidth->showEditedCB();
    shadows->showEditedCB();
    s_tonalwidth->showEditedCB();
    sh_radius->showEditedCB();
    sensihs->showEditedCB();
    blendmaskSH->showEditedCB();
    radmaskSH->showEditedCB();
    blurSHde->showEditedCB();
    chromaskSH->showEditedCB();
    gammaskSH->showEditedCB();
    slomaskSH->showEditedCB();
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
    softradiusret->showEditedCB();
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
    clarityml->showEditedCB();
    contresid->showEditedCB();
    blurcbdl->showEditedCB();
    softradiuscb->showEditedCB();
    blendmaskcb->showEditedCB();
    radmaskcb->showEditedCB();
    chromaskcb->showEditedCB();
    gammaskcb->showEditedCB();
    slomaskcb->showEditedCB();

    // Denoise
    noiselumf->showEditedCB();
    noiselumc->showEditedCB();
    noiselumdetail->showEditedCB();
    noiselequal->showEditedCB();
    noiselumf0->showEditedCB();
    noiselumf2->showEditedCB();
    noisechroc->showEditedCB();
    noisechrodetail->showEditedCB();
    adjblur->showEditedCB();
    bilateral->showEditedCB();
    sensiden->showEditedCB();

    // Set batch mode for comboBoxText
    // Color & Light
    qualitycurveMethod->append(M("GENERAL_UNCHANGED"));
    gridMethod->append(M("GENERAL_UNCHANGED"));
    // Blur & Noise
    blurMethod->append(M("GENERAL_UNCHANGED"));
    // Retinex
    retinexMethod->append(M("GENERAL_UNCHANGED"));

    // In batch mode, being able to change mask visibility is useless
    showmaskcolMethod->hide();
    showmaskexpMethod->hide();
    showmaskSHMethod->hide();
    showmaskcbMethod->hide();
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

    float R = 0.f;
    float G = 0.f;
    float B = 0.f;

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
        float x = valX - 1.f / 6.f;

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
    labgrid->setListener(tpl);
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
    gridMethodConn.block(false);
    inversConn.block(false);
    showmaskcolMethodConn.block(false);
    enaColorMaskConn.block(false);
    // Exposure
    enableexposeConn.block(false);
    inversexConn.block(false);
    showmaskexpMethodConn.block(false);
    enaExpMaskConn.block(false);
    // Shadow highlight
    enableshadhighConn.block(false);
    showmaskSHMethodConn.block(false);
    enaSHMaskConn.block(false);
    inversshConn.block(false);
    // Vibrance
    enablevibranceConn.block(false);
    pskinsconn.block(false);
    ashiftconn.block(false);
    pastsattogconn.block(false);
    // Soft Light
    enablesoftConn.block(false);
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
    enacbMaskConn.block(false);
    showmaskcbMethodConn.block(false);
    // Denoise
    enabledenoiConn.block(false);
}

void Locallab::disableListener()
{
    // printf("disableListener\n");

    FoldableToolPanel::disableListener();
    // Color & Light
    enablecolorConn.block(true);
    curvactivConn.block(true);
    qualitycurveMethodConn.block(true);
    gridMethodConn.block(true);
    inversConn.block(true);
    showmaskcolMethodConn.block(true);
    enaColorMaskConn.block(true);
    // Exposure
    enableexposeConn.block(true);
    inversexConn.block(true);
    showmaskexpMethodConn.block(true);
    enaExpMaskConn.block(true);
    // Shadow highlight
    enableshadhighConn.block(true);
    showmaskSHMethodConn.block(true);
    enaSHMaskConn.block(true);
    inversshConn.block(true);
    // Vibrance
    enablevibranceConn.block(true);
    pskinsconn.block(true);
    ashiftconn.block(true);
    pastsattogconn.block(true);
    // Soft Light
    enablesoftConn.block(true);
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
    enacbMaskConn.block(true);
    showmaskcbMethodConn.block(true);
    // Denoise
    enabledenoiConn.block(true);
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
        labgrid->setParams(pp->locallab.spots.at(index).labgridALow / LocallabParams::LABGRIDL_CORR_MAX, pp->locallab.spots.at(index).labgridBLow / LocallabParams::LABGRIDL_CORR_MAX, pp->locallab.spots.at(index).labgridAHigh / LocallabParams::LABGRIDL_CORR_MAX, pp->locallab.spots.at(index).labgridBHigh / LocallabParams::LABGRIDL_CORR_MAX, false);
        strengthgrid->setValue(pp->locallab.spots.at(index).strengthgrid);
        sensi->setValue(pp->locallab.spots.at(index).sensi);
        structcol->setValue(pp->locallab.spots.at(index).structcol);

        if (pp->locallab.spots.at(index).qualitycurveMethod == "none") {
            qualitycurveMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).qualitycurveMethod == "std") {
            qualitycurveMethod->set_active(1);
        }

        if (pp->locallab.spots.at(index).gridMethod == "one") {
            gridMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).gridMethod == "two") {
            gridMethod->set_active(1);
        }

        llshape->setCurve(pp->locallab.spots.at(index).llcurve);
        ccshape->setCurve(pp->locallab.spots.at(index).cccurve);
        LHshape->setCurve(pp->locallab.spots.at(index).LHcurve);
        HHshape->setCurve(pp->locallab.spots.at(index).HHcurve);
        invers->set_active(pp->locallab.spots.at(index).invers);
        enaColorMask->set_active(pp->locallab.spots.at(index).enaColorMask);
        CCmaskshape->setCurve(pp->locallab.spots.at(index).CCmaskcurve);
        LLmaskshape->setCurve(pp->locallab.spots.at(index).LLmaskcurve);
        HHmaskshape->setCurve(pp->locallab.spots.at(index).HHmaskcurve);
        blurcolde->setValue(pp->locallab.spots.at(index).blurcolde);
        blendmaskcol->setValue(pp->locallab.spots.at(index).blendmaskcol);
        radmaskcol->setValue(pp->locallab.spots.at(index).radmaskcol);
        chromaskcol->setValue(pp->locallab.spots.at(index).chromaskcol);
        gammaskcol->setValue(pp->locallab.spots.at(index).gammaskcol);
        slomaskcol->setValue(pp->locallab.spots.at(index).slomaskcol);
        softradiuscol->setValue(pp->locallab.spots.at(index).softradiuscol);

        // Exposure
        expexpose->setEnabled(pp->locallab.spots.at(index).expexpose);
        expcomp->setValue(pp->locallab.spots.at(index).expcomp);
        hlcompr->setValue(pp->locallab.spots.at(index).hlcompr);
        hlcomprthresh->setValue(pp->locallab.spots.at(index).hlcomprthresh);
        black->setValue(pp->locallab.spots.at(index).black);
        shcompr->setValue(pp->locallab.spots.at(index).shcompr);
        expchroma->setValue(pp->locallab.spots.at(index).expchroma);
        warm->setValue(pp->locallab.spots.at(index).warm);
        sensiex->setValue(pp->locallab.spots.at(index).sensiex);
        structexp->setValue(pp->locallab.spots.at(index).structexp);
        blurexpde->setValue(pp->locallab.spots.at(index).blurexpde);
        shapeexpos->setCurve(pp->locallab.spots.at(index).excurve);
        inversex->set_active(pp->locallab.spots.at(index).inversex);
        enaExpMask->set_active(pp->locallab.spots.at(index).enaExpMask);
        CCmaskexpshape->setCurve(pp->locallab.spots.at(index).CCmaskexpcurve);
        LLmaskexpshape->setCurve(pp->locallab.spots.at(index).LLmaskexpcurve);
        HHmaskexpshape->setCurve(pp->locallab.spots.at(index).HHmaskexpcurve);
        blendmaskexp->setValue(pp->locallab.spots.at(index).blendmaskexp);
        radmaskexp->setValue(pp->locallab.spots.at(index).radmaskexp);
        chromaskexp->setValue(pp->locallab.spots.at(index).chromaskexp);
        gammaskexp->setValue(pp->locallab.spots.at(index).gammaskexp);
        slomaskexp->setValue(pp->locallab.spots.at(index).slomaskexp);
        softradiusexp->setValue(pp->locallab.spots.at(index).softradiusexp);

        // Shadow highlight
        expshadhigh->setEnabled(pp->locallab.spots.at(index).expshadhigh);
        highlights->setValue(pp->locallab.spots.at(index).highlights);
        h_tonalwidth->setValue(pp->locallab.spots.at(index).h_tonalwidth);
        shadows->setValue(pp->locallab.spots.at(index).shadows);
        s_tonalwidth->setValue(pp->locallab.spots.at(index).s_tonalwidth);
        sh_radius->setValue(pp->locallab.spots.at(index).sh_radius);
        sensihs->setValue(pp->locallab.spots.at(index).sensihs);
        enaSHMask->set_active(pp->locallab.spots.at(index).enaSHMask);
        CCmaskSHshape->setCurve(pp->locallab.spots.at(index).CCmaskSHcurve);
        LLmaskSHshape->setCurve(pp->locallab.spots.at(index).LLmaskSHcurve);
        HHmaskSHshape->setCurve(pp->locallab.spots.at(index).HHmaskSHcurve);
        blendmaskSH->setValue(pp->locallab.spots.at(index).blendmaskSH);
        radmaskSH->setValue(pp->locallab.spots.at(index).radmaskSH);
        blurSHde->setValue(pp->locallab.spots.at(index).blurSHde);
        inverssh->set_active(pp->locallab.spots.at(index).inverssh);
        chromaskSH->setValue(pp->locallab.spots.at(index).chromaskSH);
        gammaskSH->setValue(pp->locallab.spots.at(index).gammaskSH);
        slomaskSH->setValue(pp->locallab.spots.at(index).slomaskSH);

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

        // Blur & Noise
        expblur->setEnabled(pp->locallab.spots.at(index).expblur);
        radius->setValue(pp->locallab.spots.at(index).radius);
        strength->setValue(pp->locallab.spots.at(index).strength);
        sensibn->setValue(pp->locallab.spots.at(index).sensibn);

        if (pp->locallab.spots.at(index).blurMethod == "norm") {
            blurMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).blurMethod == "inv") {
            blurMethod->set_active(1);
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
        softradiusret->setValue(pp->locallab.spots.at(index).softradiusret);

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
        clarityml->setValue(pp->locallab.spots.at(index).clarityml);
        contresid->setValue(pp->locallab.spots.at(index).contresid);
        blurcbdl->setValue(pp->locallab.spots.at(index).blurcbdl);
        softradiuscb->setValue(pp->locallab.spots.at(index).softradiuscb);
        blendmaskcb->setValue(pp->locallab.spots.at(index).blendmaskcb);
        radmaskcb->setValue(pp->locallab.spots.at(index).radmaskcb);
        chromaskcb->setValue(pp->locallab.spots.at(index).chromaskcb);
        gammaskcb->setValue(pp->locallab.spots.at(index).gammaskcb);
        slomaskcb->setValue(pp->locallab.spots.at(index).slomaskcb);
        enacbMask->set_active(pp->locallab.spots.at(index).enacbMask);
        CCmaskcbshape->setCurve(pp->locallab.spots.at(index).CCmaskcbcurve);
        LLmaskcbshape->setCurve(pp->locallab.spots.at(index).LLmaskcbcurve);
        HHmaskcbshape->setCurve(pp->locallab.spots.at(index).HHmaskcbcurve);

        // Denoise
        expdenoi->setEnabled(pp->locallab.spots.at(index).expdenoi);
        noiselumf->setValue(pp->locallab.spots.at(index).noiselumf);
        noiselumf0->setValue(pp->locallab.spots.at(index).noiselumf0);
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
                se->structexclu = spotState->structexclu;
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
                se->balan = spotState->balan;
                se->transitweak = spotState->transitweak;
                se->avoid = spotState->avoid;
                expsettings->setEditedStates(se);

                // Color & Light
                expcolor->set_inconsistent(!spotState->expcolor);
                curvactiv->set_inconsistent(multiImage && !spotState->curvactiv);
                lightness->setEditedState(spotState->lightness ? Edited : UnEdited);
                contrast->setEditedState(spotState->contrast ? Edited : UnEdited);
                chroma->setEditedState(spotState->chroma ? Edited : UnEdited);
                sensi->setEditedState(spotState->sensi ? Edited : UnEdited);
                structcol->setEditedState(spotState->structcol ? Edited : UnEdited);
                labgrid->setEdited(spotState->labgridALow || spotState->labgridBLow || spotState->labgridAHigh || spotState->labgridBHigh);
                strengthgrid->setEditedState(spotState->strengthgrid ? Edited : UnEdited);

                if (!spotState->qualitycurveMethod) {
                    qualitycurveMethod->set_active_text(M("GENERAL_UNCHANGED"));
                }

                if (!spotState->gridMethod) {
                    gridMethod->set_active_text(M("GENERAL_UNCHANGED"));
                }

                llshape->setUnChanged(!spotState->llcurve);
                ccshape->setUnChanged(!spotState->cccurve);
                LHshape->setUnChanged(!spotState->LHcurve);
                HHshape->setUnChanged(!spotState->HHcurve);
                invers->set_inconsistent(multiImage && !spotState->invers);
                enaColorMask->set_inconsistent(multiImage && !spotState->enaColorMask);
                CCmaskshape->setUnChanged(!spotState->CCmaskcurve);
                LLmaskshape->setUnChanged(!spotState->LLmaskcurve);
                HHmaskshape->setUnChanged(!spotState->HHmaskcurve);
                blurcolde->setEditedState(spotState->blurcolde ? Edited : UnEdited);
                blendmaskcol->setEditedState(spotState->blendmaskcol ? Edited : UnEdited);
                radmaskcol->setEditedState(spotState->radmaskcol ? Edited : UnEdited);
                chromaskcol->setEditedState(spotState->chromaskcol ? Edited : UnEdited);
                gammaskcol->setEditedState(spotState->gammaskcol ? Edited : UnEdited);
                slomaskcol->setEditedState(spotState->slomaskcol ? Edited : UnEdited);
                softradiuscol->setEditedState(spotState->softradiuscol ? Edited : UnEdited);

                // Exposure
                expexpose->set_inconsistent(!spotState->expexpose);
                expcomp->setEditedState(spotState->expcomp ? Edited : UnEdited);
                hlcompr->setEditedState(spotState->hlcompr ? Edited : UnEdited);
                hlcomprthresh->setEditedState(spotState->hlcomprthresh ? Edited : UnEdited);
                black->setEditedState(spotState->black ? Edited : UnEdited);
                warm->setEditedState(spotState->warm ? Edited : UnEdited);
                shcompr->setEditedState(spotState->shcompr ? Edited : UnEdited);
                expchroma->setEditedState(spotState->expchroma ? Edited : UnEdited);
                sensiex->setEditedState(spotState->sensiex ? Edited : UnEdited);
                structexp->setEditedState(spotState->structexp ? Edited : UnEdited);
                blurexpde->setEditedState(spotState->blurexpde ? Edited : UnEdited);
                shapeexpos->setUnChanged(!spotState->excurve);
                inversex->set_inconsistent(multiImage && !spotState->inversex);
                enaExpMask->set_inconsistent(multiImage && !spotState->enaExpMask);
                CCmaskexpshape->setUnChanged(!spotState->CCmaskexpcurve);
                LLmaskexpshape->setUnChanged(!spotState->LLmaskexpcurve);
                HHmaskexpshape->setUnChanged(!spotState->HHmaskexpcurve);
                blendmaskexp->setEditedState(spotState->blendmaskexp ? Edited : UnEdited);
                radmaskexp->setEditedState(spotState->radmaskexp ? Edited : UnEdited);
                chromaskexp->setEditedState(spotState->chromaskexp ? Edited : UnEdited);
                gammaskexp->setEditedState(spotState->gammaskexp ? Edited : UnEdited);
                slomaskexp->setEditedState(spotState->slomaskexp ? Edited : UnEdited);
                softradiusexp->setEditedState(spotState->softradiusexp ? Edited : UnEdited);

                // Shadow highlight
                expshadhigh->set_inconsistent(!spotState->expshadhigh);
                highlights->setEditedState(spotState->highlights ? Edited : UnEdited);
                h_tonalwidth->setEditedState(spotState->h_tonalwidth ? Edited : UnEdited);
                shadows->setEditedState(spotState->shadows ? Edited : UnEdited);
                s_tonalwidth->setEditedState(spotState->s_tonalwidth ? Edited : UnEdited);
                sh_radius->setEditedState(spotState->sh_radius ? Edited : UnEdited);
                sensihs->setEditedState(spotState->sensihs ? Edited : UnEdited);
                enaSHMask->set_inconsistent(multiImage && !spotState->enaSHMask);
                CCmaskSHshape->setUnChanged(!spotState->CCmaskSHcurve);
                LLmaskSHshape->setUnChanged(!spotState->LLmaskSHcurve);
                HHmaskSHshape->setUnChanged(!spotState->HHmaskSHcurve);
                blendmaskSH->setEditedState(spotState->blendmaskSH ? Edited : UnEdited);
                radmaskSH->setEditedState(spotState->radmaskSH ? Edited : UnEdited);
                blurSHde->setEditedState(spotState->blurSHde ? Edited : UnEdited);
                inverssh->set_inconsistent(multiImage && !spotState->inverssh);
                chromaskSH->setEditedState(spotState->chromaskSH ? Edited : UnEdited);
                gammaskSH->setEditedState(spotState->gammaskSH ? Edited : UnEdited);
                slomaskSH->setEditedState(spotState->slomaskSH ? Edited : UnEdited);

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
                softradiusret->setEditedState(spotState->softradiusret ? Edited : UnEdited);

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
                clarityml->setEditedState(spotState->clarityml ? Edited : UnEdited);
                contresid->setEditedState(spotState->contresid ? Edited : UnEdited);
                blurcbdl->setEditedState(spotState->blurcbdl ? Edited : UnEdited);
                softradiuscb->setEditedState(spotState->softradiuscb ? Edited : UnEdited);
                blendmaskcb->setEditedState(spotState->blendmaskcb ? Edited : UnEdited);
                radmaskcb->setEditedState(spotState->radmaskcb ? Edited : UnEdited);
                chromaskcb->setEditedState(spotState->chromaskcb ? Edited : UnEdited);
                gammaskcb->setEditedState(spotState->gammaskcb ? Edited : UnEdited);
                slomaskcb->setEditedState(spotState->slomaskcb ? Edited : UnEdited);
                enacbMask->set_inconsistent(multiImage && !spotState->enacbMask);
                CCmaskcbshape->setUnChanged(!spotState->CCmaskcbcurve);
                LLmaskcbshape->setUnChanged(!spotState->LLmaskcbcurve);
                HHmaskcbshape->setUnChanged(!spotState->HHmaskcbcurve);

                // Denoise
                expdenoi->set_inconsistent(!spotState->expdenoi);
                noiselumf->setEditedState(spotState->noiselumf ? Edited : UnEdited);
                noiselumf0->setEditedState(spotState->noiselumf0 ? Edited : UnEdited);
                noiselumf2->setEditedState(spotState->noiselumf2 ? Edited : UnEdited);
                noiselumc->setEditedState(spotState->noiselumc ? Edited : UnEdited);
                noiselumdetail->setEditedState(spotState->noiselumdetail ? Edited : UnEdited);
                noiselequal->setEditedState(spotState->noiselequal ? Edited : UnEdited);
                noisechrof->setEditedState(spotState->noisechrof ? Edited : UnEdited);
                noisechroc->setEditedState(spotState->noisechroc ? Edited : UnEdited);
                noisechrodetail->setEditedState(spotState->noisechrodetail ? Edited : UnEdited);
                adjblur->setEditedState(spotState->adjblur ? Edited : UnEdited);
                bilateral->setEditedState(spotState->bilateral ? Edited : UnEdited);
                sensiden->setEditedState(spotState->sensiden ? Edited : UnEdited);
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
        HCurveEditorG->show();
        curvactiv->hide();
        qualitycurveMethod->show();
        labqualcurv->show();
        expmaskcol->show();
        structcol->show();
        blurcolde->show();
        softradiuscol->show();
        showmaskcolMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
        gridFrame->hide();
    } else if (invers->get_active()) {
        sensi->show();
        llCurveEditorG->show();
        HCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        labqualcurv->hide();
        expmaskcol->hide();
        softradiuscol->hide();
        structcol->hide();
        blurcolde->show();
        gridFrame->hide();
    } else {
        sensi->show();
        llCurveEditorG->show();
        HCurveEditorG->show();
        curvactiv->hide();
        qualitycurveMethod->show();
        labqualcurv->show();
        expmaskcol->show();
        structcol->show();
        blurcolde->show();
        gridFrame->show();
        softradiuscol->show();

        if (batchMode) {
            showmaskcolMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
        }
    }

    // Update Exposure GUI according to black adjuster state (to be compliant with adjusterChanged function)
    if (multiImage && inversex->get_inconsistent()) {
        sensiex->show();
        curveEditorG->show();
    //    maskexpFrame->show();
        structexp->show();
        blurexpde->show();
        softradiusexp->show();
        showmaskexpMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
    } else if (inversex->get_active()) {
        sensiex->show();
        curveEditorG->show();
    //    maskexpFrame->hide();
        structexp->hide();
        blurexpde->show();
        softradiusexp->hide();
    } else {
        sensiex->show();
        curveEditorG->show();
    //    maskexpFrame->show();
        structexp->show();
        blurexpde->show();
        softradiusexp->show();

        if (batchMode) {
            showmaskexpMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
        }
    }

    // Update SH GUI according to black adjuster state (to be compliant with adjusterChanged function)
    if (multiImage && inversex->get_inconsistent()) {
        sensihs->show();
     //   maskSHFrame->show();
        blurSHde->show();
        showmaskSHMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
    } else if (inverssh->get_active()) {
        sensihs->show();
    //    maskSHFrame->hide();
        blurSHde->show();
    } else {
        sensihs->show();
    //    maskSHFrame->show();
        blurSHde->show();

        if (batchMode) {
            showmaskSHMethod->hide(); // Being able to change Color & Light mask visibility is useless in batch mode
        }
    }


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


    // Update Retinex GUI according to inversret button state (to be compliant with inversretChanged function)
    if (multiImage && inversret->get_inconsistent()) {
        sensih->show();
        softradiusret->show();
    } else if (inversret->get_active()) {
        sensih->show();
        softradiusret->hide();
    } else {
        sensih->show();
        softradiusret->show();
    }

    // Update Sharpening GUI according to inverssha button state (to be compliant with inversshaChanged function)
    if (multiImage && inverssha->get_inconsistent()) {
        sensisha->show();
    } else if (inverssha->get_active()) {
        sensisha->show();
    } else {
        sensisha->show();
    }
}

void Locallab::autoOpenCurve()
{
    // printf("autoOpenCurve\n");

    // TODO autoOpenCurve only considers linearity state of selected spot curve
//    llshape->openIfNonlinear();
//    ccshape->openIfNonlinear();
//    LHshape->openIfNonlinear();
//    HHshape->openIfNonlinear();
//    CCmaskshape->openIfNonlinear();
//    LLmaskshape->openIfNonlinear();
//    HHmaskshape->openIfNonlinear();
//    CCmaskexpshape->openIfNonlinear();
//    LLmaskexpshape->openIfNonlinear();
}
