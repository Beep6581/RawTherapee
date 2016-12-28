/*
 *  This file is part of RawTherapee.
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
 *
 *  2014 Jacques Desmis <jdesmis@gmail.com>
 */

#include "wavelet.h"
#include <cmath>
#include "edit.h"
#include "guiutils.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;
extern Options options;

namespace
{

    GradientMilestone makeHsvGm(double position, float h, float s, float v)
    {
        float r;
        float g;
        float b;
        Color::hsv2rgb01(h, s, v, r, g, b);
        return GradientMilestone(position, r, g, b);
    }

    std::vector<GradientMilestone> makeWholeHueRange()
    {
        std::vector<GradientMilestone> res;
        res.reserve(7);

        for (int i = 0; i < 7; ++i) {
            const float x = static_cast<float>(i) / 6.0f;
            res.push_back(makeHsvGm(x, x, 0.5f, 0.5f));
        }

        return res;
    }

}

Wavelet::Wavelet() :
    FoldableToolPanel(this, "wavelet", M("TP_WAVELET_LABEL"), true, true),
    curveEditorG(new CurveEditorGroup(options.lastWaveletCurvesDir, M("TP_WAVELET_CONTEDIT"))),
    CCWcurveEditorG(new CurveEditorGroup(options.lastWaveletCurvesDir, M("TP_WAVELET_CCURVE"))),
    curveEditorRES(new CurveEditorGroup(options.lastWaveletCurvesDir)),
    curveEditorGAM(new CurveEditorGroup(options.lastWaveletCurvesDir)),
    separatorNeutral(Gtk::manage(new Gtk::HSeparator())),
    separatoredge(Gtk::manage(new Gtk::HSeparator())),
    opaCurveEditorG(new CurveEditorGroup(options.lastWaveletCurvesDir, M("TP_WAVELET_COLORT"))),
    opacityCurveEditorG(new CurveEditorGroup(options.lastWaveletCurvesDir, M("TP_WAVELET_OPACITY"))),
    opacityCurveEditorW(new CurveEditorGroup(options.lastWaveletCurvesDir, M("TP_WAVELET_OPACITYW"))),
    opacityCurveEditorWL(new CurveEditorGroup(options.lastWaveletCurvesDir, M("TP_WAVELET_OPACITYWL"))),
    median(Gtk::manage(new Gtk::CheckButton(M("TP_WAVELET_MEDI")))),
    medianlev(Gtk::manage(new Gtk::CheckButton(M("TP_WAVELET_MEDILEV")))),
    linkedg(Gtk::manage(new Gtk::CheckButton(M("TP_WAVELET_LINKEDG")))),
    cbenab(Gtk::manage(new Gtk::CheckButton(M("TP_WAVELET_CBENAB")))),
    lipst(Gtk::manage(new Gtk::CheckButton(M("TP_WAVELET_LIPST")))),
    avoid(Gtk::manage(new Gtk::CheckButton(M("TP_WAVELET_AVOID")))),
    tmr(Gtk::manage(new Gtk::CheckButton(M("TP_WAVELET_BALCHRO")))),
    neutralchButton(Gtk::manage(new Gtk::Button(M("TP_WAVELET_NEUTRAL")))),
    rescon(Gtk::manage(new Adjuster(M("TP_WAVELET_RESCON"), -100, 100, 1, 0))),
    resconH(Gtk::manage(new Adjuster(M("TP_WAVELET_RESCONH"), -100, 100, 1, 0))),
    reschro(Gtk::manage(new Adjuster(M("TP_WAVELET_RESCHRO"), -100, 100, 1, 0))),
    tmrs(Gtk::manage(new Adjuster(M("TP_WAVELET_TMSTRENGTH"), -1.0, 2.0, 0.01, 0.0))),
    gamma(Gtk::manage(new Adjuster(M("TP_WAVELET_COMPGAMMA"), 0.4, 2.0, 0.01, 1.0))),
    sup(Gtk::manage(new Adjuster(M("TP_WAVELET_SUPE"), -100, 350, 1, 0))),
    sky(Gtk::manage(new Adjuster(M("TP_WAVELET_SKY"), -100., 100.0, 1., 0.))),
    thres(Gtk::manage(new Adjuster(M("TP_WAVELET_LEVELS"), 4, 9, 1, 7))),//3
    chroma(Gtk::manage(new Adjuster(M("TP_WAVELET_CHRO"), 1, 9, 1, 5))),
    chro(Gtk::manage(new Adjuster(M("TP_WAVELET_CHR"), 0., 100., 1., 0.))),
    contrast(Gtk::manage(new Adjuster(M("TP_WAVELET_CONTRA"), -100, 100, 1, 0))),
    thr(Gtk::manage(new Adjuster(M("TP_WAVELET_THR"), 0, 100, 1, 35))),
    thrH(Gtk::manage(new Adjuster(M("TP_WAVELET_THRH"), 0, 100, 1, 65))),
    skinprotect(Gtk::manage( new Adjuster(M("TP_WAVELET_SKIN"), -100, 100, 1, 0.) )),
    edgrad(Gtk::manage( new Adjuster(M("TP_WAVELET_EDRAD"), 0, 100, 1, 15) )),
    edgval(Gtk::manage( new Adjuster(M("TP_WAVELET_EDVAL"), 0, 100, 1, 0) )),
    edgthresh(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGTHRESH"), -50, 100, 1, 10 ))),
    strength(Gtk::manage(new Adjuster(M("TP_WAVELET_STRENGTH"), 0, 100, 1, 100))),
    balance(Gtk::manage(new Adjuster(M("TP_WAVELET_BALANCE"), -30, 100, 1, 0))),
    iter(Gtk::manage(new Adjuster(M("TP_WAVELET_ITER"), -3, 3, 1, 0))),
    hueskin(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_HUESKIN"), -314., 314., -5., 25., 170., 120., 0, false))),
    hueskin2(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_HUESKY"), -314., 314., -260., -250, -130., -140., 0, false))),
    hllev(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_HIGHLIGHT"), 0., 100., 50., 75., 100., 98., 0, false))),
    bllev(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_LOWLIGHT"), 0., 100., 0., 2., 50., 25., 0, false))),
    pastlev(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_PASTEL"), 0., 70., 0., 2., 30., 20., 0, false))),
    satlev(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_SAT"), 0., 130., 30., 45., 130., 100., 0, false))),
    edgcont(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_EDGCONT"), 0., 100., options.rtSettings.bot_left, options.rtSettings.top_left, options.rtSettings.bot_right, options.rtSettings.top_right, 0., false))),
    level0noise(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_LEVZERO"), -30., 100., 0., M("TP_WAVELET_STREN"), 1., 0., 100., 0., M("TP_WAVELET_NOIS"), 1., nullptr, false))),
    level1noise(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_LEVONE"), -30., 100., 0., M("TP_WAVELET_STREN"), 1., 0., 100., 0., M("TP_WAVELET_NOIS"), 1., nullptr, false))),
    level2noise(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_LEVTWO"), -30., 100., 0., M("TP_WAVELET_STREN"), 1., 0., 100., 0., M("TP_WAVELET_NOIS"), 1., nullptr, false))),
    level3noise(Gtk::manage(new ThresholdAdjuster(M("TP_WAVELET_LEVTHRE"), -30., 100., 0., M("TP_WAVELET_STREN"), 1., 0., 100., 0., M("TP_WAVELET_NOIS"), 1., nullptr, false))),
    threshold(Gtk::manage(new Adjuster(M("TP_WAVELET_THRESHOLD"), 1, 9, 1, 5))),
    threshold2(Gtk::manage(new Adjuster(M("TP_WAVELET_THRESHOLD2"), 1, 9, 1, 4))),
    edgedetect(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGEDETECT"), 0, 100, 1, 90))),
    edgedetectthr(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGEDETECTTHR"), 0, 100, 1, 20))),
    edgedetectthr2(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGEDETECTTHR2"), -10, 100, 1, 0))),
    edgesensi(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGESENSI"), 0, 100, 1, 60))),
    edgeampli(Gtk::manage(new Adjuster(M("TP_WAVELET_EDGEAMPLI"), 0, 100, 1, 10))),
    Lmethod(Gtk::manage(new MyComboBoxText())),
    CHmethod(Gtk::manage(new MyComboBoxText())),
    CHSLmethod(Gtk::manage(new MyComboBoxText())),
    EDmethod(Gtk::manage(new MyComboBoxText())),
    BAmethod(Gtk::manage(new MyComboBoxText())),
    NPmethod(Gtk::manage(new MyComboBoxText())),
    TMmethod(Gtk::manage(new MyComboBoxText())),
    HSmethod(Gtk::manage(new MyComboBoxText())),
    CLmethod(Gtk::manage(new MyComboBoxText())),
    Backmethod(Gtk::manage(new MyComboBoxText())),
    Tilesmethod(Gtk::manage(new MyComboBoxText())),
    daubcoeffmethod(Gtk::manage(new MyComboBoxText())),
    Dirmethod(Gtk::manage(new MyComboBoxText())),
    Medgreinf(Gtk::manage(new MyComboBoxText())),
    chanMixerHLFrame(Gtk::manage(new Gtk::Frame(M("TP_COLORTONING_HIGHLIGHT")))),
    chanMixerMidFrame(Gtk::manage(new Gtk::Frame(M("TP_COLORTONING_MIDTONES")))),
    chanMixerShadowsFrame(Gtk::manage(new Gtk::Frame(M("TP_COLORTONING_SHADOWS")))),
    wavLabels(Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER))),
    labmC(Gtk::manage(new Gtk::Label(M("TP_WAVELET_CTYPE") + ":"))),
    labmNP(Gtk::manage(new Gtk::Label(M("TP_WAVELET_NPTYPE") + ":"))),
    expchroma(new MyExpander(true, M("TP_WAVELET_LEVCH"))),
    expcontrast(new MyExpander(true, M("TP_WAVELET_LEVF"))),
    expedge(new MyExpander(true, M("TP_WAVELET_EDGE"))),
    expfinal(new MyExpander(true, M("TP_WAVELET_FINAL"))),
    expgamut(new MyExpander(false, M("TP_WAVELET_CONTR"))),
    expnoise(new MyExpander(true, M("TP_WAVELET_NOISE"))),
    expresid(new MyExpander(true, M("TP_WAVELET_RESID"))),
    expsettings(new MyExpander(false, M("TP_WAVELET_SETTINGS"))),
    exptoning(new MyExpander(true, M("TP_WAVELET_TON"))),
    neutrHBox(Gtk::manage(new Gtk::HBox()))
{
    CurveListener::setMulti(true);
    nextnlevel = 7.;

    expsettings->signal_button_release_event().connect_notify( sigc::bind( sigc::mem_fun(this, &Wavelet::foldAllButMe), expsettings) );

    expcontrast->signal_button_release_event().connect_notify( sigc::bind( sigc::mem_fun(this, &Wavelet::foldAllButMe), expcontrast) );
    enableContrastConn = expcontrast->signal_enabled_toggled().connect ( sigc::bind( sigc::mem_fun(this, &Wavelet::enableToggled), expcontrast) );

    expchroma->signal_button_release_event().connect_notify( sigc::bind ( sigc::mem_fun(this, &Wavelet::foldAllButMe), expchroma) );
    enableChromaConn = expchroma->signal_enabled_toggled().connect ( sigc::bind( sigc::mem_fun(this, &Wavelet::enableToggled), expchroma) );

    exptoning->signal_button_release_event().connect_notify( sigc::bind ( sigc::mem_fun(this, &Wavelet::foldAllButMe), exptoning) );
    enableToningConn = exptoning->signal_enabled_toggled().connect ( sigc::bind( sigc::mem_fun(this, &Wavelet::enableToggled), exptoning) );

    expnoise->signal_button_release_event().connect_notify( sigc::bind ( sigc::mem_fun(this, &Wavelet::foldAllButMe), expnoise) );
    enableNoiseConn = expnoise->signal_enabled_toggled().connect ( sigc::bind( sigc::mem_fun(this, &Wavelet::enableToggled), expnoise) );

    expedge->signal_button_release_event().connect_notify( sigc::bind ( sigc::mem_fun(this, &Wavelet::foldAllButMe), expedge) );
    enableEdgeConn = expedge->signal_enabled_toggled().connect ( sigc::bind( sigc::mem_fun(this, &Wavelet::enableToggled), expedge) );

    expgamut->signal_button_release_event().connect_notify( sigc::bind ( sigc::mem_fun(this, &Wavelet::foldAllButMe), expgamut) );

    expresid->signal_button_release_event().connect_notify( sigc::bind ( sigc::mem_fun(this, &Wavelet::foldAllButMe), expresid) );
    enableResidConn = expresid->signal_enabled_toggled().connect ( sigc::bind( sigc::mem_fun(this, &Wavelet::enableToggled), expresid) );

    expfinal->signal_button_release_event().connect_notify( sigc::bind ( sigc::mem_fun(this, &Wavelet::foldAllButMe), expfinal) );
    enableFinalConn = expfinal->signal_enabled_toggled().connect ( sigc::bind( sigc::mem_fun(this, &Wavelet::enableToggled), expfinal) );

// Wavelet Settings
    Gtk::VBox* const settingsVBox = Gtk::manage(new Gtk::VBox());
    settingsVBox->set_border_width(4);
    settingsVBox->set_spacing(2);

    strength->setAdjusterListener (this);

    thres->set_tooltip_text (M("TP_WAVELET_LEVELS_TOOLTIP"));
    thres->setAdjusterListener (this);

    Tilesmethod->append_text (M("TP_WAVELET_TILESFULL"));
    Tilesmethod->append_text (M("TP_WAVELET_TILESBIG"));
    Tilesmethod->append_text (M("TP_WAVELET_TILESLIT"));
    Tilesmethodconn = Tilesmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::TilesmethodChanged) );
    Tilesmethod->set_tooltip_text (M("TP_WAVELET_TILES_TOOLTIP"));
    Gtk::HBox* const tilesizeHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const tilesizeLabel = Gtk::manage(new Gtk::Label(M("TP_WAVELET_TILESIZE") + ":"));
    tilesizeHBox->pack_start(*tilesizeLabel, Gtk::PACK_SHRINK, 4);
    tilesizeHBox->pack_start(*Tilesmethod);

    daubcoeffmethod->set_sensitive(true);
    daubcoeffmethod->append_text (M("TP_WAVELET_DAUB2"));
    daubcoeffmethod->append_text (M("TP_WAVELET_DAUB4"));
    daubcoeffmethod->append_text (M("TP_WAVELET_DAUB6"));
    daubcoeffmethod->append_text (M("TP_WAVELET_DAUB10"));
    daubcoeffmethod->append_text (M("TP_WAVELET_DAUB14"));
    daubcoeffmethodconn = daubcoeffmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::daubcoeffmethodChanged) );
    daubcoeffmethod->set_tooltip_text (M("TP_WAVELET_DAUB_TOOLTIP"));
    Gtk::Label* const daubcoeffLabel = Gtk::manage(new Gtk::Label(M("TP_WAVELET_DAUB") + ":"));
    Gtk::HBox* const daubcoeffHBox = Gtk::manage(new Gtk::HBox());
    daubcoeffHBox->pack_start(*daubcoeffLabel, Gtk::PACK_SHRINK, 4);
    daubcoeffHBox->pack_start(*daubcoeffmethod);

    Backmethod->append_text (M("TP_WAVELET_B0"));
    Backmethod->append_text (M("TP_WAVELET_B1"));
    Backmethod->append_text (M("TP_WAVELET_B2"));
    Backmethodconn = Backmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::BackmethodChanged) );
    Gtk::HBox* const backgroundHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const backgroundLabel = Gtk::manage(new Gtk::Label(M("TP_WAVELET_BACKGROUND") + ":"));
    backgroundHBox->pack_start(*backgroundLabel, Gtk::PACK_SHRINK, 4);
    backgroundHBox->pack_start(*Backmethod);

    CLmethod->append_text (M("TP_WAVELET_LEVDIR_ONE"));
    CLmethod->append_text (M("TP_WAVELET_LEVDIR_INF"));
    CLmethod->append_text (M("TP_WAVELET_LEVDIR_SUP"));
    CLmethod->append_text (M("TP_WAVELET_LEVDIR_ALL"));
    CLmethodconn = CLmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::CLmethodChanged) );
    Gtk::HBox* const levdirMainHBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const levdirMainLabel = Gtk::manage(new Gtk::Label(M("TP_WAVELET_PROC") + ":"));
    levdirMainHBox->pack_start(*levdirMainLabel, Gtk::PACK_SHRINK, 4);
    levdirMainHBox->pack_start(*CLmethod); //same

    Lmethod->set_sensitive(false);
    Lmethod->set_sensitive(false);
    Lmethod->append_text (M("TP_WAVELET_1"));
    Lmethod->append_text (M("TP_WAVELET_2"));
    Lmethod->append_text (M("TP_WAVELET_3"));
    Lmethod->append_text (M("TP_WAVELET_4"));
    Lmethod->append_text (M("TP_WAVELET_5"));
    Lmethod->append_text (M("TP_WAVELET_6"));
    Lmethod->append_text (M("TP_WAVELET_7"));
    Lmethod->append_text (M("TP_WAVELET_8"));
    Lmethod->append_text (M("TP_WAVELET_9"));
    Lmethod->append_text (M("TP_WAVELET_SUPE"));
    Lmethod->append_text (M("TP_WAVELET_RESID"));
    Lmethod->set_active(0);
    Dirmethod->set_sensitive(false);
    Dirmethod->append_text (M("TP_WAVELET_DONE"));
    Dirmethod->append_text (M("TP_WAVELET_DTWO"));
    Dirmethod->append_text (M("TP_WAVELET_DTHR"));
    Dirmethod->append_text (M("TP_WAVELET_DALL"));
    Lmethodconn = Lmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::LmethodChanged) );
    Dirmethodconn = Dirmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::DirmethodChanged) );
    Gtk::HBox* const levdirSubHBox = Gtk::manage(new Gtk::HBox());
    levdirSubHBox->pack_start(*Lmethod);
    levdirSubHBox->pack_start(*Dirmethod, Gtk::PACK_EXPAND_WIDGET, 2); // same, but 2 not 4?

    settingsVBox->pack_start(*strength);
    settingsVBox->pack_start(*thres);
    settingsVBox->pack_start(*tilesizeHBox);
    settingsVBox->pack_start(*daubcoeffHBox);
    settingsVBox->pack_start(*backgroundHBox);
    settingsVBox->pack_start(*levdirMainHBox);
    settingsVBox->pack_start(*levdirSubHBox);

// Contrast
    Gtk::VBox* const levBox = Gtk::manage (new Gtk::VBox());
    levBox->set_border_width(4);
    levBox->set_spacing(2);


    Gtk::HBox* const buttonBox = Gtk::manage (new Gtk::HBox(true, 10));
    levBox->pack_start(*buttonBox, Gtk::PACK_SHRINK, 2);

    Gtk::Button* const contrastMinusButton = Gtk::manage (new Gtk::Button(M("TP_WAVELET_CONTRAST_MINUS")));
    buttonBox->pack_start(*contrastMinusButton);
    contrastMinusPressedConn = contrastMinusButton->signal_pressed().connect( sigc::mem_fun(*this, &Wavelet::contrastMinusPressed));

    Gtk::Button* const neutralButton = Gtk::manage (new Gtk::Button(M("TP_WAVELET_NEUTRAL")));
    buttonBox->pack_start(*neutralButton);
    neutralPressedConn = neutralButton->signal_pressed().connect( sigc::mem_fun(*this, &Wavelet::neutralPressed));

    Gtk::Button* const contrastPlusButton = Gtk::manage (new Gtk::Button(M("TP_WAVELET_CONTRAST_PLUS")));
    buttonBox->pack_start(*contrastPlusButton);
    contrastPlusPressedConn = contrastPlusButton->signal_pressed().connect( sigc::mem_fun(*this, &Wavelet::contrastPlusPressed));

    buttonBox->show_all_children();

    for(int i = 0; i < 9; i++) {
        Glib::ustring ss;

        switch( i ) {
        case 0:
            ss = Glib::ustring::compose( "%1 (%2)", (i + 1), M("TP_WAVELET_FINEST"));
            break;

        case 8:
            ss = Glib::ustring::compose( "%1 (%2)", (i + 1), M("TP_WAVELET_LARGEST"));
            break;

        default:
            ss = Glib::ustring::compose( "%1", (i + 1));
        }

        correction[i] = Gtk::manage ( new Adjuster (ss, -100, 350, 1, 0) );
        correction[i]->setAdjusterListener(this);
        levBox->pack_start(*correction[i]);
    }

    levBox->pack_start(*sup);
    sup->setAdjusterListener (this);
    wavLabels->show();
    levBox->pack_start (*wavLabels);

    Gtk::VBox* const contrastSHVBox = Gtk::manage(new Gtk::VBox);
    contrastSHVBox->set_border_width(4);
    contrastSHVBox->set_spacing(2);

    HSmethod->append_text (M("TP_WAVELET_HS1"));
    HSmethod->append_text (M("TP_WAVELET_HS2"));
    HSmethodconn = HSmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::HSmethodChanged) );

    const std::vector<GradientMilestone> milestones2 = {
        GradientMilestone(0.0, 0.0, 0.0, 0.0),
        GradientMilestone(1.0, 1.0, 1.0, 1.0)
    };

    hllev->setAdjusterListener (this);
    hllev->setBgGradient(milestones2);

    threshold->setAdjusterListener (this);
    threshold->set_tooltip_text (M("TP_WAVELET_THRESHOLD_TOOLTIP"));

    bllev->setAdjusterListener (this);
    bllev->setBgGradient(milestones2);

    threshold2->setAdjusterListener (this);
    threshold2->set_tooltip_text (M("TP_WAVELET_THRESHOLD2_TOOLTIP"));

    contrastSHVBox->pack_start(*HSmethod); //remove 2?
    contrastSHVBox->pack_start(*hllev);
    contrastSHVBox->pack_start(*threshold);
    contrastSHVBox->pack_start(*bllev);
    contrastSHVBox->pack_start(*threshold2);
    Gtk::Frame* const contrastSHFrame = Gtk::manage(new Gtk::Frame(M("TP_WAVELET_APPLYTO")));
    contrastSHFrame->add(*contrastSHVBox);
    levBox->pack_start(*contrastSHFrame);

// Chromaticity
    Gtk::VBox* const chBox = Gtk::manage (new Gtk::VBox());
    chBox->set_border_width(4);
    chBox->set_spacing(2);

    Gtk::Label* const labmch = Gtk::manage(new Gtk::Label(M("TP_WAVELET_CHTYPE") + ":"));
    Gtk::HBox* const ctboxch = Gtk::manage(new Gtk::HBox());
    ctboxch->pack_start (*labmch, Gtk::PACK_SHRINK, 1);

    CHmethod->append_text (M("TP_WAVELET_CH1"));
    CHmethod->append_text (M("TP_WAVELET_CH2"));
    CHmethod->append_text (M("TP_WAVELET_CH3"));
    CHmethodconn = CHmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::CHmethodChanged) );
    ctboxch->pack_start(*CHmethod);
    chBox->pack_start(*ctboxch);

    Gtk::HBox* const ctboxCH = Gtk::manage(new Gtk::HBox());
    ctboxCH->pack_start (*labmC, Gtk::PACK_SHRINK, 1);

    CHSLmethod->append_text (M("TP_WAVELET_CHSL"));
    CHSLmethod->append_text (M("TP_WAVELET_CHCU"));
    CHSLmethodconn = CHSLmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::CHSLmethodChanged) );
    ctboxCH->pack_start(*CHSLmethod);

    Gtk::HSeparator* const separatorChromaMethod = Gtk::manage (new Gtk::HSeparator());
    chBox->pack_start(*separatorChromaMethod, Gtk::PACK_SHRINK, 2);

    chroma->set_tooltip_text (M("TP_WAVELET_CHRO_TOOLTIP"));
    chBox->pack_start(*chroma);
    chroma->setAdjusterListener (this);

    satlev->setAdjusterListener (this);
    satlev->setBgGradient(milestones2);

    pastlev->setAdjusterListener (this);
    pastlev->setBgGradient(milestones2);

    chBox->pack_start(*pastlev);
    chBox->pack_start(*satlev);

    chro->set_tooltip_text (M("TP_WAVELET_CHR_TOOLTIP"));
    chBox->pack_start(*chro);
    chro->setAdjusterListener (this);

    Gtk::HBox* const buttonchBox = Gtk::manage (new Gtk::HBox(true, 10));
    neutralchPressedConn = neutralchButton->signal_pressed().connect( sigc::mem_fun(*this, &Wavelet::neutralchPressed));
    chBox->pack_start(*separatorNeutral, Gtk::PACK_SHRINK, 2);
    buttonchBox->pack_start(*neutralchButton);
    buttonchBox->show_all_children();
    chBox->pack_start(*buttonchBox, Gtk::PACK_SHRINK, 2);

    for(int i = 0; i < 9; i++) {
        Glib::ustring ss;

        switch( i ) {
        case 0:
            ss = Glib::ustring::compose( "%1 (%2)", (i + 1), M("TP_WAVELET_FINEST"));
            break;

        case 8:
            ss = Glib::ustring::compose( "%1 (%2)", (i + 1), M("TP_WAVELET_LARGEST"));
            break;

        default:
            ss = Glib::ustring::compose( "%1", (i + 1));
        }

        correctionch[i] = Gtk::manage ( new Adjuster (ss, -100, 100, 1, 0) );
        correctionch[i]->setAdjusterListener(this);
        chBox->pack_start(*correctionch[i]);
    }

// Toning
    Gtk::VBox* const tonBox = Gtk::manage (new Gtk::VBox());
    tonBox->set_border_width(4);
    tonBox->set_spacing(2);

    opaCurveEditorG->setCurveListener (this);

    std::vector<double> defaultCurve;

    rtengine::WaveletParams::getDefaultOpacityCurveRG(defaultCurve);
    opacityShapeRG = static_cast<FlatCurveEditor*>(opaCurveEditorG->addCurve(CT_Flat, "", nullptr, false, false));
    opacityShapeRG->setIdentityValue(0.);
    opacityShapeRG->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);

    opaCurveEditorG->curveListComplete();
    opaCurveEditorG->show();

    tonBox->pack_start( *opaCurveEditorG, Gtk::PACK_SHRINK, 2);

    opacityCurveEditorG->setCurveListener (this);

    rtengine::WaveletParams::getDefaultOpacityCurveBY(defaultCurve);
    opacityShapeBY = static_cast<FlatCurveEditor*>(opacityCurveEditorG->addCurve(CT_Flat, "", nullptr, false, false));
    opacityShapeBY->setIdentityValue(0.);
    opacityShapeBY->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);

    opacityCurveEditorG->curveListComplete();
    opacityCurveEditorG->show();

    tonBox->pack_start( *opacityCurveEditorG, Gtk::PACK_SHRINK, 2);

// Denoise and Refine
    Gtk::VBox* const noiseBox = Gtk::manage (new Gtk::VBox());
    noiseBox->set_border_width(4);
    noiseBox->set_spacing(2);

    linkedg->set_active (true);
    linkedgConn = linkedg->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::linkedgToggled) );
    noiseBox->pack_start(*linkedg);

    level0noise->setAdjusterListener (this);
    level0noise->setUpdatePolicy(RTUP_DYNAMIC);

    level1noise->setAdjusterListener (this);
    level1noise->setUpdatePolicy(RTUP_DYNAMIC);

    level2noise->setAdjusterListener (this);
    level2noise->setUpdatePolicy(RTUP_DYNAMIC);

    level3noise->setAdjusterListener (this);
    level3noise->setUpdatePolicy(RTUP_DYNAMIC);

    noiseBox->pack_start( *level0noise, Gtk::PACK_SHRINK, 0);
    noiseBox->pack_start( *level1noise, Gtk::PACK_SHRINK, 0);
    noiseBox->pack_start( *level2noise, Gtk::PACK_SHRINK, 0);
    noiseBox->pack_start( *level3noise, Gtk::PACK_SHRINK, 0);

// Edge Sharpness
    Gtk::VBox* const edgBox = Gtk::manage (new Gtk::VBox());
    edgBox->set_border_width(4);
    edgBox->set_spacing(2);

    edgval->setAdjusterListener(this);
    edgBox->pack_start(*edgval);

    edgrad->setAdjusterListener(this);
    edgBox->pack_start(*edgrad);
    edgrad->set_tooltip_markup (M("TP_WAVELET_EDRAD_TOOLTIP"));

    edgthresh->setAdjusterListener (this);
    edgthresh->set_tooltip_markup (M("TP_WAVELET_EDGTHRESH_TOOLTIP"));
    edgBox->pack_start (*edgthresh);

    Gtk::Label* const labmedgr = Gtk::manage(new Gtk::Label(M("TP_WAVELET_MEDGREINF") + ":"));
    Gtk::HBox* const edbox = Gtk::manage(new Gtk::HBox());
    edbox->pack_start (*labmedgr, Gtk::PACK_SHRINK, 1);

    Medgreinf->append_text (M("TP_WAVELET_RE1"));
    Medgreinf->append_text (M("TP_WAVELET_RE2"));
    Medgreinf->append_text (M("TP_WAVELET_RE3"));
    MedgreinfConn = Medgreinf->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::MedgreinfChanged) );
    Medgreinf->set_tooltip_markup (M("TP_WAVELET_EDGREINF_TOOLTIP"));
    edbox->pack_start(*Medgreinf);
    edgBox->pack_start(*edbox);

    Gtk::HSeparator* const separatorlc = Gtk::manage (new  Gtk::HSeparator());
    edgBox->pack_start(*separatorlc, Gtk::PACK_SHRINK, 2);

    Gtk::Label* const labmED = Gtk::manage(new Gtk::Label(M("TP_WAVELET_EDTYPE") + ":"));
    Gtk::HBox* const ctboxED = Gtk::manage(new Gtk::HBox());
    ctboxED->pack_start (*labmED, Gtk::PACK_SHRINK, 1);

    EDmethod->append_text (M("TP_WAVELET_EDSL"));
    EDmethod->append_text (M("TP_WAVELET_EDCU"));
    EDmethodconn = EDmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::EDmethodChanged) );
    ctboxED->pack_start(*EDmethod);
    edgBox->pack_start (*ctboxED);

    edgcont->setAdjusterListener (this);
    edgcont->setBgGradient(milestones2);
    edgcont->set_tooltip_markup (M("TP_WAVELET_EDGCONT_TOOLTIP"));

    // <-- Edge Sharpness  Local Contrast curve
    CCWcurveEditorG->setCurveListener (this);

    rtengine::WaveletParams::getDefaultCCWCurve(defaultCurve);
    ccshape = static_cast<FlatCurveEditor*>(CCWcurveEditorG->addCurve(CT_Flat, "", nullptr, false, false));

    ccshape->setIdentityValue(0.);
    ccshape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
    ccshape->setTooltip(M("TP_WAVELET_CURVEEDITOR_CC_TOOLTIP"));

    CCWcurveEditorG->curveListComplete();
    CCWcurveEditorG->show();
    // -->

    edgBox->pack_start (*edgcont);
    edgBox->pack_start(*CCWcurveEditorG, Gtk::PACK_SHRINK, 4);

    medianlev->set_active (true);
    medianlevConn = medianlev->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::medianlevToggled) );
    medianlev->set_tooltip_text (M("TP_WAVELET_MEDILEV_TOOLTIP"));

    Gtk::HSeparator* const separatored1 = Gtk::manage (new  Gtk::HSeparator());
    edgBox->pack_start(*separatored1, Gtk::PACK_SHRINK, 2);

    Gtk::HBox* const eddebox = Gtk::manage(new Gtk::HBox());
    edgBox->pack_start (*eddebox);
    edgBox->pack_start(*medianlev);

    edgedetect->setAdjusterListener (this);
    edgedetect->set_tooltip_text (M("TP_WAVELET_EDGEDETECT_TOOLTIP"));
    edgBox->pack_start(*edgedetect);

    edgedetectthr->setAdjusterListener (this);
    edgedetectthr->set_tooltip_text (M("TP_WAVELET_EDGEDETECTTHR_TOOLTIP"));
    edgBox->pack_start(*edgedetectthr);


    edgedetectthr2->setAdjusterListener (this);
    edgBox->pack_start(*edgedetectthr2);


    edgBox->pack_start(*separatoredge, Gtk::PACK_SHRINK, 2);

    lipst->set_active (true);
    lipstConn = lipst->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::lipstToggled) );
//  lipst->set_tooltip_text (M("TP_WAVELET_LIPST_TOOLTIP"));
    edgBox->pack_start(*lipst);

    edgesensi->setAdjusterListener (this);
    edgBox->pack_start(*edgesensi);


    edgeampli->setAdjusterListener (this);
    edgBox->pack_start(*edgeampli);


    Gtk::VBox* const ctboxES = Gtk::manage (new Gtk::VBox());

    ctboxES->set_border_width(4);
    ctboxES->set_spacing(2);

    Gtk::HBox* const ctboxNP = Gtk::manage(new Gtk::HBox());
    ctboxNP->pack_start (*labmNP, Gtk::PACK_SHRINK, 1);

    NPmethod->append_text (M("TP_WAVELET_NPNONE"));
    NPmethod->append_text (M("TP_WAVELET_NPLOW"));
    NPmethod->append_text (M("TP_WAVELET_NPHIGH"));
    NPmethodconn = NPmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::NPmethodChanged) );
    NPmethod->set_tooltip_text (M("TP_WAVELET_NPTYPE_TOOLTIP"));

    ctboxNP->pack_start(*NPmethod);
    ctboxES->pack_start(*ctboxNP);

    edgBox->pack_start(*ctboxES);

// Gamut
    Gtk::VBox* const conBox = Gtk::manage (new Gtk::VBox());
    conBox->set_border_width(4);
    conBox->set_spacing(2);

    median->set_active (true);
    medianConn = median->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::medianToggled) );
    conBox->pack_start(*median);

    hueskin->set_tooltip_markup (M("TP_WAVELET_HUESKIN_TOOLTIP"));

    //from -PI to +PI (radians) convert to hsv and draw bottombar
    const std::vector<GradientMilestone> milestones = {
        makeHsvGm(0.0000, 0.4199f, 0.5f, 0.5f), // hsv: 0.4199 rad: -3.14
        makeHsvGm(0.0540, 0.5000f, 0.5f, 0.5f), // hsv: 0.5    rad: -2.8
        makeHsvGm(0.1336, 0.6000f, 0.5f, 0.5f), // hsv: 0.60   rad: -2.3
        makeHsvGm(0.3567, 0.7500f, 0.5f, 0.5f), // hsv: 0.75   rad: -0.9
        makeHsvGm(0.4363, 0.8560f, 0.5f, 0.5f), // hsv: 0.856  rad: -0.4
        makeHsvGm(0.4841, 0.9200f, 0.5f, 0.5f), // hsv: 0.92   rad: -0.1
        makeHsvGm(0.5000, 0.9300f, 0.5f, 0.5f), // hsv: 0.93   rad:  0
        makeHsvGm(0.5366, 0.9600f, 0.5f, 0.5f), // hsv: 0.96   rad:  0.25
        makeHsvGm(0.5955, 1.0000f, 0.5f, 0.5f), // hsv: 1.     rad:  0.6
        makeHsvGm(0.6911, 0.0675f, 0.5f, 0.5f), // hsv: 0.0675 rad:  1.2
        makeHsvGm(0.7229, 0.0900f, 0.5f, 0.5f), // hsv: 0.09   rad:  1.4
        makeHsvGm(0.7707, 0.1700f, 0.5f, 0.5f), // hsv: 0.17   rad:  1.7
        makeHsvGm(0.8503, 0.2650f, 0.5f, 0.5f), // hsv: 0.265  rad:  2.1
        makeHsvGm(0.8981, 0.3240f, 0.5f, 0.5f), // hsv: 0.324  rad:  2.5
        makeHsvGm(1.0000, 0.4197f, 0.5f, 0.5f)  // hsv: 0.419  rad:  3.14
    };

    hueskin->setBgGradient(milestones);
    conBox->pack_start(*hueskin);
    hueskin->setAdjusterListener (this);

    skinprotect->setAdjusterListener(this);
    conBox->pack_start(*skinprotect);
    skinprotect->set_tooltip_markup (M("TP_WAVELET_SKIN_TOOLTIP"));

    curveEditorGAM->setCurveListener (this);

    Chshape = static_cast<FlatCurveEditor*>(curveEditorGAM->addCurve(CT_Flat, M("TP_WAVELET_CURVEEDITOR_CH")));
    Chshape->setTooltip(M("TP_WAVELET_CURVEEDITOR_CH_TOOLTIP"));
    Chshape->setCurveColorProvider(this, 5);
    curveEditorGAM->curveListComplete();
    Chshape->setBottomBarBgGradient(milestones);

    conBox->pack_start (*curveEditorGAM, Gtk::PACK_SHRINK, 4);

    avoid->set_active (true);
    avoidConn = avoid->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::avoidToggled) );
    conBox->pack_start(*avoid);

// Residual Image
    Gtk::VBox* const resBox = Gtk::manage (new Gtk::VBox());
    resBox->set_border_width(4);
    resBox->set_spacing(2);

    rescon->setAdjusterListener (this);
    resBox->pack_start(*rescon, Gtk::PACK_SHRINK);

    resBox->pack_start(*thr);
    thr->setAdjusterListener (this);

    resconH->setAdjusterListener (this);
    resBox->pack_start(*resconH, Gtk::PACK_SHRINK);

    thrH->setAdjusterListener (this);
    resBox->pack_start(*thrH, Gtk::PACK_SHRINK);

    contrast->set_tooltip_text (M("TP_WAVELET_CONTRA_TOOLTIP"));
    contrast->setAdjusterListener (this);
    resBox->pack_start(*contrast);  //keep the possibility to reinstall

    reschro->setAdjusterListener (this);
    resBox->pack_start(*reschro);


    Gtk::Label* const labmTM = Gtk::manage(new Gtk::Label(M("TP_WAVELET_TMTYPE") + ":"));
    Gtk::HBox* const ctboxTM = Gtk::manage(new Gtk::HBox());
    ctboxTM->pack_start (*labmTM, Gtk::PACK_SHRINK, 1);

    Gtk::HSeparator* const separatorR0 = Gtk::manage (new  Gtk::HSeparator());
    resBox->pack_start(*separatorR0, Gtk::PACK_SHRINK, 2);

    TMmethod->append_text (M("TP_WAVELET_COMPCONT"));
    TMmethod->append_text (M("TP_WAVELET_COMPTM"));
    TMmethodconn = TMmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::TMmethodChanged) );
    ctboxTM->pack_start(*TMmethod);
    resBox->pack_start (*ctboxTM);

    tmrs->set_tooltip_text (M("TP_WAVELET_TMSTRENGTH_TOOLTIP"));

    resBox->pack_start(*tmrs);
    tmrs->setAdjusterListener (this);

    gamma->set_tooltip_text (M("TP_WAVELET_COMPGAMMA_TOOLTIP"));
    resBox->pack_start(*gamma);
    gamma->setAdjusterListener (this);

    Gtk::HSeparator* const separatorR1 = Gtk::manage (new  Gtk::HSeparator());
    resBox->pack_start(*separatorR1, Gtk::PACK_SHRINK, 2);

    hueskin2->set_tooltip_markup (M("TP_WAVELET_HUESKY_TOOLTIP"));
    hueskin2->setBgGradient(milestones);
    resBox->pack_start(*hueskin2);
    hueskin2->setAdjusterListener (this);

    sky->set_tooltip_text (M("TP_WAVELET_SKY_TOOLTIP"));
    sky->setAdjusterListener (this);

    resBox->pack_start(*sky);

    // whole hue range
    const std::vector<GradientMilestone> milestones3 = makeWholeHueRange();

    curveEditorRES->setCurveListener (this);

    hhshape = static_cast<FlatCurveEditor*>(curveEditorRES->addCurve(CT_Flat, M("TP_WAVELET_CURVEEDITOR_HH")));
    hhshape->setTooltip(M("TP_WAVELET_CURVEEDITOR_HH_TOOLTIP"));
    hhshape->setCurveColorProvider(this, 5);
    curveEditorRES->curveListComplete();
    hhshape->setBottomBarBgGradient(milestones3);

    resBox->pack_start (*curveEditorRES, Gtk::PACK_SHRINK, 4);

    // Toning and Color Balance
    Gtk::HSeparator* const separatorCB = Gtk::manage (new  Gtk::HSeparator());

    Gtk::VBox* const chanMixerHLBox = Gtk::manage (new Gtk::VBox());
    Gtk::VBox* const chanMixerMidBox = Gtk::manage (new Gtk::VBox());
    Gtk::VBox* const chanMixerShadowsBox = Gtk::manage (new Gtk::VBox());

    cbenab->set_active (true);
    cbenabConn = cbenab->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::cbenabToggled) );
    cbenab->set_tooltip_text (M("TP_WAVELET_CB_TOOLTIP"));

    Gtk::Image* const iblueR   = Gtk::manage (new RTImage ("ajd-wb-temp1.png"));
    Gtk::Image* const iyelL    = Gtk::manage (new RTImage ("ajd-wb-temp2.png"));
    Gtk::Image* const imagL    = Gtk::manage (new RTImage ("ajd-wb-green1.png"));
    Gtk::Image* const igreenR  = Gtk::manage (new RTImage ("ajd-wb-green2.png"));

    Gtk::Image* const  iblueRm  = Gtk::manage (new RTImage ("ajd-wb-temp1.png"));
    Gtk::Image* const  iyelLm   = Gtk::manage (new RTImage ("ajd-wb-temp2.png"));
    Gtk::Image* const  imagLm   = Gtk::manage (new RTImage ("ajd-wb-green1.png"));
    Gtk::Image* const  igreenRm = Gtk::manage (new RTImage ("ajd-wb-green2.png"));

    Gtk::Image* const iblueRh  = Gtk::manage (new RTImage ("ajd-wb-temp1.png"));
    Gtk::Image* const iyelLh   = Gtk::manage (new RTImage ("ajd-wb-temp2.png"));
    Gtk::Image* const imagLh   = Gtk::manage (new RTImage ("ajd-wb-green1.png"));
    Gtk::Image* const igreenRh = Gtk::manage (new RTImage ("ajd-wb-green2.png"));

    greenhigh = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., igreenRh, imagLh));
    bluehigh = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., iblueRh, iyelLh));
    greenmed = Gtk::manage (new Adjuster ("", -100., 100., 1., 0.,  igreenRm, imagLm));
    bluemed = Gtk::manage (new Adjuster ("", -100., 100., 1., 0. , iblueRm , iyelLm));
    greenlow = Gtk::manage (new Adjuster ("", -100., 100., 1., 0.,  igreenR, imagL));
    bluelow = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., iblueR , iyelL));

    chanMixerHLBox->pack_start (*greenhigh);
    chanMixerHLBox->pack_start (*bluehigh);
    chanMixerMidBox->pack_start (*greenmed);
    chanMixerMidBox->pack_start (*bluemed);
    chanMixerShadowsBox->pack_start (*greenlow);
    chanMixerShadowsBox->pack_start (*bluelow);

    greenlow->setAdjusterListener (this);
    bluelow->setAdjusterListener (this);
    greenmed->setAdjusterListener (this);
    bluemed->setAdjusterListener (this);
    greenhigh->setAdjusterListener (this);
    bluehigh->setAdjusterListener (this);

    resBox->pack_start(*separatorCB, Gtk::PACK_SHRINK);

    chanMixerHLFrame->add(*chanMixerHLBox);
    chanMixerMidFrame->add(*chanMixerMidBox);
    chanMixerShadowsFrame->add(*chanMixerShadowsBox);
    resBox->pack_start(*cbenab);

    resBox->pack_start(*chanMixerHLFrame, Gtk::PACK_SHRINK);
    resBox->pack_start(*chanMixerMidFrame, Gtk::PACK_SHRINK);
    resBox->pack_start(*chanMixerShadowsFrame, Gtk::PACK_SHRINK);

    // Reset sliders
    neutrHBox->set_border_width (2);

    //RTImage *resetImg = Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
    //neutral->set_image(*resetImg);
    Gtk::Button* const neutral = Gtk::manage(new Gtk::Button(M("TP_COLORTONING_NEUTRAL")));
    neutral->set_tooltip_text (M("TP_COLORTONING_NEUTRAL_TIP"));
    neutralconn = neutral->signal_pressed().connect( sigc::mem_fun(*this, &Wavelet::neutral_pressed) );
    neutral->show();
    neutrHBox->pack_start (*neutral, Gtk::PACK_EXPAND_WIDGET);

    resBox->pack_start (*neutrHBox);

// Final Touchup
    Gtk::VBox* const ctboxBA = Gtk::manage (new Gtk::VBox());

    ctboxBA->set_border_width(4);
    ctboxBA->set_spacing(2);

    //Gtk::HSeparator *separatorfin = Gtk::manage (new  Gtk::HSeparator());
    //ctboxBA->pack_start(*separatorfin, Gtk::PACK_SHRINK, 2);
    Gtk::Label* const labmBA = Gtk::manage(new Gtk::Label(M("TP_WAVELET_BATYPE") + ":"));
    Gtk::HBox* const ctboxFI = Gtk::manage(new Gtk::HBox());
    ctboxFI->pack_start (*labmBA, Gtk::PACK_SHRINK, 1);

    BAmethod->append_text (M("TP_WAVELET_BANONE"));
    BAmethod->append_text (M("TP_WAVELET_BASLI"));
    BAmethod->append_text (M("TP_WAVELET_BACUR"));
    BAmethodconn = BAmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::BAmethodChanged) );
    ctboxFI->pack_start(*BAmethod);
    ctboxBA->pack_start(*ctboxFI);

    balance->setAdjusterListener (this);
    balance->set_tooltip_text (M("TP_WAVELET_BALANCE_TOOLTIP"));

    opacityCurveEditorW->setCurveListener (this);

    rtengine::WaveletParams::getDefaultOpacityCurveW(defaultCurve);
    opacityShape = static_cast<FlatCurveEditor*>(opacityCurveEditorW->addCurve(CT_Flat, "", nullptr, false, false));
    opacityShape->setIdentityValue(0.);
    opacityShape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
    opacityShape->setBottomBarBgGradient(milestones2);

    // This will add the reset button at the end of the curveType buttons
    opacityCurveEditorW->curveListComplete();
    opacityCurveEditorW->show();

    iter->setAdjusterListener (this);
    iter->set_tooltip_text (M("TP_WAVELET_ITER_TOOLTIP"));

    Gtk::HSeparator* const separatorbalend = Gtk::manage (new  Gtk::HSeparator());

    opacityCurveEditorWL->setCurveListener (this);

    rtengine::WaveletParams::getDefaultOpacityCurveWL(defaultCurve);
    opacityShapeWL = static_cast<FlatCurveEditor*>(opacityCurveEditorWL->addCurve(CT_Flat, "", nullptr, false, false));
    opacityShapeWL->setIdentityValue(0.);
    opacityShapeWL->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
    opacityShapeWL->setTooltip(M("TP_WAVELET_OPACITYWL_TOOLTIP"));

    // This will add the reset button at the end of the curveType buttons
    opacityCurveEditorWL->curveListComplete();
    opacityCurveEditorWL->show();

    curveEditorG->setCurveListener (this);

    clshape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, M("TP_WAVELET_CURVEEDITOR_CL")));
    clshape->setTooltip(M("TP_WAVELET_CURVEEDITOR_CL_TOOLTIP"));
    clshape->setBottomBarBgGradient(milestones2);
    clshape->setLeftBarBgGradient(milestones2);

    curveEditorG->curveListComplete();

    tmr->set_active (true);
    tmr->set_tooltip_text (M("TP_WAVELET_BALCHRO_TOOLTIP"));
    tmrConn = tmr->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::tmrToggled) );

    Gtk::VBox* const finalBox = Gtk::manage (new Gtk::VBox());
    finalBox->set_border_width(4);
    finalBox->set_spacing(2);

    finalBox->pack_start (*ctboxBA);
    finalBox->pack_start(*balance);

    finalBox->pack_start( *opacityCurveEditorW, Gtk::PACK_SHRINK, 2);

    finalBox->pack_start(*iter);

    finalBox->pack_start(*tmr);
    finalBox->pack_start(*separatorbalend, Gtk::PACK_SHRINK, 2);
    finalBox->pack_start( *opacityCurveEditorWL, Gtk::PACK_SHRINK, 2);

    finalBox->pack_start (*curveEditorG, Gtk::PACK_SHRINK, 4);

//-----------------------------

    expsettings->add(*settingsVBox);
    pack_start (*expsettings);

    expcontrast->add(*levBox);
    pack_start (*expcontrast);

    expchroma->add(*chBox);
    pack_start (*expchroma);

    exptoning->add(*tonBox);
    pack_start (*exptoning);

    expnoise->add(*noiseBox);
    pack_start (*expnoise);

    expedge->add(*edgBox);
    pack_start (*expedge);

    expgamut->add(*conBox);
    pack_start (*expgamut);

    expresid->add(*resBox);
    pack_start(*expresid);

    expfinal->add(*finalBox);
    pack_start(*expfinal);
}

Wavelet::~Wavelet ()
{
    delete opaCurveEditorG;
    delete opacityCurveEditorG;
    delete CCWcurveEditorG;
    delete curveEditorRES;
    delete curveEditorGAM;
    delete curveEditorG;
    delete opacityCurveEditorW;
    delete opacityCurveEditorWL;
}

int wavUpdateUI (void* data)
{
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
    (static_cast<Wavelet*>(data))->wavComputed_ ();
    return 0;
}

void Wavelet::wavChanged (double nlevel)
{
    nextnlevel = nlevel;
    g_idle_add (wavUpdateUI, this);
}
bool Wavelet::wavComputed_ ()
{

    disableListener ();
    enableListener ();
    updatewavLabel ();
    return false;
}
void Wavelet::updatewavLabel ()
{
    if (!batchMode) {
        float lv;
        lv = nextnlevel;
        wavLabels->set_text(
            Glib::ustring::compose(M("TP_WAVELET_LEVLABEL"),
                                   Glib::ustring::format(std::fixed, std::setprecision(0), lv))
        );
    }
}

// Will only reset the chanel mixer
// WARNING!  In mutiImage mode, and for sliders in ADD mode, this will reset the slider to 0, but not to the default value as in SET mode.
void Wavelet::neutral_pressed ()
{
    disableListener();
    greenlow->resetValue(false);
    bluelow->resetValue(false);
    greenmed->resetValue(false);
    bluemed->resetValue(false);
    greenhigh->resetValue(false);
    bluehigh->resetValue(false);

    enableListener();

    if (listener && getEnabled()) {
        listener->panelChanged (EvWavNeutral, M("ADJUSTER_RESET_TO_DEFAULT"));
    }
}

void Wavelet::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    /*****************************************************************************************************
     *
     *                                        Disable listeners and signals
     *
     *****************************************************************************************************/

    disableListener ();
    Lmethodconn.block(true);
    CLmethodconn.block(true);
    Backmethodconn.block(true);
    Tilesmethodconn.block(true);
    daubcoeffmethodconn.block(true);
    Dirmethodconn.block(true);
    CHmethodconn.block(true);
    CHSLmethodconn.block(true);
    EDmethodconn.block(true);
    NPmethodconn.block(true);
    BAmethodconn.block(true);
    TMmethodconn.block(true);
    HSmethodconn.block(true);
    MedgreinfConn.block(true);
    cbenabConn.block(true);
    enableChromaConn.block(true);
    enableContrastConn.block(true);
    enableEdgeConn.block(true);
    enableFinalConn.block(true);
    enableNoiseConn.block(true);
    enableResidConn.block(true);
    enableToningConn.block(true);

    /*****************************************************************************************************
     *
     *                               Set the GUI reflecting the given ProcParams
     *
     *****************************************************************************************************/

    //HSmethod->set_active (1);   // Note: default values are controlled in rtengine::ProcParams::SetDefaults
    if (pp->wavelet.HSmethod == "without") {
        HSmethod->set_active (0);
    } else if (pp->wavelet.HSmethod == "with") {
        HSmethod->set_active (1);
    }

    //CHmethod->set_active (1);
    if (pp->wavelet.CHmethod == "without") {
        CHmethod->set_active (0);
    } else if (pp->wavelet.CHmethod == "with") {
        CHmethod->set_active (1);
    } else if (pp->wavelet.CHmethod == "link") {
        CHmethod->set_active (2);
    }

    //Medgreinf->set_active (1);
    if (pp->wavelet.Medgreinf == "more") {
        Medgreinf->set_active (0);
    } else if (pp->wavelet.Medgreinf == "none") {
        Medgreinf->set_active (1);
    } else if (pp->wavelet.Medgreinf == "less") {
        Medgreinf->set_active (2);
    }

    //CHSLmethod->set_active (1);
    if (pp->wavelet.CHSLmethod == "SL") {
        CHSLmethod->set_active (0);
    } else if (pp->wavelet.CHSLmethod == "CU") {
        CHSLmethod->set_active (1);
    }

    //EDmethod->set_active (1);
    if (pp->wavelet.EDmethod == "SL") {
        EDmethod->set_active (0);
    } else if (pp->wavelet.EDmethod == "CU") {
        EDmethod->set_active (1);
    }

    if (pp->wavelet.NPmethod == "none") {
        NPmethod->set_active (0);
    } else if (pp->wavelet.NPmethod == "low") {
        NPmethod->set_active (1);
    } else if (pp->wavelet.NPmethod == "high") {
        NPmethod->set_active (2);
    }

    //BAmethod->set_active (0);
    if (pp->wavelet.BAmethod == "none") {
        BAmethod->set_active (0);
    } else if (pp->wavelet.BAmethod == "sli") {
        BAmethod->set_active (1);
    } else if (pp->wavelet.BAmethod == "cur") {
        BAmethod->set_active (2);
    }

    //TMmethod->set_active (1);
    if (pp->wavelet.TMmethod == "cont") {
        TMmethod->set_active (0);
    } else if (pp->wavelet.TMmethod == "tm") {
        TMmethod->set_active (1);
    }

//  else if (pp->wavelet.TMmethod=="both")
//      TMmethod->set_active (2);

    //Backmethod->set_active (3);
    if (pp->wavelet.Backmethod == "black") {
        Backmethod->set_active (0);
    } else if (pp->wavelet.Backmethod == "grey") {
        Backmethod->set_active (1);
    } else if (pp->wavelet.Backmethod == "resid") {
        Backmethod->set_active (2);
    }

    //CLmethod->set_active (3);
    if (pp->wavelet.CLmethod == "one") {
        CLmethod->set_active (0);
    } else if (pp->wavelet.CLmethod == "inf") {
        CLmethod->set_active (1);
    } else if (pp->wavelet.CLmethod == "sup") {
        CLmethod->set_active (2);
    } else if (pp->wavelet.CLmethod == "all") {
        CLmethod->set_active (3);
    }

    //Tilesmethod->set_active (2);
    if (pp->wavelet.Tilesmethod == "full") {
        Tilesmethod->set_active (0);
    } else if (pp->wavelet.Tilesmethod == "big") {
        Tilesmethod->set_active (1);
    } else if (pp->wavelet.Tilesmethod == "lit") {
        Tilesmethod->set_active (2);
    }

    //daubcoeffmethod->set_active (4);
    if (pp->wavelet.daubcoeffmethod == "2_") {
        daubcoeffmethod->set_active (0);
    } else if (pp->wavelet.daubcoeffmethod == "4_") {
        daubcoeffmethod->set_active (1);
    } else if (pp->wavelet.daubcoeffmethod == "6_") {
        daubcoeffmethod->set_active (2);
    } else if (pp->wavelet.daubcoeffmethod == "10_") {
        daubcoeffmethod->set_active (3);
    } else if (pp->wavelet.daubcoeffmethod == "14_") {
        daubcoeffmethod->set_active (4);
    }

    //Dirmethod->set_active (3);
    if (pp->wavelet.Dirmethod == "one") {
        Dirmethod->set_active (0);
    } else if (pp->wavelet.Dirmethod == "two") {
        Dirmethod->set_active (1);
    } else if (pp->wavelet.Dirmethod == "thr") {
        Dirmethod->set_active (2);
    } else if (pp->wavelet.Dirmethod == "all") {
        Dirmethod->set_active (3);
    }

    int selectedLevel = atoi(pp->wavelet.Lmethod.data()) - 1;
    Lmethod->set_active (selectedLevel == -1 ? 4 : selectedLevel);

    ccshape->setCurve (pp->wavelet.ccwcurve);
    opacityShapeRG->setCurve (pp->wavelet.opacityCurveRG);
    opacityShapeBY->setCurve (pp->wavelet.opacityCurveBY);
    opacityShape->setCurve (pp->wavelet.opacityCurveW);
    opacityShapeWL->setCurve (pp->wavelet.opacityCurveWL);
    hhshape->setCurve (pp->wavelet.hhcurve);
    Chshape->setCurve (pp->wavelet.Chcurve);
    clshape->setCurve (pp->wavelet.wavclCurve);
    expcontrast->setEnabled (pp->wavelet.expcontrast);
    expchroma->setEnabled (pp->wavelet.expchroma);
    expedge->setEnabled (pp->wavelet.expedge);
    expresid->setEnabled (pp->wavelet.expresid);
    expfinal->setEnabled (pp->wavelet.expfinal);
    exptoning->setEnabled (pp->wavelet.exptoning);
    expnoise->setEnabled (pp->wavelet.expnoise);

    setEnabled(pp->wavelet.enabled);
    avoidConn.block (true);
    avoid->set_active (pp->wavelet.avoid);
    avoidConn.block (false);
    tmrConn.block (true);
    tmr->set_active (pp->wavelet.tmr);
    tmrConn.block (false);
    medianConn.block (true);
    median->set_active (pp->wavelet.median);
    medianConn.block (false);
    medianlevConn.block (true);
    medianlev->set_active (pp->wavelet.medianlev);
    medianlevConn.block (false);
    linkedgConn.block (true);
    linkedg->set_active (pp->wavelet.linkedg);
    linkedgConn.block (false);
    cbenabConn.block (true);
    cbenab->set_active (pp->wavelet.cbenab);
    cbenabConn.block (false);

    lipstConn.block (true);
    lipst->set_active (pp->wavelet.lipst);
    lipstConn.block (false);
    //edgreinfConn.block (true);
    //edgreinf->set_active (pp->wavelet.edgreinf);
    //edgreinfConn.block (false);
    //lastedgreinf = pp->wavelet.edgreinf;
    lastmedian = pp->wavelet.median;
    lastmedianlev = pp->wavelet.medianlev;
    lastlinkedg = pp->wavelet.linkedg;
    lastcbenab = pp->wavelet.cbenab;
    lastlipst = pp->wavelet.lipst;
    lastavoid = pp->wavelet.avoid;
    lasttmr = pp->wavelet.tmr;
    rescon->setValue (pp->wavelet.rescon);
    resconH->setValue (pp->wavelet.resconH);
    reschro->setValue (pp->wavelet.reschro);
    tmrs->setValue (pp->wavelet.tmrs);
    gamma->setValue (pp->wavelet.gamma);
    sup->setValue (pp->wavelet.sup);
    sky->setValue (pp->wavelet.sky);
    thres->setValue (pp->wavelet.thres);
    chroma->setValue (pp->wavelet.chroma);
    chro->setValue (pp->wavelet.chro);
    contrast->setValue (pp->wavelet.contrast);
    edgrad->setValue (pp->wavelet.edgrad);
    edgval->setValue (pp->wavelet.edgval);
    edgthresh->setValue (pp->wavelet.edgthresh);
    thr->setValue (pp->wavelet.thr);
    thrH->setValue (pp->wavelet.thrH);
    skinprotect->setValue(pp->wavelet.skinprotect);
    hueskin->setValue<int>(pp->wavelet.hueskin);
    hueskin2->setValue<int>(pp->wavelet.hueskin2);
    threshold->setValue(pp->wavelet.threshold);
    threshold2->setValue(pp->wavelet.threshold2);
    edgedetect->setValue(pp->wavelet.edgedetect);
    edgedetectthr->setValue(pp->wavelet.edgedetectthr);
    edgedetectthr2->setValue(pp->wavelet.edgedetectthr2);
    edgesensi->setValue(pp->wavelet.edgesensi);
    edgeampli->setValue(pp->wavelet.edgeampli);
    hllev->setValue<int>(pp->wavelet.hllev);
    bllev->setValue<int>(pp->wavelet.bllev);
    pastlev->setValue<int>(pp->wavelet.pastlev);
    satlev->setValue<int>(pp->wavelet.satlev);
    edgcont->setValue<int>(pp->wavelet.edgcont);

    greenlow->setValue (pp->wavelet.greenlow);
    bluelow->setValue (pp->wavelet.bluelow);
    greenmed->setValue (pp->wavelet.greenmed);
    bluemed->setValue (pp->wavelet.bluemed);
    greenhigh->setValue (pp->wavelet.greenhigh);
    bluehigh->setValue (pp->wavelet.bluehigh);

    level0noise->setValue<double>(pp->wavelet.level0noise);
    level1noise->setValue<double>(pp->wavelet.level1noise);
    level2noise->setValue<double>(pp->wavelet.level2noise);
    level3noise->setValue<double>(pp->wavelet.level3noise);
    strength->setValue(pp->wavelet.strength);
    balance->setValue(pp->wavelet.balance);
    iter->setValue(pp->wavelet.iter);

    for (int i = 0; i < 9; i++) {
        correction[i]->setValue(pp->wavelet.c[i]);
    }

    for (int i = 0; i < 9; i++) {
        correctionch[i]->setValue(pp->wavelet.ch[i]);
    }

    /*****************************************************************************************************
     *
     *           Set the inconsistent state (for combobox, select the "GENERAL_UNCHANGED" entry)
     *
     *****************************************************************************************************/

    if (pedited) {
        if (!pedited->wavelet.Lmethod) {
            Lmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.CLmethod) {
            CLmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.Backmethod) {
            Backmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.Tilesmethod) {
            Tilesmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.daubcoeffmethod) {
            daubcoeffmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.Dirmethod) {
            Dirmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.CHmethod) {
            CHmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.CHSLmethod) {
            CHSLmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.EDmethod) {
            EDmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.NPmethod) {
            NPmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.BAmethod) {
            BAmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.TMmethod) {
            TMmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.HSmethod) {
            HSmethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->wavelet.Medgreinf) {
            Medgreinf->set_active_text(M("GENERAL_UNCHANGED"));
        }

        set_inconsistent (multiImage && !pedited->wavelet.enabled);
        ccshape->setUnChanged  (!pedited->wavelet.ccwcurve);
        expcontrast->set_inconsistent   (!pedited->wavelet.expcontrast);
        expchroma->set_inconsistent   (!pedited->wavelet.expchroma);
        expedge->set_inconsistent   (!pedited->wavelet.expedge);
        expresid->set_inconsistent   (!pedited->wavelet.expresid);
        expfinal->set_inconsistent   (!pedited->wavelet.expfinal);
        exptoning->set_inconsistent   (!pedited->wavelet.exptoning);
        expnoise->set_inconsistent   (!pedited->wavelet.expnoise);
        opacityShapeRG->setCurve (pp->wavelet.opacityCurveRG);
        opacityShapeBY->setCurve (pp->wavelet.opacityCurveBY);
        opacityShape->setCurve (pp->wavelet.opacityCurveW);
        opacityShapeWL->setCurve (pp->wavelet.opacityCurveWL);
        hhshape->setUnChanged  (!pedited->wavelet.hhcurve);
        Chshape->setUnChanged  (!pedited->wavelet.Chcurve);
        clshape->setUnChanged  (!pedited->wavelet.wavclCurve);
        avoid->set_inconsistent (!pedited->wavelet.avoid);
        tmr->set_inconsistent (!pedited->wavelet.tmr);
        edgthresh->setEditedState (pedited->wavelet.edgthresh ? Edited : UnEdited);
        rescon->setEditedState (pedited->wavelet.rescon ? Edited : UnEdited);
        resconH->setEditedState (pedited->wavelet.resconH ? Edited : UnEdited);
        reschro->setEditedState (pedited->wavelet.reschro ? Edited : UnEdited);
        tmrs->setEditedState (pedited->wavelet.tmrs ? Edited : UnEdited);
        gamma->setEditedState (pedited->wavelet.gamma ? Edited : UnEdited);
        sup->setEditedState (pedited->wavelet.sup ? Edited : UnEdited);
        sky->setEditedState (pedited->wavelet.sky ? Edited : UnEdited);
        thres->setEditedState (pedited->wavelet.thres ? Edited : UnEdited);
        balance->setEditedState (pedited->wavelet.balance ? Edited : UnEdited);
        iter->setEditedState (pedited->wavelet.iter ? Edited : UnEdited);
        threshold->setEditedState (pedited->wavelet.threshold ? Edited : UnEdited);
        threshold2->setEditedState (pedited->wavelet.threshold2 ? Edited : UnEdited);
        edgedetect->setEditedState (pedited->wavelet.edgedetect ? Edited : UnEdited);
        edgedetectthr->setEditedState (pedited->wavelet.edgedetectthr ? Edited : UnEdited);
        edgedetectthr2->setEditedState (pedited->wavelet.edgedetectthr2 ? Edited : UnEdited);
        edgesensi->setEditedState (pedited->wavelet.edgesensi ? Edited : UnEdited);
        edgeampli->setEditedState (pedited->wavelet.edgeampli ? Edited : UnEdited);
        chroma->setEditedState (pedited->wavelet.chroma ? Edited : UnEdited);
        chro->setEditedState (pedited->wavelet.chro ? Edited : UnEdited);

        greenlow->setEditedState (pedited->wavelet.greenlow ? Edited : UnEdited);
        bluelow->setEditedState (pedited->wavelet.bluelow ? Edited : UnEdited);
        greenmed->setEditedState (pedited->wavelet.greenmed ? Edited : UnEdited);
        bluemed->setEditedState (pedited->wavelet.bluemed ? Edited : UnEdited);
        greenhigh->setEditedState (pedited->wavelet.greenhigh ? Edited : UnEdited);
        bluehigh->setEditedState (pedited->wavelet.bluehigh ? Edited : UnEdited);

        median->set_inconsistent (!pedited->wavelet.median);
        medianlev->set_inconsistent (!pedited->wavelet.medianlev);
        linkedg->set_inconsistent (!pedited->wavelet.linkedg);
//      edgreinf->set_inconsistent (!pedited->wavelet.edgreinf);
        cbenab->set_inconsistent (!pedited->wavelet.cbenab);
        lipst->set_inconsistent (!pedited->wavelet.lipst);
        contrast->setEditedState (pedited->wavelet.contrast ? Edited : UnEdited);
        edgrad->setEditedState (pedited->wavelet.edgrad ? Edited : UnEdited);
        edgval->setEditedState (pedited->wavelet.edgval ? Edited : UnEdited);
        thr->setEditedState (pedited->wavelet.thr ? Edited : UnEdited);
        thrH->setEditedState (pedited->wavelet.thrH ? Edited : UnEdited);
        skinprotect->setEditedState (pedited->wavelet.skinprotect ? Edited : UnEdited);
        hueskin->setEditedState     (pedited->wavelet.hueskin ? Edited : UnEdited);
        hueskin2->setEditedState    (pedited->wavelet.hueskin2 ? Edited : UnEdited);
        hllev->setEditedState   (pedited->wavelet.hllev ? Edited : UnEdited);
        bllev->setEditedState   (pedited->wavelet.bllev ? Edited : UnEdited);
        pastlev->setEditedState     (pedited->wavelet.pastlev ? Edited : UnEdited);
        satlev->setEditedState  (pedited->wavelet.satlev ? Edited : UnEdited);
        strength->setEditedState(pedited->wavelet.strength ? Edited : UnEdited);
        edgcont->setEditedState     (pedited->wavelet.edgcont ? Edited : UnEdited);
        level0noise->setEditedState     (pedited->wavelet.level0noise ? Edited : UnEdited);
        level1noise->setEditedState     (pedited->wavelet.level1noise ? Edited : UnEdited);
        level2noise->setEditedState     (pedited->wavelet.level2noise ? Edited : UnEdited);
        level3noise->setEditedState     (pedited->wavelet.level3noise ? Edited : UnEdited);

        for(int i = 0; i < 9; i++) {
            correction[i]->setEditedState (pedited->wavelet.c[i] ? Edited : UnEdited);
        }

        for(int i = 0; i < 9; i++) {
            correctionch[i]->setEditedState (pedited->wavelet.ch[i] ? Edited : UnEdited);
        }
    }

    /*****************************************************************************************************
     *
     *        Update the GUI, all at once if not in Batch editing (in this case, display EVERYTHING)
     *
     *****************************************************************************************************/

    if (!batchMode) {
        int y;
        y = thres->getValue();
        int z;

        for(z = y; z < 9; z++) {
            correction[z]->hide();
        }

        for(z = 0; z < y; z++) {
            correction[z]->show();
        }

        if (pp->wavelet.CHSLmethod == "SL") {
            for(z = y; z < 9; z++) {
                correctionch[z]->hide();
            }

            for(z = 0; z < y; z++) {
                correctionch[z]->show();
            }
        }

        //adjusterUpdateUI(tmrs);
        HSmethodUpdateUI();
        CHmethodUpdateUI();
        //MedgreinfUpdateUI();
        //CHSLmethodUpdateUI();
        EDmethodUpdateUI();
        NPmethodUpdateUI();
        BAmethodUpdateUI();
        TMmethodUpdateUI();
        //BackmethodUpdateUI();
        CLmethodUpdateUI();
        lipstUpdateUI ();
        //TilesmethodUpdateUI();
        //daubcoeffmethodUpdateUI();
        //DirmethodUpdateUI();
        //LmethodUpdateUI();
        enabledUpdateUI ();
        medianlevUpdateUI ();
        cbenabUpdateUI ();

        if(z == 9) {
            sup->show();
        } else {
            sup->hide();
        }
    }

    /*****************************************************************************************************
     *
     *                                        Enable listeners and signals
     *
     *****************************************************************************************************/

    Lmethodconn.block(false);
    CLmethodconn.block(false);
    Backmethodconn.block(false);
    Tilesmethodconn.block(false);
    daubcoeffmethodconn.block(false);
    CHmethodconn.block(false);
    CHSLmethodconn.block(false);
    EDmethodconn.block(false);
    NPmethodconn.block(false);
    BAmethodconn.block(false);
    TMmethodconn.block(false);
    HSmethodconn.block(false);
    Dirmethodconn.block(false);
    MedgreinfConn.block(false);
    enableChromaConn.block(false);
    enableContrastConn.block(false);
    enableEdgeConn.block(false);
    enableFinalConn.block(false);
    enableNoiseConn.block(false);
    enableResidConn.block(false);
    enableToningConn.block(false);

    enableListener ();
}

void Wavelet::setEditProvider  (EditDataProvider *provider)
{
    ccshape->setEditProvider(provider);
    opacityShapeRG->setEditProvider(provider);
    opacityShapeBY->setEditProvider(provider);
    opacityShape->setEditProvider(provider);
    opacityShapeWL->setEditProvider(provider);
    hhshape->setEditProvider(provider);
    Chshape->setEditProvider(provider);
    clshape->setEditProvider(provider);
}

void Wavelet::autoOpenCurve ()
{
    ccshape->openIfNonlinear();
    //opacityShapeRG->openIfNonlinear();
    //opacityShapeBY->openIfNonlinear();
}

void Wavelet::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->wavelet.enabled        = getEnabled();
    pp->wavelet.avoid          = avoid->get_active ();
    pp->wavelet.tmr            = tmr->get_active ();
    pp->wavelet.rescon         = rescon->getValue();
    pp->wavelet.resconH        = resconH->getValue();
    pp->wavelet.reschro        = reschro->getValue();
    pp->wavelet.tmrs           = tmrs->getValue();
    pp->wavelet.gamma          = gamma->getValue();
    pp->wavelet.sup            = sup->getValue();
    pp->wavelet.sky            = sky->getValue();
    pp->wavelet.thres          = thres->getValue();
    pp->wavelet.chroma         = chroma->getValue();
    pp->wavelet.chro           = chro->getValue();
    pp->wavelet.median         = median->get_active ();
    pp->wavelet.medianlev      = medianlev->get_active ();
    pp->wavelet.linkedg        = linkedg->get_active ();
//  pp->wavelet.edgreinf       = edgreinf->get_active ();
    pp->wavelet.cbenab         = cbenab->get_active ();
    pp->wavelet.lipst          = lipst->get_active ();
    pp->wavelet.contrast       = contrast->getValue();
    pp->wavelet.edgrad         = edgrad->getValue();
    pp->wavelet.edgval         = edgval->getValue();
    pp->wavelet.edgthresh      = edgthresh->getValue();
    pp->wavelet.thr            = thr->getValue();
    pp->wavelet.thrH           = thrH->getValue();
    pp->wavelet.hueskin        = hueskin->getValue<int> ();
    pp->wavelet.hueskin2       = hueskin2->getValue<int> ();
    pp->wavelet.skinprotect    = skinprotect->getValue();
    pp->wavelet.threshold      = threshold->getValue();
    pp->wavelet.threshold2     = threshold2->getValue();
    pp->wavelet.edgedetect     = edgedetect->getValue();
    pp->wavelet.edgedetectthr  = edgedetectthr->getValue();
    pp->wavelet.edgedetectthr2 = edgedetectthr2->getValue();
    pp->wavelet.edgesensi     = edgesensi->getValue();
    pp->wavelet.edgeampli     = edgeampli->getValue();
    pp->wavelet.hllev          = hllev->getValue<int> ();
    pp->wavelet.bllev          = bllev->getValue<int> ();
    pp->wavelet.edgcont        = edgcont->getValue<int> ();
    pp->wavelet.level0noise    = level0noise->getValue<double> ();
    pp->wavelet.level1noise    = level1noise->getValue<double> ();
    pp->wavelet.level2noise    = level2noise->getValue<double> ();
    pp->wavelet.level3noise    = level3noise->getValue<double> ();
    pp->wavelet.ccwcurve       = ccshape->getCurve ();
    pp->wavelet.opacityCurveRG = opacityShapeRG->getCurve ();
    pp->wavelet.opacityCurveBY = opacityShapeBY->getCurve ();
    pp->wavelet.opacityCurveW  = opacityShape->getCurve ();
    pp->wavelet.opacityCurveWL = opacityShapeWL->getCurve ();
    pp->wavelet.hhcurve        = hhshape->getCurve ();
    pp->wavelet.Chcurve        = Chshape->getCurve ();
    pp->wavelet.pastlev        = pastlev->getValue<int> ();
    pp->wavelet.satlev         = satlev->getValue<int> ();
    pp->wavelet.strength       = (int) strength->getValue();
    pp->wavelet.balance        = (int) balance->getValue();

    pp->wavelet.greenlow       = greenlow->getValue ();
    pp->wavelet.bluelow        = bluelow->getValue ();
    pp->wavelet.greenmed       = greenmed->getValue ();
    pp->wavelet.bluemed        = bluemed->getValue ();
    pp->wavelet.greenhigh      = greenhigh->getValue ();
    pp->wavelet.bluehigh       = bluehigh->getValue ();
    pp->wavelet.expcontrast    = expcontrast->getEnabled();
    pp->wavelet.expchroma      = expchroma->getEnabled();
    pp->wavelet.expedge        = expedge->getEnabled();
    pp->wavelet.expresid       = expresid->getEnabled();
    pp->wavelet.expfinal       = expfinal->getEnabled();
    pp->wavelet.exptoning      = exptoning->getEnabled();
    pp->wavelet.expnoise       = expnoise->getEnabled();

    pp->wavelet.iter = (int) iter->getValue();
    pp->wavelet.wavclCurve = clshape->getCurve ();

    for (int i = 0; i < 9; i++) {
        pp->wavelet.c[i] = (int) correction[i]->getValue();
    }

    for (int i = 0; i < 9; i++) {
        pp->wavelet.ch[i] = (int) correctionch[i]->getValue();
    }

    if (pedited) {
        pedited->wavelet.enabled         = !get_inconsistent();
        pedited->wavelet.avoid           = !avoid->get_inconsistent();
        pedited->wavelet.tmr             = !tmr->get_inconsistent();
        pedited->wavelet.median          = !median->get_inconsistent();
        pedited->wavelet.medianlev       = !medianlev->get_inconsistent();
        pedited->wavelet.linkedg         = !linkedg->get_inconsistent();
        pedited->wavelet.cbenab          = !cbenab->get_inconsistent();
        pedited->wavelet.lipst           = !lipst->get_inconsistent();
        pedited->wavelet.Medgreinf       =  Medgreinf->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.Lmethod         = Lmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.CLmethod        = CLmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.Backmethod      = Backmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.Tilesmethod     = Tilesmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.daubcoeffmethod = daubcoeffmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.CHmethod        = CHmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.CHSLmethod      = CHSLmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.EDmethod        = EDmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.NPmethod        = NPmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.BAmethod        = BAmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.TMmethod        = TMmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.HSmethod        = HSmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.Dirmethod       = Dirmethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->wavelet.edgthresh       = edgthresh->getEditedState();
        pedited->wavelet.rescon          = rescon->getEditedState();
        pedited->wavelet.resconH         = resconH->getEditedState();
        pedited->wavelet.reschro         = reschro->getEditedState();
        pedited->wavelet.tmrs            = tmrs->getEditedState();
        pedited->wavelet.gamma           = gamma->getEditedState();
        pedited->wavelet.sup             = sup->getEditedState();
        pedited->wavelet.sky             = sky->getEditedState();
        pedited->wavelet.thres           = thres->getEditedState();
        pedited->wavelet.threshold       = threshold->getEditedState();
        pedited->wavelet.threshold2      = threshold2->getEditedState();
        pedited->wavelet.edgedetect      = edgedetect->getEditedState();
        pedited->wavelet.edgedetectthr   = edgedetectthr->getEditedState();
        pedited->wavelet.edgedetectthr2  = edgedetectthr2->getEditedState();
        pedited->wavelet.edgesensi       = edgesensi->getEditedState();
        pedited->wavelet.edgeampli       = edgeampli->getEditedState();
        pedited->wavelet.chroma          = chroma->getEditedState();
        pedited->wavelet.chro            = chro->getEditedState();
        pedited->wavelet.contrast        = contrast->getEditedState();
        pedited->wavelet.edgrad          = edgrad->getEditedState();
        pedited->wavelet.edgval          = edgval->getEditedState();
        pedited->wavelet.thr             = thr->getEditedState();
        pedited->wavelet.thrH            = thrH->getEditedState();
        pedited->wavelet.hueskin         = hueskin->getEditedState ();
        pedited->wavelet.hueskin2        = hueskin2->getEditedState ();
        pedited->wavelet.skinprotect     = skinprotect->getEditedState();
        pedited->wavelet.hllev           = hllev->getEditedState ();
        pedited->wavelet.ccwcurve        = !ccshape->isUnChanged ();
        pedited->wavelet.edgcont         = edgcont->getEditedState ();
        pedited->wavelet.level0noise     = level0noise->getEditedState ();
        pedited->wavelet.level1noise     = level1noise->getEditedState ();
        pedited->wavelet.level2noise     = level2noise->getEditedState ();
        pedited->wavelet.level3noise     = level3noise->getEditedState ();
        pedited->wavelet.opacityCurveRG  = !opacityShapeRG->isUnChanged ();
        pedited->wavelet.opacityCurveBY  = !opacityShapeBY->isUnChanged ();
        pedited->wavelet.opacityCurveW   = !opacityShape->isUnChanged ();
        pedited->wavelet.opacityCurveWL  = !opacityShapeWL->isUnChanged ();
        pedited->wavelet.hhcurve         = !hhshape->isUnChanged ();
        pedited->wavelet.Chcurve         = !Chshape->isUnChanged ();
        pedited->wavelet.bllev           = bllev->getEditedState ();
        pedited->wavelet.pastlev         = pastlev->getEditedState ();
        pedited->wavelet.satlev          = satlev->getEditedState ();
        pedited->wavelet.strength        = strength->getEditedState ();
        pedited->wavelet.greenlow        = greenlow->getEditedState ();
        pedited->wavelet.bluelow         = bluelow->getEditedState ();
        pedited->wavelet.greenmed        = greenmed->getEditedState ();
        pedited->wavelet.bluemed         = bluemed->getEditedState ();
        pedited->wavelet.greenhigh       = greenhigh->getEditedState ();
        pedited->wavelet.bluehigh        = bluehigh->getEditedState ();
        pedited->wavelet.balance         = balance->getEditedState ();
        pedited->wavelet.iter            = iter->getEditedState ();
        pedited->wavelet.wavclCurve      = !clshape->isUnChanged ();
        pedited->wavelet.expcontrast     = !expcontrast->get_inconsistent();
        pedited->wavelet.expchroma       = !expchroma->get_inconsistent();
        pedited->wavelet.expedge         = !expedge->get_inconsistent();
        pedited->wavelet.expresid        = !expresid->get_inconsistent();
        pedited->wavelet.expfinal        = !expfinal->get_inconsistent();
        pedited->wavelet.exptoning       = !exptoning->get_inconsistent();
        pedited->wavelet.expnoise        = !expnoise->get_inconsistent();

        for(int i = 0; i < 9; i++) {
            pedited->wavelet.c[i]        = correction[i]->getEditedState();
        }

        for(int i = 0; i < 9; i++) {
            pedited->wavelet.ch[i]       = correctionch[i]->getEditedState();
        }

    }

    if (CHmethod->get_active_row_number() == 0) {
        pp->wavelet.CHmethod = "without";
    } else if (CHmethod->get_active_row_number() == 1) {
        pp->wavelet.CHmethod = "with";
    } else if (CHmethod->get_active_row_number() == 2) {
        pp->wavelet.CHmethod = "link";
    }

    if (Medgreinf->get_active_row_number() == 0) {
        pp->wavelet.Medgreinf = "more";
    } else if (Medgreinf->get_active_row_number() == 1) {
        pp->wavelet.Medgreinf = "none";
    } else if (Medgreinf->get_active_row_number() == 2) {
        pp->wavelet.Medgreinf = "less";
    }

    if (CHSLmethod->get_active_row_number() == 0) {
        pp->wavelet.CHSLmethod = "SL";
    } else if (CHSLmethod->get_active_row_number() == 1) {
        pp->wavelet.CHSLmethod = "CU";
    }

    if (EDmethod->get_active_row_number() == 0) {
        pp->wavelet.EDmethod = "SL";
    } else if (EDmethod->get_active_row_number() == 1) {
        pp->wavelet.EDmethod = "CU";
    }

    if (NPmethod->get_active_row_number() == 0) {
        pp->wavelet.NPmethod = "none";
    } else if (NPmethod->get_active_row_number() == 1) {
        pp->wavelet.NPmethod = "low";
    } else if (NPmethod->get_active_row_number() == 2) {
        pp->wavelet.NPmethod = "high";
    }

    if (BAmethod->get_active_row_number() == 0) {
        pp->wavelet.BAmethod = "none";
    } else if (BAmethod->get_active_row_number() == 1) {
        pp->wavelet.BAmethod = "sli";
    } else if (BAmethod->get_active_row_number() == 2) {
        pp->wavelet.BAmethod = "cur";
    }

    if (TMmethod->get_active_row_number() == 0) {
        pp->wavelet.TMmethod = "cont";
    } else if (TMmethod->get_active_row_number() == 1) {
        pp->wavelet.TMmethod = "tm";
    }

    //  else if (TMmethod->get_active_row_number()==2)
    //      pp->wavelet.TMmethod = "both";

    if (HSmethod->get_active_row_number() == 0) {
        pp->wavelet.HSmethod = "without";
    } else if (HSmethod->get_active_row_number() == 1) {
        pp->wavelet.HSmethod = "with";
    }

    if (CLmethod->get_active_row_number() == 0) {
        pp->wavelet.CLmethod = "one";
    } else if (CLmethod->get_active_row_number() == 1) {
        pp->wavelet.CLmethod = "inf";
    } else if (CLmethod->get_active_row_number() == 2) {
        pp->wavelet.CLmethod = "sup";
    } else if (CLmethod->get_active_row_number() == 3) {
        pp->wavelet.CLmethod = "all";
    }

    if (Backmethod->get_active_row_number() == 0) {
        pp->wavelet.Backmethod = "black";
    } else if (Backmethod->get_active_row_number() == 1) {
        pp->wavelet.Backmethod = "grey";
    } else if (Backmethod->get_active_row_number() == 2) {
        pp->wavelet.Backmethod = "resid";
    }

    if (Tilesmethod->get_active_row_number() == 0) {
        pp->wavelet.Tilesmethod = "full";
    } else if (Tilesmethod->get_active_row_number() == 1) {
        pp->wavelet.Tilesmethod = "big";
    } else if (Tilesmethod->get_active_row_number() == 2) {
        pp->wavelet.Tilesmethod = "lit";
    }

    if (daubcoeffmethod->get_active_row_number() == 0) {
        pp->wavelet.daubcoeffmethod = "2_";
    } else if (daubcoeffmethod->get_active_row_number() == 1) {
        pp->wavelet.daubcoeffmethod = "4_";
    } else if (daubcoeffmethod->get_active_row_number() == 2) {
        pp->wavelet.daubcoeffmethod = "6_";
    } else if (daubcoeffmethod->get_active_row_number() == 3) {
        pp->wavelet.daubcoeffmethod = "10_";
    } else if (daubcoeffmethod->get_active_row_number() == 4) {
        pp->wavelet.daubcoeffmethod = "14_";
    }

    if (Dirmethod->get_active_row_number() == 0) {
        pp->wavelet.Dirmethod = "one";
    } else if (Dirmethod->get_active_row_number() == 1) {
        pp->wavelet.Dirmethod = "two";
    } else if (Dirmethod->get_active_row_number() == 2) {
        pp->wavelet.Dirmethod = "thr";
    } else if (Dirmethod->get_active_row_number() == 3) {
        pp->wavelet.Dirmethod = "all";
    }

    char lMethod[3]; // one additional char to avoid buffer overrun if someone increases number of levels > 9
    sprintf(lMethod, "%d", Lmethod->get_active_row_number() + 1);
    pp->wavelet.Lmethod = lMethod;
}

void Wavelet::curveChanged (CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == ccshape) {
            listener->panelChanged (EvWavCCCurve, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == opacityShapeRG) {
            listener->panelChanged (EvWavColor, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == opacityShapeBY) {
            listener->panelChanged (EvWavOpac, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == opacityShape) {
            listener->panelChanged (EvWavopacity, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == opacityShapeWL) {
            listener->panelChanged (EvWavopacityWL, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == hhshape) {
            listener->panelChanged (EvWavHHCurve, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == Chshape) {
            listener->panelChanged (EvWavCHCurve, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == clshape) {
            listener->panelChanged (EvWavCLCurve, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void Wavelet::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    for (int i = 0; i < 9; i++) {
        correction[i]->setDefault(defParams->wavelet.c[i]);
    }

    for (int i = 0; i < 9; i++) {
        correctionch[i]->setDefault(defParams->wavelet.ch[i]);
    }

    strength->setDefault(defParams->wavelet.strength );
    balance->setDefault(defParams->wavelet.balance );
    iter->setDefault(defParams->wavelet.iter );
    rescon->setDefault (defParams->wavelet.rescon);
    resconH->setDefault (defParams->wavelet.resconH);
    reschro->setDefault (defParams->wavelet.reschro);
    tmrs->setDefault (defParams->wavelet.tmrs);
    gamma->setDefault (defParams->wavelet.gamma);
    sup->setDefault (defParams->wavelet.sup);
    sky->setDefault (defParams->wavelet.sky);
    thres->setDefault (defParams->wavelet.thres);
    threshold->setDefault (defParams->wavelet.threshold);
    threshold2->setDefault (defParams->wavelet.threshold2);
    edgedetect->setDefault (defParams->wavelet.edgedetect);
    edgedetectthr->setDefault (defParams->wavelet.edgedetectthr);
    edgedetectthr2->setDefault (defParams->wavelet.edgedetectthr2);
    edgesensi->setDefault (defParams->wavelet.edgesensi);
    edgeampli->setDefault (defParams->wavelet.edgeampli);
    chroma->setDefault (defParams->wavelet.chroma);
    chro->setDefault (defParams->wavelet.chro);
    contrast->setDefault (defParams->wavelet.contrast);
    edgrad->setDefault (defParams->wavelet.edgrad);
    edgval->setDefault (defParams->wavelet.edgval);
    edgthresh->setDefault (defParams->wavelet.edgthresh);
    thr->setDefault (defParams->wavelet.thr);
    thrH->setDefault (defParams->wavelet.thrH);
    hueskin->setDefault<int> (defParams->wavelet.hueskin);
    hueskin2->setDefault<int> (defParams->wavelet.hueskin2);
    hllev->setDefault<int> (defParams->wavelet.hllev);
    bllev->setDefault<int> (defParams->wavelet.bllev);
    pastlev->setDefault<int> (defParams->wavelet.pastlev);
    satlev->setDefault<int> (defParams->wavelet.satlev);
    edgcont->setDefault<int> (defParams->wavelet.edgcont);
    level0noise->setDefault<double> (defParams->wavelet.level0noise);
    level1noise->setDefault<double> (defParams->wavelet.level1noise);
    level2noise->setDefault<double> (defParams->wavelet.level2noise);
    level3noise->setDefault<double> (defParams->wavelet.level3noise);

    greenlow->setDefault (defParams->wavelet.greenlow);
    bluelow->setDefault (defParams->wavelet.bluelow);
    greenmed->setDefault (defParams->wavelet.greenmed);
    bluemed->setDefault (defParams->wavelet.bluemed);
    greenhigh->setDefault (defParams->wavelet.greenhigh);
    bluehigh->setDefault (defParams->wavelet.bluehigh);

    if (pedited) {
        greenlow->setDefaultEditedState (pedited->wavelet.greenlow ? Edited : UnEdited);
        bluelow->setDefaultEditedState (pedited->wavelet.bluelow ? Edited : UnEdited);
        greenmed->setDefaultEditedState (pedited->wavelet.greenmed ? Edited : UnEdited);
        bluemed->setDefaultEditedState (pedited->wavelet.bluemed ? Edited : UnEdited);
        greenhigh->setDefaultEditedState (pedited->wavelet.greenhigh ? Edited : UnEdited);
        bluehigh->setDefaultEditedState (pedited->wavelet.bluehigh ? Edited : UnEdited);

        rescon->setDefault (defParams->wavelet.rescon);
        resconH->setDefault (defParams->wavelet.resconH);
        reschro->setDefault (defParams->wavelet.reschro);
        tmrs->setDefault (defParams->wavelet.tmrs);
        gamma->setDefault (defParams->wavelet.gamma);
        sup->setDefault (defParams->wavelet.sup);
        sky->setDefaultEditedState(pedited->wavelet.sky ? Edited : UnEdited);
        thres->setDefaultEditedState(pedited->wavelet.thres ? Edited : UnEdited);
        threshold->setDefaultEditedState(pedited->wavelet.threshold ? Edited : UnEdited);
        threshold2->setDefaultEditedState(pedited->wavelet.threshold2 ? Edited : UnEdited);
        edgedetect->setDefaultEditedState(pedited->wavelet.edgedetect ? Edited : UnEdited);
        edgedetectthr->setDefaultEditedState(pedited->wavelet.edgedetectthr ? Edited : UnEdited);
        edgedetectthr2->setDefaultEditedState(pedited->wavelet.edgedetectthr2 ? Edited : UnEdited);
        edgesensi->setDefaultEditedState(pedited->wavelet.edgesensi ? Edited : UnEdited);
        edgeampli->setDefaultEditedState(pedited->wavelet.edgeampli ? Edited : UnEdited);
        chroma->setDefaultEditedState(pedited->wavelet.chroma ? Edited : UnEdited);
        chro->setDefaultEditedState(pedited->wavelet.chro ? Edited : UnEdited);
        contrast->setDefaultEditedState(pedited->wavelet.contrast ? Edited : UnEdited);
        edgrad->setDefaultEditedState(pedited->wavelet.edgrad ? Edited : UnEdited);
        edgval->setDefaultEditedState(pedited->wavelet.edgval ? Edited : UnEdited);
        edgthresh->setDefault (defParams->wavelet.edgthresh);
        thr->setDefaultEditedState(pedited->wavelet.thr ? Edited : UnEdited);
        thrH->setDefaultEditedState(pedited->wavelet.thrH ? Edited : UnEdited);
        skinprotect->setDefaultEditedState(pedited->wavelet.skinprotect ? Edited : UnEdited);
        hueskin->setDefaultEditedState(pedited->wavelet.hueskin ? Edited : UnEdited);
        hueskin2->setDefaultEditedState(pedited->wavelet.hueskin2 ? Edited : UnEdited);
        hllev->setDefaultEditedState(pedited->wavelet.hllev ? Edited : UnEdited);
        bllev->setDefaultEditedState(pedited->wavelet.bllev ? Edited : UnEdited);
        pastlev->setDefaultEditedState(pedited->wavelet.pastlev ? Edited : UnEdited);
        satlev->setDefaultEditedState(pedited->wavelet.satlev ? Edited : UnEdited);
        edgcont->setDefaultEditedState(pedited->wavelet.edgcont ? Edited : UnEdited);
        strength->setDefaultEditedState(pedited->wavelet.strength ? Edited : UnEdited);
        balance->setDefaultEditedState(pedited->wavelet.balance ? Edited : UnEdited);
        iter->setDefaultEditedState(pedited->wavelet.iter ? Edited : UnEdited);
        level0noise->setDefaultEditedState(pedited->wavelet.level0noise ? Edited : UnEdited);
        level1noise->setDefaultEditedState(pedited->wavelet.level1noise ? Edited : UnEdited);
        level2noise->setDefaultEditedState(pedited->wavelet.level2noise ? Edited : UnEdited);
        level3noise->setDefaultEditedState(pedited->wavelet.level3noise ? Edited : UnEdited);

        for (int i = 0; i < 9; i++) {
            correction[i]->setDefaultEditedState(pedited->wavelet.c[i] ? Edited : UnEdited);
        }

        for (int i = 0; i < 9; i++) {
            correctionch[i]->setDefaultEditedState(pedited->wavelet.ch[i] ? Edited : UnEdited);
        }
    } else {
        rescon->setDefaultEditedState(Irrelevant);
        resconH->setDefaultEditedState(Irrelevant);
        reschro->setDefaultEditedState(Irrelevant);
        tmrs->setDefaultEditedState(Irrelevant);
        gamma->setDefaultEditedState(Irrelevant);
        sup->setDefaultEditedState(Irrelevant);
        sky->setDefaultEditedState(Irrelevant);
        thres->setDefaultEditedState(Irrelevant);
        threshold->setDefaultEditedState(Irrelevant);
        threshold2->setDefaultEditedState(Irrelevant);
        edgedetect->setDefaultEditedState(Irrelevant);
        edgedetectthr->setDefaultEditedState(Irrelevant);
        edgedetectthr2->setDefaultEditedState(Irrelevant);
        edgesensi->setDefaultEditedState(Irrelevant);
        edgeampli->setDefaultEditedState(Irrelevant);
        chroma->setDefaultEditedState(Irrelevant);
        chro->setDefaultEditedState(Irrelevant);
        contrast->setDefaultEditedState(Irrelevant);
        edgrad->setDefaultEditedState(Irrelevant);
        edgval->setDefaultEditedState(Irrelevant);
        edgthresh->setDefaultEditedState(Irrelevant);
        thr->setDefaultEditedState(Irrelevant);
        thrH->setDefaultEditedState(Irrelevant);
        skinprotect->setDefaultEditedState(Irrelevant);
        hueskin->setDefaultEditedState (Irrelevant);
        hueskin2->setDefaultEditedState (Irrelevant);
        hllev->setDefaultEditedState (Irrelevant);
        bllev->setDefaultEditedState (Irrelevant);
        edgcont->setDefaultEditedState (Irrelevant);
        level0noise->setDefaultEditedState (Irrelevant);
        level1noise->setDefaultEditedState (Irrelevant);
        level2noise->setDefaultEditedState (Irrelevant);
        level3noise->setDefaultEditedState (Irrelevant);
        pastlev->setDefaultEditedState (Irrelevant);
        satlev->setDefaultEditedState (Irrelevant);
        strength->setDefaultEditedState (Irrelevant);
        balance->setDefaultEditedState (Irrelevant);
        iter->setDefaultEditedState (Irrelevant);

        for (int i = 0; i < 9; i++) {
            correction[i]->setDefaultEditedState(Irrelevant);
        }

        for (int i = 0; i < 9; i++) {
            correctionch[i]->setDefaultEditedState(Irrelevant);
        }
    }
}
void Wavelet::adjusterChanged (ThresholdAdjuster* a, double newBottom, double newTop)
{
    if (listener && (multiImage || getEnabled()) ) {
        if(a == level0noise) {
            listener->panelChanged (EvWavlev0nois,
                                    Glib::ustring::compose(Glib::ustring(M("TP_WAVELET_NOIS") + ": %1" + "\n" + M("TP_WAVELET_STREN") + ": %2"), int(newTop), int(newBottom)));
        } else if(a == level1noise) {
            listener->panelChanged (EvWavlev1nois,
                                    Glib::ustring::compose(Glib::ustring(M("TP_WAVELET_NOIS") + ": %1" + "\n" + M("TP_WAVELET_STREN") + ": %2"), int(newTop), int(newBottom)));
        } else if(a == level2noise) {
            listener->panelChanged (EvWavlev2nois,
                                    Glib::ustring::compose(Glib::ustring(M("TP_WAVELET_NOIS") + ": %1" + "\n" + M("TP_WAVELET_STREN") + ": %2"), int(newTop), int(newBottom)));
        } else if(a == level3noise) {
            listener->panelChanged (EvWavlev3nois,
                                    Glib::ustring::compose(Glib::ustring(M("TP_WAVELET_NOIS") + ": %1" + "\n" + M("TP_WAVELET_STREN") + ": %2"), int(newTop), int(newBottom)));
        }

    }
}


void Wavelet::adjusterChanged2 (ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
    if (listener && (multiImage || getEnabled()) ) {
        if(a == hueskin) {
            listener->panelChanged (EvWavHueskin, hueskin->getHistoryString());
        } else if(a == hueskin2) {
            listener->panelChanged (EvWavHueskin2, hueskin2->getHistoryString());
        } else if(a == hllev) {
            listener->panelChanged (EvWavlhl, hllev->getHistoryString());
        } else if(a == bllev) {
            listener->panelChanged (EvWavlbl, bllev->getHistoryString());
        } else if(a == pastlev) {
            listener->panelChanged (EvWavpast, pastlev->getHistoryString());
        } else if(a == satlev) {
            listener->panelChanged (EvWavsat, satlev->getHistoryString());
        } else if(a == edgcont) {
            listener->panelChanged (EvWavedgcont, edgcont->getHistoryString());
        }
    }
}

void Wavelet::HSmethodUpdateUI()
{
    if (!batchMode) {
        if(HSmethod->get_active_row_number() == 0) { //without
            hllev->hide();
            bllev->hide();
            threshold->hide();
            threshold2->hide();
        } else { //with
            hllev->show();
            bllev->show();
            threshold->show();
            threshold2->show();
        }
    }
}

void Wavelet::HSmethodChanged()
{
    HSmethodUpdateUI();

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavHSmet, HSmethod->get_active_text ());
    }
}

void Wavelet::CHmethodUpdateUI()
{
    if (!batchMode) {
        if(CHmethod->get_active_row_number() == 0) {
            CHSLmethod->show();
            pastlev->hide();
            satlev->hide();
            chroma->hide();
            chro->hide();
            labmC->show();
            separatorNeutral->hide();
            neutralchButton->show();
            int y = thres->getValue();
            int z;

            for(z = y; z < 9; z++) {
                correctionch[z]->hide();
            }

            for(z = 0; z < y; z++) {
                correctionch[z]->show();
            }
        } else if(CHmethod->get_active_row_number() == 1) {
            CHSLmethod->show();
            pastlev->show();
            satlev->show();
            chroma->show();
            chro->hide();
            labmC->show();
            separatorNeutral->show();
            neutralchButton->show();
            int y = thres->getValue();
            int z;

            for(z = y; z < 9; z++) {
                correctionch[z]->hide();
            }

            for(z = 0; z < y; z++) {
                correctionch[z]->show();
            }
        } else {
            chro->show();
            pastlev->hide();
            satlev->hide();
            chroma->hide();
            CHSLmethod->hide();
            labmC->hide();
            separatorNeutral->hide();
            neutralchButton->hide();

            for (int i = 0; i < 9; i++) {
                correctionch[i]->hide();
            }
        }
    }
}

void Wavelet::CHmethodChanged()
{
    CHmethodUpdateUI();

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavCHmet, CHmethod->get_active_text ());
    }
}

/*
void Wavelet::CHSLmethodChangedUI() {
    if (!batchMode) {
        if(CHSLmethod->get_active_row_number()==0 && CHmethod->get_active_row_number() != 2) {//SL
            //CLVcurveEditorG->hide();
            neutralchButton->show();
                int y=thres->getValue();
                int z;
                for(z=y;z<9;z++) correctionch[z]->hide();
                for(z=0;z<y;z++) correctionch[z]->show();
        }
        else if(CHSLmethod->get_active_row_number()==1) {//CU
            //CLVcurveEditorG->show();
            neutralchButton->hide();
            for (int i = 0; i < 9; i++) {
                correctionch[i]->hide();
            }
        }
    }
}
*/

void Wavelet::CHSLmethodChanged()
{
    /*
        CHSLmethodChangedUI();
        if (listener && (multiImage||enabled->get_active()) ) {
            listener->panelChanged (EvWavCHSLmet, CHSLmethod->get_active_text ());
        }
    */
}

void Wavelet::EDmethodUpdateUI()
{
    if (!batchMode) {
        if(EDmethod->get_active_row_number() == 0 ) { //SL
            CCWcurveEditorG->hide();
            edgcont->show();
        } else if(EDmethod->get_active_row_number() == 1) { //CU
            CCWcurveEditorG->show();
            edgcont->hide();
        }
    }
}
void Wavelet::EDmethodChanged()
{
    EDmethodUpdateUI();

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavEDmet, EDmethod->get_active_text ());
    }
}

void Wavelet::NPmethodUpdateUI()
{
    if (!batchMode) {

    }
}
void Wavelet::NPmethodChanged()
{
    NPmethodUpdateUI();

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavNPmet, NPmethod->get_active_text ());
    }
}



void Wavelet::BAmethodUpdateUI()
{
    if (!batchMode) {
        if(BAmethod->get_active_row_number() == 0 ) { //none
            balance->hide();
            opacityCurveEditorW->hide();
            iter->hide();

        } else if(BAmethod->get_active_row_number() == 1) { //sli
            opacityCurveEditorW->hide();
            balance->show();
            iter->show();
        } else if(BAmethod->get_active_row_number() == 2) { //CU
            opacityCurveEditorW->show();
            balance->hide();
            iter->show();
        }
    }
}
void Wavelet::BAmethodChanged()
{
    BAmethodUpdateUI();

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavBAmet, BAmethod->get_active_text ());
    }
}


void Wavelet::TMmethodUpdateUI()
{
    /*
    if (!batchMode) {
        if(TMmethod->get_active_row_number()==0 ) {//
            tmr->hide();
        }
        else {
            tmr->show();
        }
    }
    */
}

void Wavelet::TMmethodChanged()
{
    TMmethodUpdateUI();

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavTMmet, TMmethod->get_active_text ());
    }
}

/*
void Wavelet::BackmethodUpdateUI() {
    if (!batchMode) {
    }
}
*/

void Wavelet::BackmethodChanged()
{
    //BackmethodUpdateUI();
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavBackmet, Backmethod->get_active_text ());
    }
}

void Wavelet::CLmethodUpdateUI()
{
    if (!batchMode) {
        if (CLmethod->get_active_row_number() == 0) {
            CLmethod->set_active (0);
            Lmethod->set_sensitive(true);
            Dirmethod->set_sensitive(true);
        } else if (CLmethod->get_active_row_number() == 1) {
            CLmethod->set_active (1);
            Lmethod->set_sensitive(true);
            Dirmethod->set_sensitive(true);
        } else if (CLmethod->get_active_row_number() == 2) {
            CLmethod->set_active (2);
            Lmethod->set_sensitive(true);
            Dirmethod->set_sensitive(true);
        } else if (CLmethod->get_active_row_number() == 3) {
            CLmethod->set_active (3);
            Lmethod->set_sensitive(false);
            Dirmethod->set_sensitive(false);
        }
    }
}

void Wavelet::CLmethodChanged()
{
    CLmethodUpdateUI();

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavCLmet, CLmethod->get_active_text ());
    }
}

/*
void Wavelet::TilesmethodUpdateUI() {
    if (!batchMode) {
    }
}
*/

void Wavelet::TilesmethodChanged()
{
    //TilesmethodUpdateUI();
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavTilesmet, Tilesmethod->get_active_text ());
    }
}

/*
void Wavelet::daubcoeffmethodUpdateUI() {
    if (!batchMode) {
    }}
*/
void Wavelet::daubcoeffmethodChanged()
{
    //daubcoeffmethodUpdateUI();
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavdaubcoeffmet, daubcoeffmethod->get_active_text ());
    }
}

/*
void Wavelet::MedgreinfUpdateUI() {
    if (!batchMode) {
    }
}
*/
void Wavelet::MedgreinfChanged()
{
    //MedgreinfUpdateUI();
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavedgreinf, Medgreinf->get_active_text ());
    }
}

/*
void Wavelet::DirmethodUpdateUI() {
    if (!batchMode) {
    }
}
*/
void Wavelet::DirmethodChanged()
{
    //DirmethodUpdateUI();
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavDirmeto, Dirmethod->get_active_text ());
    }
}

/*
void Wavelet::LmethodUpdateUI() {
    if (!batchMode) {
    }
}
*/
void Wavelet::LmethodChanged()
{
    //LmethodUpdateUI();
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvWavLmet, Lmethod->get_active_text ());
    }
}

void Wavelet::setBatchMode (bool batchMode)
{
    Lmethod->append_text (M("GENERAL_UNCHANGED"));
    CLmethod->append_text (M("GENERAL_UNCHANGED"));
    Backmethod->append_text (M("GENERAL_UNCHANGED"));
    Tilesmethod->append_text (M("GENERAL_UNCHANGED"));
    daubcoeffmethod->append_text (M("GENERAL_UNCHANGED"));
    CHmethod->append_text (M("GENERAL_UNCHANGED"));
    Medgreinf->append_text (M("GENERAL_UNCHANGED"));
    CHSLmethod->append_text (M("GENERAL_UNCHANGED"));
    EDmethod->append_text (M("GENERAL_UNCHANGED"));
    NPmethod->append_text (M("GENERAL_UNCHANGED"));
    BAmethod->append_text (M("GENERAL_UNCHANGED"));
    TMmethod->append_text (M("GENERAL_UNCHANGED"));
    HSmethod->append_text (M("GENERAL_UNCHANGED"));
    Dirmethod->append_text (M("GENERAL_UNCHANGED"));
    CCWcurveEditorG->setBatchMode (batchMode);
    opaCurveEditorG->setBatchMode (batchMode);
    opacityCurveEditorG->setBatchMode (batchMode);
    opacityCurveEditorW->setBatchMode (batchMode);
    opacityCurveEditorWL->setBatchMode (batchMode);
    curveEditorRES->setBatchMode (batchMode);
    curveEditorGAM->setBatchMode (batchMode);
    rescon->showEditedCB ();
    resconH->showEditedCB ();
    reschro->showEditedCB ();
    tmrs->showEditedCB ();
    gamma->showEditedCB ();
    sup->showEditedCB ();
    sky->showEditedCB ();
    thres->showEditedCB ();
    threshold->showEditedCB ();
    threshold2->showEditedCB ();
    edgedetect->showEditedCB ();
    edgedetectthr->showEditedCB ();
    edgedetectthr2->showEditedCB ();
    edgesensi->showEditedCB ();
    edgeampli->showEditedCB ();
    chroma->showEditedCB ();
    chro->showEditedCB ();
    contrast->showEditedCB ();
    edgrad->showEditedCB ();
    edgval->showEditedCB ();
    edgthresh->showEditedCB ();
    thr->showEditedCB ();
    thrH->showEditedCB ();
    skinprotect->showEditedCB();
    hueskin->showEditedCB ();
    hueskin2->showEditedCB ();
    hllev->showEditedCB ();
    bllev->showEditedCB ();
    pastlev->showEditedCB ();
    satlev->showEditedCB ();
    edgcont->showEditedCB ();
    strength->showEditedCB ();
    balance->showEditedCB ();
    iter->showEditedCB ();
    level0noise->showEditedCB ();
    level1noise->showEditedCB ();
    level2noise->showEditedCB ();
    level3noise->showEditedCB ();

    ToolPanel::setBatchMode (batchMode);

    for (int i = 0; i < 9; i++) {
        correction[i]->showEditedCB();
    }

    for (int i = 0; i < 9; i++) {
        correctionch[i]->showEditedCB();
    }
}

void Wavelet::adjusterUpdateUI (Adjuster* a)
{
    /*
        if (!batchMode) {
            if (a == tmrs ) {
                float tm;
                tm=tmrs->getValue();
                if(tm==0.f) tmr->hide();
                else tmr->show();
            }
            else if (a == gamma ) {
                float tm;
                tm=tmrs->getValue();
                if(tm==0.f) tmr->hide();
                else tmr->show();
                );
            }
        }
    */
}

void Wavelet::adjusterChanged (Adjuster* a, double newval)
{
    if (listener && (multiImage || getEnabled()) ) {
        if (a == edgthresh) {
            listener->panelChanged (EvWavtiles, edgthresh->getTextValue());
        } else if (a == rescon ) {
            listener->panelChanged (EvWavrescon, rescon->getTextValue());
        } else if (a == resconH ) {
            listener->panelChanged (EvWavresconH, resconH->getTextValue());
        } else if (a == reschro ) {
            listener->panelChanged (EvWavreschro, reschro->getTextValue());
        } else if (a == tmrs ) {
            adjusterUpdateUI(a);
            listener->panelChanged (EvWavtmrs, tmrs->getTextValue());
        } else if (a == gamma ) {
            adjusterUpdateUI(a);
            listener->panelChanged (EvWavgamma, gamma->getTextValue());
        } else if (a == sky ) {
            listener->panelChanged (EvWavsky, sky->getTextValue());
        } else if (a == sup ) {
            listener->panelChanged (EvWavsup, sup->getTextValue());
        } else if (a == chroma ) {
            listener->panelChanged (EvWavchroma, chroma->getTextValue());
        } else if (a == chro ) {
            listener->panelChanged (EvWavchro, chro->getTextValue());
        } else if (a == contrast ) {
            listener->panelChanged (EvWavunif, contrast->getTextValue());
        } else if (a == thr ) {
            listener->panelChanged (EvWavthr, thr->getTextValue());
        } else if (a == thrH ) {
            listener->panelChanged (EvWavthrH, thrH->getTextValue());
        } else if (a == threshold ) {
            listener->panelChanged (EvWavThreshold, threshold->getTextValue());
        } else if (a == threshold2 ) {
            listener->panelChanged (EvWavThreshold2, threshold2->getTextValue());
        } else if (a == edgedetect ) {
            listener->panelChanged (EvWavedgedetect, edgedetect->getTextValue());
        } else if (a == edgedetectthr ) {
            listener->panelChanged (EvWavedgedetectthr, edgedetectthr->getTextValue());
        } else if (a == edgedetectthr2 ) {
            listener->panelChanged (EvWavedgedetectthr2, edgedetectthr2->getTextValue());
        } else if (a == edgesensi ) {
            listener->panelChanged (EvWavedgesensi, edgesensi->getTextValue());
        } else if (a == edgeampli ) {
            listener->panelChanged (EvWavedgeampli, edgeampli->getTextValue());
        } else if (a == edgrad ) {
            listener->panelChanged (EvWavedgrad, edgrad->getTextValue());
        } else if (a == edgval ) {
            listener->panelChanged (EvWavedgval, edgval->getTextValue());
        } else if (a == thres ) {
            int y = thres->getValue();
            int z;

            for(z = y; z < 9; z++) {
                correction[z]->hide();
            }

            for(z = 0; z < y; z++) {
                correction[z]->show();
            }

            for(z = y; z < 9; z++) {
                correctionch[z]->hide();
            }

            for(z = 0; z < y; z++) {
                correctionch[z]->show();
            }

            if(z == 9) {
                sup->show();
            } else {
                sup->hide();
            }

            listener->panelChanged (EvWavthres, thres->getTextValue());
        } else if (a == skinprotect) {
            listener->panelChanged (EvWavSkin, skinprotect->getTextValue());
        } else if (a == strength) {
            listener->panelChanged (EvWavStrength, strength->getTextValue());
        } else if (a == balance) {
            listener->panelChanged (EvWavbalance, balance->getTextValue());
        } else if (a == iter) {
            listener->panelChanged (EvWaviter, iter->getTextValue());
        } else if (a == greenhigh ) {
            listener->panelChanged (EvWavgreenhigh, greenhigh->getTextValue());
        } else if (a == bluehigh ) {
            listener->panelChanged (EvWavbluehigh, bluehigh->getTextValue());
        } else if (a == greenmed ) {
            listener->panelChanged (EvWavgreenmed, greenmed->getTextValue());
        } else if (a == bluemed ) {
            listener->panelChanged (EvWavbluemed, bluemed->getTextValue());
        } else if (a == greenlow ) {
            listener->panelChanged (EvWavgreenlow, greenlow->getTextValue());
        } else if (a == bluelow ) {
            listener->panelChanged (EvWavbluelow, bluelow->getTextValue());
        }

        if ((a == correction[0] || a == correction[1] || a == correction[2] || a == correction[3] || a == correction[4] || a == correction[5] || a == correction[6] || a == correction[7] || a == correction[8])) {
            listener->panelChanged (EvWavelet,
                                    Glib::ustring::compose("%1, %2, %3, %4, %5, %6, %7, %8, %9",
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[0]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[1]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[2]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[3]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[4]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[5]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[6]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[7]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correction[8]->getValue()))
                                   );
        } else if (a == correctionch[0] || a == correctionch[1] || a == correctionch[2] || a == correctionch[3] || a == correctionch[4] || a == correctionch[5] || a == correctionch[6] || a == correctionch[7] || a == correctionch[8] ) {
            listener->panelChanged (EvWaveletch,
                                    Glib::ustring::compose("%1, %2, %3, %4, %5, %6, %7, %8, %9",
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[0]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[1]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[2]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[3]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[4]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[5]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[6]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[7]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(0), correctionch[8]->getValue()))
                                   );
        }
    }
}

void Wavelet::enabledUpdateUI ()
{
    if (!batchMode) {
        int y = thres->getValue();
        int z;

        for(z = y; z < 9; z++) {
            correction[z]->hide();
        }

        for(z = 0; z < y; z++) {
            correction[z]->show();
        }

        if(z == 9) {
            sup->show();
        } else {
            sup->hide();
        }

//      adjusterUpdateUI(tmrs);
    }
}

void Wavelet::enabledChanged ()
{
    enabledUpdateUI();

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvWavEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvWavEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvWavEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Wavelet::medianToggled ()
{

    if (multiImage) {
        if (median->get_inconsistent()) {
            median->set_inconsistent (false);
            medianConn.block (true);
            median->set_active (false);
            medianConn.block (false);
        } else if (lastmedian) {
            median->set_inconsistent (true);
        }

        lastmedian = median->get_active ();
    }

    if (listener && (multiImage || getEnabled ())) {
        if (median->get_inconsistent()) {
            listener->panelChanged (EvWavmedian, M("GENERAL_UNCHANGED"));
        } else if (median->get_active ()) {
            listener->panelChanged (EvWavmedian, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvWavmedian, M("GENERAL_DISABLED"));
        }
    }
}

void Wavelet::medianlevUpdateUI ()
{
    if (!batchMode) {
        if (medianlev->get_active ()) {
            edgedetect->show();
            lipst->show();
            separatoredge->show();

            edgedetectthr->show();
            edgedetectthr2->show();

            //  edgesensi->show();
            //  edgeampli->show();
            //  NPmethod->show();
            //  labmNP->show();
            if (lipst->get_active ()) {
                edgesensi->show();
                edgeampli->show();
                NPmethod->show();
                labmNP->show();
            }

            else {
                edgesensi->hide();
                edgeampli->hide();
                NPmethod->hide();
                labmNP->hide();
            }
        } else {
            edgedetect->hide();
            lipst->hide();
            edgedetectthr->hide();
            edgedetectthr2->hide();
            edgesensi->hide();
            edgeampli->hide();
            separatoredge->hide();
            NPmethod->hide();
            labmNP->hide();

        }
    }
}

void Wavelet::medianlevToggled ()
{

    if (multiImage) {
        if (medianlev->get_inconsistent()) {
            medianlev->set_inconsistent (false);
            medianlevConn.block (true);
            medianlev->set_active (false);
            medianlevConn.block (false);
        } else if (lastmedianlev) {
            medianlev->set_inconsistent (true);
        }

        lastmedianlev = medianlev->get_active ();
    }

    medianlevUpdateUI();

    if (listener && (multiImage || getEnabled ())) {
        if (medianlev->get_inconsistent()) {
            listener->panelChanged (EvWavmedianlev, M("GENERAL_UNCHANGED"));
        } else if (medianlev->get_active () ) {
            listener->panelChanged (EvWavmedianlev, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvWavmedianlev, M("GENERAL_DISABLED"));
        }
    }
}

void Wavelet::linkedgToggled ()
{

    if (multiImage) {
        if (linkedg->get_inconsistent()) {
            linkedg->set_inconsistent (false);
            linkedgConn.block (true);
            linkedg->set_active (false);
            linkedgConn.block (false);
        } else if (lastlinkedg) {
            linkedg->set_inconsistent (true);
        }

        lastlinkedg = linkedg->get_active ();
    }

    if (listener && (multiImage || getEnabled ())) {
        if (linkedg->get_inconsistent()) {
            listener->panelChanged (EvWavlinkedg, M("GENERAL_UNCHANGED"));
        } else if (linkedg->get_active () ) {
            listener->panelChanged (EvWavlinkedg, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvWavlinkedg, M("GENERAL_DISABLED"));
        }
    }
}

void Wavelet::cbenabUpdateUI ()
{
    if (!batchMode) {
        if(cbenab->get_active ()) {
            chanMixerHLFrame->show();
            chanMixerMidFrame->show();
            chanMixerShadowsFrame->show();
            neutrHBox->show();
        } else {
            chanMixerHLFrame->hide();
            chanMixerMidFrame->hide();
            chanMixerShadowsFrame->hide();
            neutrHBox->hide();
        }
    }
}

void Wavelet::cbenabToggled ()
{
    if (multiImage) {
        if (cbenab->get_inconsistent()) {
            cbenab->set_inconsistent (false);
            cbenabConn.block (true);
            cbenab->set_active (false);
            cbenabConn.block (false);
        } else if (lastcbenab) {
            cbenab->set_inconsistent (true);
        }

        lastcbenab = cbenab->get_active ();
    }

    cbenabUpdateUI();

    if (listener && (multiImage || getEnabled ())) {
        if (cbenab->get_inconsistent() ) {
            listener->panelChanged (EvWavcbenab, M("GENERAL_UNCHANGED"));
        } else if (cbenab->get_active () ) {
            listener->panelChanged (EvWavcbenab, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvWavcbenab, M("GENERAL_DISABLED"));
        }
    }
}



void Wavelet::lipstUpdateUI ()
{
    if (!batchMode) {
        if (lipst->get_active ()) {
            NPmethod->show();
            edgesensi->show();
            edgeampli->show();
            labmNP->show();
        } else  {
            NPmethod->hide();
            edgesensi->hide();
            edgeampli->hide();
            labmNP->hide();

        }
    }
}


void Wavelet::lipstToggled ()
{

    if (multiImage) {
        if (lipst->get_inconsistent()) {
            lipst->set_inconsistent (false);
            lipstConn.block (true);
            lipst->set_active (false);
            lipstConn.block (false);
        } else if (lastlipst) {
            lipst->set_inconsistent (true);
        }

        lastlipst = lipst->get_active ();
    }

    lipstUpdateUI();

    if (listener && (multiImage || getEnabled ())) {
        if (lipst->get_inconsistent()) {
            listener->panelChanged (EvWavlipst, M("GENERAL_UNCHANGED"));
        } else if (lipst->get_active ()) {
            listener->panelChanged (EvWavlipst, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvWavlipst, M("GENERAL_DISABLED"));
        }
    }
}

/*
void Wavelet::edgreinfToggled () {

    if (batchMode) {
        if (edgreinf->get_inconsistent()) {
            edgreinf->set_inconsistent (false);
            edgreinfConn.block (true);
            edgreinf->set_active (false);
            edgreinfConn.block (false);
        }
        else if (lastedgreinf)
            edgreinf->set_inconsistent (true);

        lastedgreinf = edgreinf->get_active ();
    }

    if (listener) {
        if (enabled->get_inconsistent ())
            listener->panelChanged (EvWavedgreinf, M("GENERAL_UNCHANGED"));
        else if (enabled->get_active ())
            listener->panelChanged (EvWavedgreinf, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvWavedgreinf, M("GENERAL_DISABLED"));
    }
}
*/
void Wavelet::avoidToggled ()
{

    if (multiImage) {
        if (avoid->get_inconsistent()) {
            avoid->set_inconsistent (false);
            avoidConn.block (true);
            avoid->set_active (false);
            avoidConn.block (false);
        } else if (lastavoid) {
            avoid->set_inconsistent (true);
        }

        lastavoid = avoid->get_active ();
    }

    if (listener && (multiImage || getEnabled ())) {
        if (avoid->get_inconsistent()) {
            listener->panelChanged (EvWavavoid, M("GENERAL_UNCHANGED"));
        } else if (avoid->get_active ()) {
            listener->panelChanged (EvWavavoid, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvWavavoid, M("GENERAL_DISABLED"));
        }
    }
}

void Wavelet::tmrToggled ()
{

    if (multiImage) {
        if (tmr->get_inconsistent()) {
            tmr->set_inconsistent (false);
            tmrConn.block (true);
            tmr->set_active (false);
            tmrConn.block (false);
        } else if (lasttmr) {
            tmr->set_inconsistent (true);
        }

        lasttmr = tmr->get_active ();
    }

    if (listener && (multiImage || getEnabled ())) {
        if (tmr->get_inconsistent()) {
            listener->panelChanged (EvWavtmr, M("GENERAL_UNCHANGED"));
        } else if (tmr->get_active ()) {
            listener->panelChanged (EvWavtmr, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvWavtmr, M("GENERAL_DISABLED"));
        }
    }
}


void Wavelet::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R, G, B;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve

        Color::hsv2rgb01(float(valX), float(valY), 0.5f, R, G, B);
    }
    /*  else if (callerId == 2) {    // cc - bottom bar

        //  float value = (1.f - 0.7f) * float(valX) + 0.7f;
            float value = (1.f - 0.7f) * float(valX) + 0.7f;
            // whole hue range
            // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
            Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
        }
        */
    else if (callerId == 4) {    // LH - bottom bar
        Color::hsv2rgb01(float(valX), 0.5f, float(valY), R, G, B);
    } else if (callerId == 5) {  // HH - bottom bar
        float h = float((valY - 0.5) * 0.3 + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01(h, 0.5f, 0.5f, R, G, B);
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
}
void Wavelet::setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool threshold2add, bool thresadd, bool chroadd, bool chromaadd, bool contrastadd, bool skinadd, bool reschroadd, bool tmrsadd, bool resconadd, bool resconHadd, bool thradd, bool thrHadd, bool skyadd, bool edgradadd, bool edgvaladd, bool strengthadd,  bool gammaadd, bool edgedetectadd, bool edgedetectthradd , bool edgedetectthr2add)
{

    for (int i = 0; i < 9; i++) {
        correction[i]->setAddMode(multiplieradd);
    }

    threshold->setAddMode(thresholdadd);
    skinprotect->setAddMode(skinadd);
    threshold2->setAddMode(threshold2add);
    thres->setAddMode(thresadd);
    chro->setAddMode(chroadd);
    chroma->setAddMode(chromaadd);
    contrast->setAddMode(contrastadd);
    rescon->setAddMode(resconadd);
    resconH->setAddMode(resconHadd);
    reschro->setAddMode(reschroadd);
    tmrs->setAddMode(tmrsadd);
    thr->setAddMode(thradd);
    thrH->setAddMode(thrHadd);
    sky->setAddMode(skyadd);
    edgrad->setAddMode(edgradadd);
    edgval->setAddMode(edgvaladd);
    strength->setAddMode(strengthadd);
    gamma->setAddMode(gammaadd);
    edgedetect->setAddMode(edgedetectadd);
    edgedetectthr->setAddMode(edgedetectthradd);
    edgedetectthr2->setAddMode(edgedetectthr2add);
}


void Wavelet::neutralPressed ()
{
    for (int i = 0; i < 9; i++) {
        correction[i]->setValue(0);
        adjusterChanged(correction[i], 0);
    }
}

void Wavelet::neutralchPressed ()
{

    for (int i = 0; i < 9; i++) {
        correctionch[i]->setValue(0);
        adjusterChanged(correctionch[i], 0);
    }
}


void Wavelet::contrastPlusPressed ()
{

    for (int i = 0; i < 9; i++) {
        int inc = 1 * (9 - i);
        correction[i]->setValue(correction[i]->getValue() + inc);
        adjusterChanged(correction[i], correction[i]->getValue());
    }
}


void Wavelet::contrastMinusPressed ()
{

    for (int i = 0; i < 9; i++) {
        int inc = -1 * (9 - i);
        correction[i]->setValue(correction[i]->getValue() + inc);
        adjusterChanged(correction[i], correction[i]->getValue());
    }
}

void Wavelet::foldAllButMe (GdkEventButton* event, MyExpander *expander)
{
    if (event->button == 3) {
        expsettings->set_expanded(expsettings == expander);
        expcontrast->set_expanded(expcontrast == expander);
        expchroma->set_expanded(expchroma == expander);
        exptoning->set_expanded(exptoning == expander);
        expnoise->set_expanded(expnoise == expander);
        expedge->set_expanded(expedge == expander);
        expgamut->set_expanded(expgamut == expander);
        expresid->set_expanded(expresid == expander);
        expfinal->set_expanded(expfinal == expander);
    }
}

void Wavelet::enableToggled(MyExpander *expander)
{
    if (listener) {
        rtengine::ProcEvent event = NUMOFEVENTS;

        if (expander == expcontrast) {
            event = EvWavenacont;
        } else if (expander == expchroma) {
            event = EvWavenachrom;
        } else if (expander == exptoning) {
            event = EvWavenatoning;
        } else if (expander == expnoise) {
            event = EvWavenanoise;
        } else if (expander == expedge) {
            event = EvWavenaedge;
        } else if (expander == expresid) {
            event = EvWavenares;
        } else if (expander == expfinal) {
            event = EvWavenafin;
        } else
            // unknown expander, returning !
        {
            return;
        }

        if (expander->get_inconsistent()) {
            listener->panelChanged (event, M("GENERAL_UNCHANGED"));
        } else if (expander->getEnabled()) {
            listener->panelChanged (event, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (event, M("GENERAL_DISABLED"));
        }
    }
}

void Wavelet::writeOptions(std::vector<int> &tpOpen)
{
    tpOpen.push_back (expsettings->get_expanded ());
    tpOpen.push_back (expcontrast->get_expanded ());
    tpOpen.push_back (expchroma->get_expanded ());
    tpOpen.push_back (exptoning->get_expanded ());
    tpOpen.push_back (expnoise->get_expanded ());
    tpOpen.push_back (expedge->get_expanded ());
    tpOpen.push_back (expgamut->get_expanded ());
    tpOpen.push_back (expresid->get_expanded ());
    tpOpen.push_back (expfinal->get_expanded ());
}

void Wavelet::updateToolState(std::vector<int> &tpOpen)
{
    if(tpOpen.size() == 9) {
        expsettings->set_expanded(tpOpen.at(0));
        expcontrast->set_expanded(tpOpen.at(1));
        expchroma->set_expanded(tpOpen.at(2));
        exptoning->set_expanded(tpOpen.at(3));
        expnoise->set_expanded(tpOpen.at(4));
        expedge->set_expanded(tpOpen.at(5));
        expgamut->set_expanded(tpOpen.at(6));
        expresid->set_expanded(tpOpen.at(7));
        expfinal->set_expanded(tpOpen.at(8));
    }
}

