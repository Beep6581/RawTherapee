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



Locallab::Locallab ():
    FoldableToolPanel (this, "locallab", M ("TP_LOCALLAB_LABEL"), false, true),
    EditSubscriber (ET_OBJECTS), lastObject (-1), draggedPointOldAngle (-1000.),
    expcolor (new MyExpander (true, M ("TP_LOCALLAB_COFR"))),
    expblur (new MyExpander (true, M ("TP_LOCALLAB_BLUFR"))),
    exptonemap (new MyExpander (true, M ("TP_LOCALLAB_TM"))),
    expreti (new MyExpander (true, M ("TP_LOCALLAB_RETI"))),
    expsharp (new MyExpander (true, M ("TP_LOCALLAB_SHARP"))),
    expcbdl (new MyExpander (true, M ("TP_LOCALLAB_CBDL"))),
    expdenoi (new MyExpander (true, M ("TP_LOCALLAB_DENOIS"))),
    expsettings (new MyExpander (false, M ("TP_LOCALLAB_SETTINGS"))),

    llCurveEditorG (new CurveEditorGroup (options.lastlocalCurvesDir, M ("TP_LOCALLAB_LUM"))),
    LocalcurveEditorgainT (new CurveEditorGroup (options.lastlocalCurvesDir, M ("TP_LOCALLAB_TRANSMISSIONGAIN"))),
    LocalcurveEditorgainTrab (new CurveEditorGroup (options.lastlocalCurvesDir, M ("TP_LOCALLAB_TRANSMISSIONGAINRAB"))),


    anbspot (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_ANBSPOT"), 0, 1, 1, 0))),
    locX (Gtk::manage (new Adjuster (M ("TP_LOCAL_WIDTH"), 0, 1500, 1, 250))),
    locXL (Gtk::manage (new Adjuster (M ("TP_LOCAL_WIDTH_L"), 0, 1500, 1, 250))),
    degree (Gtk::manage (new Adjuster (M ("TP_LOCAL_DEGREE"), -180, 180, 1, 0))),
    locY (Gtk::manage (new Adjuster (M ("TP_LOCAL_HEIGHT"), 0, 1500, 1, 250))),
    locYT (Gtk::manage (new Adjuster (M ("TP_LOCAL_HEIGHT_T"), 0, 1500, 1, 250))),
    centerX (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CENTER_X"), -1000, 1000, 1, 0))),
    centerY (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CENTER_Y"), -1000, 1000, 1, 0))),
    circrad (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CIRCRADIUS"), 4, 150, 1, 18))),
    thres (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_THRES"), 1, 35, 1, 18))),
    proxi (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_PROXI"), 0, 60, 1, 20))),
    lightness (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_LIGHTNESS"), -100, 100, 1, 0))),
    contrast (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0))),
    chroma (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CHROMA"), -100, 150, 1, 0))),
    sensi (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    radius (Gtk::manage ( new Adjuster (M ("TP_LOCALLAB_RADIUS"), 1, 100, 1, 1) )),
    strength (Gtk::manage ( new Adjuster (M ("TP_LOCALLAB_STRENGTH"), 0, 100, 1, 0) )),
    sensibn (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SENSIBN"), 0, 100, 1, 60))),
    transit (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_TRANSIT"), 5, 95, 1, 60))),
    stren (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_STREN"), -50, 100, 1, 0))),
    gamma (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_GAM"), 80, 150, 1, 100))),
    estop (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_ESTOP"), 10, 400, 1, 140))),
    scaltm (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SCALTM"), 1, 100, 1, 3))),
    rewei (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_REWEI"), 0, 9, 1, 0))),
    sensitm (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    str (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_STR"), 0, 100, 1, 0))),
    neigh (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_NEIGH"), 14, 150, 1, 50))),
    vart (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_VART"), 50, 500, 1, 200))),
    chrrt (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CHRRT"), 0, 100, 1, 0))),
    sensih (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SENSIH"), 0, 100, 1, 19))),
    retrab (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_RETRAB"), 0, 10000, 1, 500))),
    threshold (Gtk::manage ( new Adjuster (M ("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 100, 1, 20) )),
    sensicb (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SENSICB"), 0, 100, 1, 19))),
    sharradius (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SHARRADIUS"), 42, 250, 1, 4))),
    sharamount (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 75))),
    shardamping (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 75))),
    shariter (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sensisha (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    noiselumf (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_NOISELUMFINE"), 0, 100, 1, 0))),
    noiselumc (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_NOISELUMCOARSE"), 0, 100, 1, 0))),
    noisechrof (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_NOISECHROFINE"), 0, 100, 1, 0))),
    noisechroc (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_NOISECHROCOARSE"), 0, 100, 1, 0))),
    hueref (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_HUEREF"), -3.15, 3.15, 0.01, 0))),
    chromaref (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CHROMAREF"), 0, 200, 0.01, 0))),
    lumaref (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_LUMAMAREF"), 0, 100, 0.01, 0))),

    Smethod (Gtk::manage (new MyComboBoxText ())),
    qualityMethod (Gtk::manage (new MyComboBoxText ())),
    retinexMethod (Gtk::manage (new MyComboBoxText ())),
    qualitycurveMethod (Gtk::manage (new MyComboBoxText ())),

    shapeFrame (Gtk::manage (new Gtk::Frame (M ("TP_LOCALLAB_SHFR")))),
    artifFrame (Gtk::manage (new Gtk::Frame (M ("TP_LOCALLAB_ARTIF")))),
    superFrame (Gtk::manage (new Gtk::Frame ())),

    artifVBox (Gtk::manage (new Gtk::VBox ())),
    shapeVBox (Gtk::manage (new Gtk::VBox ())),
    tmBox (Gtk::manage (new Gtk::VBox())),
    retiBox (Gtk::manage (new Gtk::VBox())),
    colorVBox (Gtk::manage ( new Gtk::VBox())),
    blurrVBox (Gtk::manage ( new Gtk::VBox())),
    sharpVBox (Gtk::manage ( new Gtk::VBox())),
    cbdlVBox (Gtk::manage ( new Gtk::VBox())),
    denoisVBox (Gtk::manage ( new Gtk::VBox())),
    superVBox (Gtk::manage (new Gtk::VBox ())),


    labmdh (Gtk::manage (new Gtk::Label (M ("TP_LOCRETI_METHOD") + ":"))),
    labqual (Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_QUAL_METHOD") + ":"))),
    labqualcurv (Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_QUALCURV_METHOD") + ":"))),
    labmS (Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_STYPE") + ":"))),

    ctboxS (Gtk::manage (new Gtk::HBox ())),
    dhbox (Gtk::manage (new Gtk::HBox ())),
    qualbox (Gtk::manage (new Gtk::HBox ())),
    qualcurvbox (Gtk::manage (new Gtk::HBox ())),

    avoid (Gtk::manage (new Gtk::CheckButton (M ("TP_LOCALLAB_AVOID")))),
    activlum (Gtk::manage (new Gtk::CheckButton (M ("TP_LOCALLAB_ACTIV")))),
    invers (Gtk::manage (new Gtk::CheckButton (M ("TP_LOCALLAB_INVERS")))),
    curvactiv (Gtk::manage (new Gtk::CheckButton (M ("TP_LOCALLAB_CURV")))),
    inversrad (Gtk::manage (new Gtk::CheckButton (M ("TP_LOCALLAB_INVERS")))),
    inversret (Gtk::manage (new Gtk::CheckButton (M ("TP_LOCALLAB_INVERS")))),
    inverssha (Gtk::manage (new Gtk::CheckButton (M ("TP_LOCALLAB_INVERS"))))


{
    CurveListener::setMulti (true);
    ProcParams params;
    editHBox = Gtk::manage (new Gtk::HBox());
    edit = Gtk::manage (new Gtk::ToggleButton());
    edit->add (*Gtk::manage (new RTImage ("editmodehand.png")));
    edit->set_tooltip_text (M ("EDIT_OBJECT_TOOLTIP"));
    editConn = edit->signal_toggled().connect ( sigc::mem_fun (*this, &Locallab::editToggled) );
    editHBox->pack_start (*edit, Gtk::PACK_SHRINK, 0);
    pack_start (*editHBox, Gtk::PACK_SHRINK, 0);
    int realnbspot;


    realnbspot = options.rtSettings.nspot;

    nbspot = Gtk::manage (new Adjuster (M ("TP_LOCALLAB_NBSPOT"), 1, realnbspot, 1, 1));

    if (options.rtSettings.locdelay) {

        if (nbspot->delay < 200) {
            nbspot->delay = 200;
        }
    }


    nbspot->setAdjusterListener (this);
    nbspot->set_tooltip_text (M ("TP_LOCALLAB_NBSPOT_TOOLTIP"));


    anbspot->setAdjusterListener (this);
    anbspot->set_tooltip_text (M ("TP_LOCALLAB_ANBSPOT_TOOLTIP"));

    shapeFrame->set_label_align (0.025, 0.5);

    expsettings->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Locallab::foldAllButMe), expsettings) );


    expcolor->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Locallab::foldAllButMe), expcolor) );
    enablecolorConn = expcolor->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Locallab::enableToggled), expcolor) );

    expblur->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Locallab::foldAllButMe), expblur) );
    enableblurConn = expblur->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Locallab::enableToggled), expblur) );

    exptonemap->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Locallab::foldAllButMe), exptonemap) );
    enabletonemapConn = exptonemap->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Locallab::enableToggled), exptonemap) );

    expreti->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Locallab::foldAllButMe), expreti) );
    enableretiConn = expreti->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Locallab::enableToggled), expreti) );

    expsharp->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Locallab::foldAllButMe), expsharp) );
    enablesharpConn = expsharp->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Locallab::enableToggled), expsharp) );

    expcbdl->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Locallab::foldAllButMe), expcbdl) );
    enablecbdlConn = expcbdl->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Locallab::enableToggled), expcbdl) );

    expdenoi->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Locallab::foldAllButMe), expdenoi) );
    enabledenoiConn = expdenoi->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Locallab::enableToggled), expdenoi) );


    ctboxS->pack_start (*labmS, Gtk::PACK_SHRINK, 4);
    ctboxS->set_tooltip_markup (M ("TP_LOCALLAB_STYPE_TOOLTIP"));

    Smethod->append (M ("TP_LOCALLAB_IND"));
    Smethod->append (M ("TP_LOCALLAB_SYM"));
    Smethod->append (M ("TP_LOCALLAB_INDSL"));
    Smethod->append (M ("TP_LOCALLAB_SYMSL"));
    Smethod->set_active (0);
    Smethodconn = Smethod->signal_changed().connect ( sigc::mem_fun (*this, &Locallab::SmethodChanged) );

    //locX->set_tooltip_text (M("TP_LOCAL_WIDTH_TOOLTIP"));
    locX->setAdjusterListener (this);

    //locX->set_tooltip_text (M("TP_LOCAL_WIDTH_TOOLTIP"));
    locXL->setAdjusterListener (this);

    //degree->set_tooltip_text (M("TP_LOCAL_DEGREE_TOOLTIP"));
    degree->setAdjusterListener (this);

    //locY->set_tooltip_text (M("TP_LOCAL_HEIGHT_TOOLTIP"));
    locY->setAdjusterListener (this);

    //locY->set_tooltip_text (M("TP_LOCAL_HEIGHT_TOOLTIP"));
    locYT->setAdjusterListener (this);

    //centerX->set_tooltip_text (M("TP_LOCALLAB_CENTER_X_TOOLTIP"));
    centerX->setAdjusterListener (this);

    //centerY->set_tooltip_text (M("TP_LOCALLAB_CENTER_Y_TOOLTIP"));
    centerY->setAdjusterListener (this);

    circrad->setAdjusterListener (this);


    qualityMethod->append (M ("TP_LOCALLAB_STD"));
    qualityMethod->append (M ("TP_LOCALLAB_ENH"));
    qualityMethod->append (M ("TP_LOCALLAB_ENHDEN"));
    qualityMethod->set_active (0);
    qualityMethodConn = qualityMethod->signal_changed().connect ( sigc::mem_fun (*this, &Locallab::qualityMethodChanged) );
    qualityMethod->set_tooltip_markup (M ("TP_LOCALLAB_METHOD_TOOLTIP"));

    thres->setAdjusterListener (this);

    proxi->setAdjusterListener (this);
    std::vector<GradientMilestone> milestones;
    std::vector<double> defaultCurve;
    std::vector<double> defaultCurve2;
    std::vector<double> defaultCurve2rab;
    std::vector<double> defaultCurve3;
    std::vector<double> defaultCurve4;

    irg   = Gtk::manage (new RTImage ("Chanmixer-RG.png"));

    qualitycurveMethod->append (M ("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append (M ("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->append (M ("TP_LOCALLAB_CURVENH"));
    qualitycurveMethod->set_active (0);
    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect ( sigc::mem_fun (*this, &Locallab::qualitycurveMethodChanged) );
    qualitycurveMethod->set_tooltip_markup (M ("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));



    llCurveEditorG->setCurveListener (this);

    rtengine::LocallabParams::getDefaultLLCurve (defaultCurve);
    llshape = static_cast<DiagonalCurveEditor*> (llCurveEditorG->addCurve (CT_Diagonal, "L(L)"));
    llshape->setResetCurve (DiagonalCurveType (defaultCurve.at (0)), defaultCurve);
    llshape->setTooltip (M ("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    milestones.push_back ( GradientMilestone (0., 0., 0., 0.) );
    milestones.push_back ( GradientMilestone (1., 1., 1., 1.) );
    llshape->setBottomBarBgGradient (milestones);
    llshape->setLeftBarBgGradient (milestones);

    rtengine::LocallabParams::getDefaultCCCurve (defaultCurve4);
    ccshape = static_cast<DiagonalCurveEditor*> (llCurveEditorG->addCurve (CT_Diagonal, "C(C)"));
    ccshape->setResetCurve (DiagonalCurveType (defaultCurve4.at (0)), defaultCurve4);
    ccshape->setTooltip (M ("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    milestones.push_back ( GradientMilestone (0., 0., 0., 0.) );
    milestones.push_back ( GradientMilestone (1., 1., 1., 1.) );
    ccshape->setBottomBarBgGradient (milestones);
    ccshape->setLeftBarBgGradient (milestones);

    rtengine::LocallabParams::getDefaultLHCurve (defaultCurve3);

    LHshape = static_cast<FlatCurveEditor*> (llCurveEditorG->addCurve (CT_Flat, "L(H)", nullptr, false, true));

    LHshape->setIdentityValue (0.);
    LHshape->setResetCurve (FlatCurveType (defaultCurve3.at (0)), defaultCurve3);
    LHshape->setTooltip (M ("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    LHshape->setCurveColorProvider (this, 1);
    milestones.clear();

    for (int i = 0; i < 7; i++) {
        float R, G, B;
        float x = float (i) * (1.0f / 6.0);
        Color::hsv2rgb01 (x, 0.5f, 0.5f, R, G, B);
        milestones.push_back ( GradientMilestone (double (x), double (R), double (G), double (B)) );
    }

    LHshape->setBottomBarBgGradient (milestones);


    llCurveEditorG->curveListComplete();



    //lightness->set_tooltip_text (M("TP_LOCALLAB_LIGHTNESS_TOOLTIP"));
    lightness->setAdjusterListener (this);

    //contrast->set_tooltip_text (M("TP_LOCALLAB_CONTRAST_TOOLTIP"));
    contrast->setAdjusterListener (this);

    //chroma->set_tooltip_text (M("TP_LOCALLAB_CHROMA_TOOLTIP"));
    chroma->setAdjusterListener (this);

    sensi->set_tooltip_text (M ("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensi->setAdjusterListener (this);

    //radius->set_tooltip_text (M("TP_LOCALLAB_RADIUS_TOOLTIP"));
    radius->setAdjusterListener (this);
    //radius->set_tooltip_text (M("TP_LOCALLAB_RADIUS_TOOLTIP"));
    strength->setAdjusterListener (this);


    sensibn->set_tooltip_text (M ("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensibn->setAdjusterListener (this);

    activlum->set_active (false);
    activlumConn  = activlum->signal_toggled().connect ( sigc::mem_fun (*this, &Locallab::activlumChanged) );

    transit->set_tooltip_text (M ("TP_LOCALLAB_TRANSIT_TOOLTIP"));
    transit->setAdjusterListener (this);

    invers->set_active (false);
    inversConn  = invers->signal_toggled().connect ( sigc::mem_fun (*this, &Locallab::inversChanged) );

    curvactiv->set_active (false);
    curvactivConn  = curvactiv->signal_toggled().connect ( sigc::mem_fun (*this, &Locallab::curvactivChanged) );

    inversrad->set_active (false);
    inversradConn  = inversrad->signal_toggled().connect ( sigc::mem_fun (*this, &Locallab::inversradChanged) );

    inversret->set_active (false);
    inversretConn  = inversret->signal_toggled().connect ( sigc::mem_fun (*this, &Locallab::inversretChanged) );

//tone mapping local

    stren->setAdjusterListener (this);

    gamma->setAdjusterListener (this);

    estop->setAdjusterListener (this);

    scaltm->setAdjusterListener (this);

    rewei->setAdjusterListener (this);

    sensitm->set_tooltip_text (M ("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensitm->setAdjusterListener (this);

//end TM


//retinex local

    dhbox->pack_start (*labmdh, Gtk::PACK_SHRINK, 1);

    retinexMethod->append (M ("TP_RETINEX_LOW"));
    retinexMethod->append (M ("TP_RETINEX_UNIFORM"));
    retinexMethod->append (M ("TP_RETINEX_HIGH"));
    retinexMethod->set_active (0);
    retinexMethodConn = retinexMethod->signal_changed().connect ( sigc::mem_fun (*this, &Locallab::retinexMethodChanged) );
    retinexMethod->set_tooltip_markup (M ("TP_LOCRETI_METHOD_TOOLTIP"));

    str->setAdjusterListener (this);
    neigh->setAdjusterListener (this);
    vart->setAdjusterListener (this);
    chrrt->setAdjusterListener (this);
    sensih->set_tooltip_text (M ("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensih->setAdjusterListener (this);
    retrab->setAdjusterListener (this);


    LocalcurveEditorgainT->setCurveListener (this);
    rtengine::LocallabParams::getDefaultLocalgainCurveT (defaultCurve2);


    cTgainshape = static_cast<FlatCurveEditor*> (LocalcurveEditorgainT->addCurve (CT_Flat, "", nullptr, false, false));

    cTgainshape->setIdentityValue (0.);
    cTgainshape->setResetCurve (FlatCurveType (defaultCurve2.at (0)), defaultCurve2);
    cTgainshape->setTooltip (M ("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));

    LocalcurveEditorgainTrab->setCurveListener (this);

    rtengine::LocallabParams::getDefaultLocalgainCurveTrab (defaultCurve2rab);


    cTgainshaperab = static_cast<FlatCurveEditor*> (LocalcurveEditorgainTrab->addCurve (CT_Flat, "", nullptr, false, false));


    cTgainshaperab->setIdentityValue (0.);
    cTgainshaperab->setResetCurve (FlatCurveType (defaultCurve2rab.at (0)), defaultCurve2rab);
    cTgainshaperab->setTooltip (M ("TP_RETINEX_GAINTRANSMISSIONRAB_TOOLTIP"));

    LocalcurveEditorgainT->curveListComplete();
    LocalcurveEditorgainT->show();
    LocalcurveEditorgainTrab->curveListComplete();
    LocalcurveEditorgainTrab->show();


// end reti

    avoid->set_active (false);
    avoidConn  = avoid->signal_toggled().connect ( sigc::mem_fun (*this, &Locallab::avoidChanged) );
    shapeVBox->pack_start (*nbspot);
    pack_start (*anbspot);

    hueref->setAdjusterListener (this);
    chromaref->setAdjusterListener (this);
    lumaref->setAdjusterListener (this);

    pack_start (*hueref);
    pack_start (*chromaref);
    pack_start (*lumaref);

    anbspot->hide();//keep anbspot  - i used it to test diffrent algo...
    hueref->hide();
    chromaref->hide();
    lumaref->hide();
    ctboxS->pack_start (*Smethod);
    shapeVBox->pack_start (*ctboxS);

    shapeVBox->pack_start (*locX);
    shapeVBox->pack_start (*locXL);
    //pack_start (*degree);
    shapeVBox->pack_start (*locY);
    shapeVBox->pack_start (*locYT);
    shapeVBox->pack_start (*centerX);
    shapeVBox->pack_start (*centerY);
    shapeVBox->pack_start (*circrad);
    qualbox->pack_start (*labqual, Gtk::PACK_SHRINK, 4);
    qualbox->pack_start (*qualityMethod);
    shapeVBox->pack_start (*qualbox);
    shapeVBox->pack_start (*transit);

    artifFrame->set_label_align (0.025, 0.5);
    artifFrame->set_tooltip_text (M ("TP_LOCALLAB_ARTIF_TOOLTIP"));


    artifVBox->pack_start (*thres);
    artifVBox->pack_start (*proxi);
    artifFrame->add (*artifVBox);
    shapeVBox->pack_start (*artifFrame);

    expsettings->add (*shapeVBox);
    expsettings->setLevel (2);
    pack_start (*expsettings);



    Gtk::HBox * buttonBox1 = Gtk::manage (new Gtk::HBox (true, 10));

    Gtk::Button * lumacontrastMinusButton = Gtk::manage (new Gtk::Button (M ("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")));
    buttonBox1->pack_start (*lumacontrastMinusButton);
    lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect ( sigc::mem_fun (*this, &Locallab::lumacontrastMinusPressed));

    Gtk::Button * lumaneutralButton = Gtk::manage (new Gtk::Button (M ("TP_DIRPYREQUALIZER_LUMANEUTRAL")));
    buttonBox1->pack_start (*lumaneutralButton);
    lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect ( sigc::mem_fun (*this, &Locallab::lumaneutralPressed));

    Gtk::Button * lumacontrastPlusButton = Gtk::manage (new Gtk::Button (M ("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS")));
    buttonBox1->pack_start (*lumacontrastPlusButton);
    lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect ( sigc::mem_fun (*this, &Locallab::lumacontrastPlusPressed));

    cbdlVBox->pack_start (*buttonBox1);

    for (int i = 0; i < 5; i++) {
        Glib::ustring ss;
        ss = Glib::ustring::format (i);

        if     (i == 0) {
            ss += Glib::ustring::compose (" (%1)", M ("TP_DIRPYREQUALIZER_LUMAFINEST"));
        } else if (i == 4) {
            ss += Glib::ustring::compose (" (%1)", M ("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        }

        multiplier[i] = Gtk::manage ( new Adjuster (ss, 0, 400, 1, 100) );
        multiplier[i]->setAdjusterListener (this);
        cbdlVBox->pack_start (*multiplier[i]);
    }

    Gtk::HSeparator *separator3 = Gtk::manage (new  Gtk::HSeparator());
    cbdlVBox->pack_start (*separator3, Gtk::PACK_SHRINK, 2);

    threshold->setAdjusterListener (this);
    cbdlVBox->pack_start (*threshold);

    sensicb->set_tooltip_text (M ("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensicb->setAdjusterListener (this);
    cbdlVBox->pack_start (*sensicb);

    sharradius->setAdjusterListener (this);

    sharamount->setAdjusterListener (this);

    shardamping->setAdjusterListener (this);

    shariter->setAdjusterListener (this);


    sensisha->set_tooltip_text (M ("TP_LOCALLAB_SENSIS_TOOLTIP"));
    sensisha->setAdjusterListener (this);

    inverssha->set_active (false);
    inversshaConn  = inverssha->signal_toggled().connect ( sigc::mem_fun (*this, &Locallab::inversshaChanged) );

    sharpVBox->pack_start (*sharradius);
    sharpVBox->pack_start (*sharamount);
    sharpVBox->pack_start (*shardamping);
    sharpVBox->pack_start (*shariter);
    sharpVBox->pack_start (*sensisha);
    sharpVBox->pack_start (*inverssha);


    noiselumf->setAdjusterListener (this);

    noiselumc->setAdjusterListener (this);

    noisechrof->setAdjusterListener (this);

    noisechroc->setAdjusterListener (this);

    denoisVBox->pack_start (*noiselumf);
    denoisVBox->pack_start (*noiselumc);
    denoisVBox->pack_start (*noisechrof);
    denoisVBox->pack_start (*noisechroc);

    neutrHBox1 = Gtk::manage (new Gtk::HBox ());

    neutral1 = Gtk::manage (new Gtk::Button (M ("TP_LOCALLAB_NEUTRAL")));
    RTImage *resetImg1 = Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
    neutral1->set_image (*resetImg1);
    neutral1->set_tooltip_text (M ("TP_LOCALLAB_NEUTRAL_TIP"));
    neutralconn1 = neutral1->signal_pressed().connect ( sigc::mem_fun (*this, &Locallab::neutral_pressed) );
    neutral1->show();
    neutrHBox1->pack_start (*neutral1);
    pack_start (*neutrHBox1);

    superFrame->set_label_align (0.025, 0.5);
    Gtk::VBox *superVBox = Gtk::manage ( new Gtk::VBox());
    superFrame->set_label_widget (*curvactiv);


    superVBox->pack_start (*lightness);
    superVBox->pack_start (*contrast);
    superFrame->add (*superVBox);
    colorVBox->pack_start (*superFrame);

    colorVBox->pack_start (*chroma);
    colorVBox->pack_start (*sensi);

    qualcurvbox->pack_start (*labqualcurv, Gtk::PACK_SHRINK, 4);
    qualcurvbox->pack_start (*qualitycurveMethod);

    colorVBox->pack_start (*qualcurvbox);


    colorVBox->pack_start (*llCurveEditorG, Gtk::PACK_SHRINK, 2);
    colorVBox->pack_start (*invers);

    expcolor->add (*colorVBox);
    expcolor->setLevel (2);
    pack_start (*expcolor);

    blurrVBox->pack_start (*radius);
    blurrVBox->pack_start (*strength);
    blurrVBox->pack_start (*sensibn);
    blurrVBox->pack_start (*activlum);

    blurrVBox->pack_start (*inversrad);
    expblur->add (*blurrVBox);
    expblur->setLevel (2);
    pack_start (*expblur);

    tmBox->pack_start (*stren);
    tmBox->pack_start (*gamma);
    tmBox->pack_start (*estop);
    tmBox->pack_start (*scaltm);
    tmBox->pack_start (*rewei);
    tmBox->pack_start (*sensitm);

    exptonemap->add (*tmBox);
    exptonemap->setLevel (2);
    pack_start (*exptonemap);


    retiBox->pack_start (*retinexMethod);
    retiBox->pack_start (*str);
    retiBox->pack_start (*chrrt);
    retiBox->pack_start (*neigh);
    retiBox->pack_start (*vart);
    retiBox->pack_start (*sensih);
    retiBox->pack_start (*retrab);

    retiBox->pack_start (*LocalcurveEditorgainTrab, Gtk::PACK_SHRINK, 4);

    retiBox->pack_start (*LocalcurveEditorgainT, Gtk::PACK_SHRINK, 4);
    retiBox->pack_start (*inversret);

    expreti->add (*retiBox);
    expreti->setLevel (2);
    pack_start (*expreti);


    expsharp->add (*sharpVBox);
    expsharp->setLevel (2);
    pack_start (*expsharp);

    expcbdl->add (*cbdlVBox);
    expcbdl->setLevel (2);
    pack_start (*expcbdl);

    expdenoi->add (*denoisVBox);
    expdenoi->setLevel (2);
    pack_start (*expdenoi);


//    pack_start (*transit);
    pack_start (*avoid);//keep avoid clor shift in case of

    neutrHBox = Gtk::manage (new Gtk::HBox ());

    neutral = Gtk::manage (new Gtk::Button (M ("TP_LOCALLAB_NEUTRAL")));
    RTImage *resetImg = Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
    neutral->set_image (*resetImg);
    neutral->set_tooltip_text (M ("TP_LOCALLAB_NEUTRAL_TIP"));
    neutralconn = neutral->signal_pressed().connect ( sigc::mem_fun (*this, &Locallab::neutral_pressed) );
    neutral->show();
    neutrHBox->pack_start (*neutral);
    pack_start (*neutrHBox);


    // Instantiating the Editing geometry; positions will be initialized later
    Line  *hLine, *vLine, *locYLine[2], *locXLine[2];
    Circle *centerCircle;
//   Arcellipse *oneellipse;

    Beziers *onebeziers[3];
    Beziers *twobeziers[3];
    Beziers *thrbeziers[3];
    Beziers *foubeziers[3];
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

        twobeziers[0] = new Beziers();
        twobeziers[0]->datum = Geometry::IMAGE;
        twobeziers[0]->innerLineWidth = innw;

        twobeziers[1] = new Beziers();
        twobeziers[1]->datum = Geometry::IMAGE;
        twobeziers[1]->innerLineWidth = innw;

        twobeziers[2] = new Beziers();
        twobeziers[2]->datum = Geometry::IMAGE;
        twobeziers[2]->innerLineWidth = innw;

        thrbeziers[0] = new Beziers();
        thrbeziers[0]->datum = Geometry::IMAGE;
        thrbeziers[0]->innerLineWidth = innw;

        thrbeziers[1] = new Beziers();
        thrbeziers[1]->datum = Geometry::IMAGE;
        thrbeziers[1]->innerLineWidth = innw;

        thrbeziers[2] = new Beziers();
        thrbeziers[2]->datum = Geometry::IMAGE;
        thrbeziers[2]->innerLineWidth = innw;

        foubeziers[0] = new Beziers();
        foubeziers[0]->datum = Geometry::IMAGE;
        foubeziers[0]->innerLineWidth = innw;

        foubeziers[1] = new Beziers();
        foubeziers[1]->datum = Geometry::IMAGE;
        foubeziers[1]->innerLineWidth = innw;

        foubeziers[2] = new Beziers();
        foubeziers[2]->datum = Geometry::IMAGE;
        foubeziers[2]->innerLineWidth = innw;
    }

    // oneellipse->radiusInImageSpace = true;
    // oneellipse->radius = locX->getValue();
    // oneellipse->filled = false;

    EditSubscriber::visibleGeometry.push_back ( locXLine[0] );
    EditSubscriber::visibleGeometry.push_back ( locXLine[1] );
    EditSubscriber::visibleGeometry.push_back ( locYLine[0] );
    EditSubscriber::visibleGeometry.push_back ( locYLine[1] );
    EditSubscriber::visibleGeometry.push_back ( centerCircle );

    if (options.showdelimspot) {
        EditSubscriber::visibleGeometry.push_back ( onebeziers[0] );
        EditSubscriber::visibleGeometry.push_back ( onebeziers[1] );
        EditSubscriber::visibleGeometry.push_back ( onebeziers[2] );
        EditSubscriber::visibleGeometry.push_back ( twobeziers[0] );
        EditSubscriber::visibleGeometry.push_back ( twobeziers[1] );
        EditSubscriber::visibleGeometry.push_back ( twobeziers[2] );
        EditSubscriber::visibleGeometry.push_back ( thrbeziers[0] );
        EditSubscriber::visibleGeometry.push_back ( thrbeziers[1] );
        EditSubscriber::visibleGeometry.push_back ( thrbeziers[2] );
        EditSubscriber::visibleGeometry.push_back ( foubeziers[0] );
        EditSubscriber::visibleGeometry.push_back ( foubeziers[1] );
        EditSubscriber::visibleGeometry.push_back ( foubeziers[2] );
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

        twobeziers[0] = new Beziers();
        twobeziers[0]->datum = Geometry::IMAGE;
        twobeziers[0]->innerLineWidth = innw;

        twobeziers[1] = new Beziers();
        twobeziers[1]->datum = Geometry::IMAGE;
        twobeziers[1]->innerLineWidth = innw;

        twobeziers[2] = new Beziers();
        twobeziers[2]->datum = Geometry::IMAGE;
        twobeziers[2]->innerLineWidth = innw;

        thrbeziers[0] = new Beziers();
        thrbeziers[0]->datum = Geometry::IMAGE;
        thrbeziers[0]->innerLineWidth = innw;

        thrbeziers[1] = new Beziers();
        thrbeziers[1]->datum = Geometry::IMAGE;
        thrbeziers[1]->innerLineWidth = innw;

        thrbeziers[2] = new Beziers();
        thrbeziers[2]->datum = Geometry::IMAGE;
        thrbeziers[2]->innerLineWidth = innw;

        foubeziers[0] = new Beziers();
        foubeziers[0]->datum = Geometry::IMAGE;
        foubeziers[0]->innerLineWidth = innw;

        foubeziers[1] = new Beziers();
        foubeziers[1]->datum = Geometry::IMAGE;
        foubeziers[1]->innerLineWidth = innw;

        foubeziers[2] = new Beziers();
        foubeziers[2]->datum = Geometry::IMAGE;
        foubeziers[2]->innerLineWidth = innw;
    }

//   oneellipse->radiusInImageSpace = true;
//   oneellipse->radius = 10;//locX->getValue();
//    oneellipse->filled = false;

    EditSubscriber::mouseOverGeometry.push_back ( locXLine[0] );
    EditSubscriber::mouseOverGeometry.push_back ( locXLine[1] );

    EditSubscriber::mouseOverGeometry.push_back ( locYLine[0] );
    EditSubscriber::mouseOverGeometry.push_back ( locYLine[1] );

    EditSubscriber::mouseOverGeometry.push_back ( centerCircle );

    if (options.showdelimspot) {
        EditSubscriber::mouseOverGeometry.push_back ( onebeziers[0] );
        EditSubscriber::mouseOverGeometry.push_back ( onebeziers[1] );
        EditSubscriber::mouseOverGeometry.push_back ( onebeziers[2] );
        EditSubscriber::mouseOverGeometry.push_back ( twobeziers[0] );
        EditSubscriber::mouseOverGeometry.push_back ( twobeziers[1] );
        EditSubscriber::mouseOverGeometry.push_back ( twobeziers[2] );
        EditSubscriber::mouseOverGeometry.push_back ( thrbeziers[0] );
        EditSubscriber::mouseOverGeometry.push_back ( thrbeziers[1] );
        EditSubscriber::mouseOverGeometry.push_back ( thrbeziers[2] );
        EditSubscriber::mouseOverGeometry.push_back ( foubeziers[0] );
        EditSubscriber::mouseOverGeometry.push_back ( foubeziers[1] );
        EditSubscriber::mouseOverGeometry.push_back ( foubeziers[2] );
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
void Locallab::foldAllButMe (GdkEventButton* event, MyExpander *expander)
{
    if (event->button == 3) {
        expsettings->set_expanded (expsettings == expander);
        expcolor->set_expanded (expcolor == expander);
        expblur->set_expanded (expblur == expander);
        exptonemap->set_expanded (exptonemap == expander);
        expreti->set_expanded (expreti == expander);
        expsharp->set_expanded (expsharp == expander);
        expcbdl->set_expanded (expcbdl == expander);
        expdenoi->set_expanded (expdenoi == expander);

    }
}

void Locallab::enableToggled (MyExpander *expander)
{
    if (listener) {
        rtengine::ProcEvent event = NUMOFEVENTS;

        if (expander == expcolor) {
            event = EvLocenacolor;
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
            listener->panelChanged (event, M ("GENERAL_UNCHANGED"));
        } else if (expander->getEnabled()) {
            listener->panelChanged (event, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (event, M ("GENERAL_DISABLED"));
        }

    }
}

void Locallab::writeOptions (std::vector<int> &tpOpen)
{
    tpOpen.push_back (expsettings->get_expanded ());
    tpOpen.push_back (expcolor->get_expanded ());
    tpOpen.push_back (expblur->get_expanded ());
    tpOpen.push_back (exptonemap->get_expanded ());
    tpOpen.push_back (expreti->get_expanded ());
    tpOpen.push_back (expsharp->get_expanded ());
    tpOpen.push_back (expcbdl->get_expanded ());
    tpOpen.push_back (expdenoi->get_expanded ());

}

void Locallab::updateToolState (std::vector<int> &tpOpen)
{
    if (tpOpen.size() == 8) {
        expsettings->set_expanded (tpOpen.at (0));
        expcolor->set_expanded (tpOpen.at (1));
        expblur->set_expanded (tpOpen.at (2));
        exptonemap->set_expanded (tpOpen.at (3));
        expreti->set_expanded (tpOpen.at (4));
        expsharp->set_expanded (tpOpen.at (5));
        expcbdl->set_expanded (tpOpen.at (6));
        expdenoi->set_expanded (tpOpen.at (7));
    }
}



void Locallab::neutral_pressed ()
{
    Smethod->set_active (0);
    locX->resetValue (false);
    locXL->resetValue (false);
    locY->resetValue (false);
    locYT->resetValue (false);
    centerX->resetValue (false);
    centerY->resetValue (false);
    circrad->resetValue (false);
    qualityMethod->set_active (0);
    qualitycurveMethod->set_active (0);
    thres->resetValue (false);
    proxi->resetValue (false);
    lightness->resetValue (false);
    chroma->resetValue (false);
    contrast->resetValue (false);
    sensi->resetValue (false);
    radius->resetValue (false);
    strength->resetValue (false);
    transit->resetValue (false);
    sensibn->resetValue (false);
    invers->set_active (false);
    curvactiv->set_active (false);
    inversrad->set_active (false);
    inversret->set_active (false);
    stren->resetValue (false);
    gamma->resetValue (false);
    estop->resetValue (false);
    scaltm->resetValue (false);
    rewei->resetValue (false);
    sensitm->resetValue (false);
    retinexMethod->set_active (2);
    str->resetValue (false);
    neigh->resetValue (false);
    vart->resetValue (false);
    chrrt->resetValue (false);
    sensih->resetValue (false);
    retrab->resetValue (false);
//    cTgainshape->reset();
//  cTgainshape->setCurve (creti);
    avoid->set_active (false);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->resetValue (false);
    }

    threshold->resetValue (false);
    sensicb->resetValue (false);
    sharradius->resetValue (false);
    sharamount->resetValue (false);
    shardamping->resetValue (false);
    shariter->resetValue (false);
    sensisha->resetValue (false);
    inverssha->set_active (false);
    noiselumf->resetValue (false);
    noiselumc->resetValue (false);
    noisechrof->resetValue (false);
    noisechroc->resetValue (false);


}


void Locallab::lumaneutralPressed ()
{

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue (100);
        adjusterChanged (multiplier[i], 100);
    }
}


void Locallab::lumacontrastPlusPressed ()
{

    for (int i = 0; i < 5; i++) {
        float inc = (5 - i);
        multiplier[i]->setValue (multiplier[i]->getValue() + inc);
        adjusterChanged (multiplier[i], multiplier[i]->getValue());
    }
}


void Locallab::lumacontrastMinusPressed ()
{

    for (int i = 0; i < 5; i++) {
        float inc = - (5 - i);
        multiplier[i]->setValue (multiplier[i]->getValue() + inc);
        adjusterChanged (multiplier[i], multiplier[i]->getValue());
    }
}



void Locallab::autoOpenCurve ()
{
    cTgainshape->openIfNonlinear();
    //  llshape->openIfNonlinear();
    //  LHshape->openIfNonlinear();

}


int localChangedUI (void* data)
{

    GThreadLock lock;
    (static_cast<Locallab*> (data))->localComputed_ ();

    return 0;
}

int localretChangedUI (void* data)
{

    GThreadLock lock;
    (static_cast<Locallab*> (data))->localretComputed_ ();

    return 0;
}

bool Locallab::localretComputed_ ()
{
    disableListener ();

    //Reticurv
//update GUI and MIP specially for curve

    int *s_datc;
    s_datc = new int[70];
    int siz;
    //printf("nexts=%s\n", nextstr2.c_str());
    ImProcFunctions::strcurv_data (nextstr2, s_datc, siz);
    std::vector<double>   creti;

    for (int j = 0; j < siz; j++) {
        creti.push_back ((double) (s_datc[j]) / 1000.);
    }

    delete [] s_datc;

    cTgainshape->setCurve (creti);

    int *s_datcl;
    s_datcl = new int[70];
    int sizl;
    ImProcFunctions::strcurv_data (nextll_str2, s_datcl, sizl);
    std::vector<double>   cll;

    for (int j = 0; j < sizl; j++) {
        cll.push_back ((double) (s_datcl[j]) / 1000.);
    }

    delete [] s_datcl;

    llshape->setCurve (cll);


    int *s_datcc;
    s_datcc = new int[70];
    int sizc;
    ImProcFunctions::strcurv_data (nextcc_str2, s_datcc, sizc);
    std::vector<double>   ccc;

    for (int j = 0; j < sizc; j++) {
        ccc.push_back ((double) (s_datcc[j]) / 1000.);
    }

    delete [] s_datcc;

    ccshape->setCurve (ccc);

    int *s_datch;
    s_datch = new int[70];
    int sizh;
    ImProcFunctions::strcurv_data (nextlh_str2, s_datch, sizh);
    std::vector<double>   clh;

    for (int j = 0; j < sizh; j++) {
        clh.push_back ((double) (s_datch[j]) / 1000.);
    }

    delete [] s_datch;

    LHshape->setCurve (clh);


    enableListener ();

    //update all sliders by this strange process!
    if (anbspot->getValue() == 0) {
        anbspot->setValue (1);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged (anbspot, 1);

    } else if (anbspot->getValue() == 1) {
        anbspot->setValue (0);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged (anbspot, 0);

    }

    //update all curves
    std::vector<double>   cretirab;
    cretirab = cTgainshaperab->getCurve ();

    if (cretirab.at (5) == 0.70) {
        cretirab.at (5) = 0.9;
        cTgainshaperab->setCurve (cretirab);

        curveChanged (cTgainshaperab);
    } else if (cretirab.at (5) == 0.90) {
        cretirab.at (5) = 0.7;
        cTgainshaperab->setCurve (cretirab);
        curveChanged (cTgainshaperab);

    }


//    printf("G2 anbspot=%i\n", anbspot->getValue());

    if (listener) { //for all sliders
        listener->panelChanged (Evlocallabanbspot, "");//anbspot->getTextValue());
    }

    if (listener) {//for curve
        listener->panelChanged (EvlocallabCTgainCurverab, M (""));
    }

    if (listener) {//for curve
        listener->panelChanged (EvlocallabCTgainCurve, M (""));
    }

    if (listener) {//for curve
        listener->panelChanged (Evlocallabllshape, M (""));
    }

    if (listener) {//for curve
        listener->panelChanged (Evlocallabccshape, M (""));
    }

    if (listener) {//for curve
        listener->panelChanged (EvlocallabLHshape, M (""));
    }


}

bool Locallab::localComputed_ ()
{
//update GUI and MIP
    disableListener ();

    //size spot
    circrad->setValue (nextdatasp[2]);
    //center and cursor
    locX->setValue (nextdatasp[3]);
    locY->setValue (nextdatasp[4]);
    locYT->setValue (nextdatasp[5]);
    locXL->setValue (nextdatasp[6]);
    centerX->setValue (nextdatasp[7]);
    centerY->setValue (nextdatasp[8]);

    //sliders
    lightness->setValue (nextdatasp[9]);
    contrast->setValue (nextdatasp[10]);
    chroma->setValue (nextdatasp[11]);
    sensi->setValue (nextdatasp[12]);
    transit->setValue (nextdatasp[13]);

    //inverse
    if (nextdatasp[14] == 0) {
        invers->set_active (false);
    } else {
        invers->set_active (true);
    }

    //method cursor
    if (nextdatasp[15] == 0) {
        Smethod->set_active (0);
    } else if (nextdatasp[15] == 1) {
        Smethod->set_active (1);
    } else if (nextdatasp[15] == 2) {
        Smethod->set_active (2);
    } else if (nextdatasp[15] == 3) {
        Smethod->set_active (3);
    }

    //sliders blurr
    radius->setValue (nextdatasp[17]);
    strength->setValue (nextdatasp[18]);
    sensibn->setValue (nextdatasp[19]);

    //inverse
    if (nextdatasp[20] == 0) {
        inversrad->set_active (false);
    } else {
        inversrad->set_active (true);
    }

    //sliders retinex
    str->setValue (nextdatasp[21]);
    chrrt->setValue (nextdatasp[22]);
    neigh->setValue (nextdatasp[23]);
    vart->setValue (nextdatasp[24]);
    sensih->setValue (nextdatasp[25]);

    //inverse
    if (nextdatasp[26] == 0) {
        inversret->set_active (false);
    } else {
        inversret->set_active (true);
    }

    //method retinex
    if (nextdatasp[27] == 0) {
        retinexMethod->set_active (0);
    } else if (nextdatasp[27] == 1) {
        retinexMethod->set_active (1);
    } else if (nextdatasp[27] == 2) {
        retinexMethod->set_active (2);
    }

    //sharpening
    sharradius->setValue (nextdatasp[28]);
    sharamount->setValue (nextdatasp[29]);
    shardamping->setValue (nextdatasp[30]);
    shariter->setValue (nextdatasp[31]);
    sensisha->setValue (nextdatasp[32]);

    if (nextdatasp[33] == 0) {
        inverssha->set_active (false);
    } else {
        inverssha->set_active (true);
    }

    if (nextdatasp[34] == 0) {
        qualityMethod->set_active (0);
    } else if (nextdatasp[34] == 1) {
        qualityMethod->set_active (1);
    } else if (nextdatasp[34] == 2) {
        qualityMethod->set_active (2);
    }

    thres->setValue (nextdatasp[35]);
    proxi->setValue (nextdatasp[36]);

    //denoise
    noiselumf->setValue (nextdatasp[37]);
    noiselumc->setValue (nextdatasp[38]);
    noisechrof->setValue (nextdatasp[39]);
    noisechroc->setValue (nextdatasp[40]);

    //cbdl
    multiplier[0]->setValue (nextdatasp[41]);
    multiplier[1]->setValue (nextdatasp[42]);
    multiplier[2]->setValue (nextdatasp[43]);
    multiplier[3]->setValue (nextdatasp[44]);
    multiplier[4]->setValue (nextdatasp[45]);
    threshold->setValue (nextdatasp[46]);
    sensicb->setValue (nextdatasp[47]);

    //blur luma
    if (nextdatasp[48] == 0) {
        activlum->set_active (false);
    } else {
        activlum->set_active (true);
    }

//TM
    stren->setValue (nextdatasp[49]);
    gamma->setValue (nextdatasp[50]);
    estop->setValue (nextdatasp[51]);
    scaltm->setValue (nextdatasp[52]);
    rewei->setValue (nextdatasp[53]);
    sensitm->setValue (nextdatasp[54]);
    //  usleep(10000);

    //Reticurv
    retrab->setValue (nextdatasp[55]);

    //curvactiv
    if (nextdatasp[56] == 0) {
        curvactiv->set_active (false);
    } else {
        curvactiv->set_active (true);
    }

    if (nextdatasp[57] == 0) {
        qualitycurveMethod->set_active (0);
    } else if (nextdatasp[57] == 1) {
        qualitycurveMethod->set_active (1);
    } else if (nextdatasp[57] == 2) {
        qualitycurveMethod->set_active (2);
    }

    double intermed = 0.01 * (double) nextdatasp[58];
    hueref->setValue (intermed);
    chromaref->setValue (nextdatasp[59]);
    lumaref->setValue (nextdatasp[60]);

    int *s_datc;
    s_datc = new int[70];
    int siz;
    ImProcFunctions::strcurv_data (nextstr, s_datc, siz);


    std::vector<double>   creti;

    for (int j = 0; j < siz; j++) {
        creti.push_back ((double) (s_datc[j]) / 1000.);
    }

    delete [] s_datc;

    cTgainshape->setCurve (creti);

    //LLcurv
    int *s_datcl;
    s_datcl = new int[70];
    int sizl;
    ImProcFunctions::strcurv_data (nextll_str, s_datcl, sizl);


    std::vector<double>   cll;

    for (int j = 0; j < sizl; j++) {
        cll.push_back ((double) (s_datcl[j]) / 1000.);
    }

    delete [] s_datcl;
    llshape->setCurve (cll);

    //CCcurv
    int *s_datcc;
    s_datcc = new int[70];
    int sizc;
    ImProcFunctions::strcurv_data (nextcc_str, s_datcc, sizc);


    std::vector<double>   ccc;

    for (int j = 0; j < sizc; j++) {
        ccc.push_back ((double) (s_datcc[j]) / 1000.);
    }

    delete [] s_datcc;
    ccshape->setCurve (ccc);


    //LHcurv
    int *s_datch;
    s_datch = new int[70];
    int sizh;
    ImProcFunctions::strcurv_data (nextlh_str, s_datch, sizh);


    std::vector<double>   clh;

    for (int j = 0; j < sizh; j++) {
        clh.push_back ((double) (s_datch[j]) / 1000.);
    }

    delete [] s_datch;
    LHshape->setCurve (clh);


    //  usleep(10000);


    enableListener ();

    //update all sliders by this strange process!
    if (anbspot->getValue() == 0) {
        anbspot->setValue (1);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged (anbspot, 1);

    } else if (anbspot->getValue() == 1) {
        anbspot->setValue (0);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged (anbspot, 0);

    }


    //update all curves
    std::vector<double>   cretirab;
    cretirab = cTgainshaperab->getCurve ();

    if (cretirab.at (5) == 0.70) {
        cretirab.at (5) = 0.9;
        cTgainshaperab->setCurve (cretirab);

        curveChanged (cTgainshaperab);
    } else if (cretirab.at (5) == 0.90) {
        cretirab.at (5) = 0.7;
        cTgainshaperab->setCurve (cretirab);
        curveChanged (cTgainshaperab);

    }

    //

//   printf("G1 maj anbspot=%i  cretirab=%f\n", anbspot->getValue(), cretirab.at(5));


    //add events for each cases
    if (listener) { //for all sliders
        listener->panelChanged (Evlocallabanbspot, "");//anbspot->getTextValue());
    }

    if (listener) {//for curve
        listener->panelChanged (EvlocallabCTgainCurverab, M (""));
    }

    if (listener) {//for inverse color
        listener->panelChanged (Evlocallabinvers, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for curvactiv
        listener->panelChanged (Evlocallabcurvactiv, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for inverse blurr
        listener->panelChanged (Evlocallabinversrad, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for quality method
        listener->panelChanged (EvlocallabqualityMethod, qualityMethod->get_active_text ());
    }

    if (listener) {//for quality method
        listener->panelChanged (EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text ());
    }

    if (listener) {//for inverse retinex
        listener->panelChanged (Evlocallabinversret, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for inverse sharpen
        listener->panelChanged (Evlocallabinverssha, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for Smethod : position of mouse cursor
        listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
    }

    if (listener) {//for retinex method
        listener->panelChanged (EvlocallabretinexMethod, retinexMethod->get_active_text ());
    }

    if (listener) {//for curve reti
        listener->panelChanged (EvlocallabCTgainCurve, M (""));
    }

    if (listener) {//for curve LL
        listener->panelChanged (Evlocallabllshape, M (""));
    }

    if (listener) {//for curve LH
        listener->panelChanged (EvlocallabLHshape, M (""));
    }

    if (listener) {//for curve LH
        listener->panelChanged (Evlocallabccshape, M (""));
    }

    return false;
}

void Locallab::localChanged  (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, int sp, int maxdat)
{
    for (int i = 2; i < 61; i++) {
        nextdatasp[i] = datasp[i][sp];
    }

    nextstr = datastr;
    nextll_str = ll_str;
    nextlh_str = lh_str;
    nextcc_str = cc_str;

    nextlength = maxdat;
    g_idle_add (localChangedUI, this);
}

void Locallab::localretChanged  (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, int sp, int maxdat)
{
    nextlength = maxdat;
    nextstr2 = datastr;
    nextll_str2 = ll_str;
    nextlh_str2 = lh_str;
    nextcc_str2 = cc_str;

    g_idle_add (localretChangedUI, this);
}


void Locallab::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    enablecolorConn.block (true);
    enableblurConn.block (true);
    enabletonemapConn.block (true);
    enableretiConn.block (true);
    enablesharpConn.block (true);
    enablecbdlConn.block (true);
    enabledenoiConn.block (true);


    if (pedited) {
        degree->setEditedState (pedited->locallab.degree ? Edited : UnEdited);
        locY->setEditedState (pedited->locallab.locY ? Edited : UnEdited);
        locX->setEditedState (pedited->locallab.locX ? Edited : UnEdited);
        locYT->setEditedState (pedited->locallab.locYT ? Edited : UnEdited);
        locXL->setEditedState (pedited->locallab.locXL ? Edited : UnEdited);
        centerX->setEditedState (pedited->locallab.centerX ? Edited : UnEdited);
        centerY->setEditedState (pedited->locallab.centerY ? Edited : UnEdited);
        circrad->setEditedState (pedited->locallab.circrad ? Edited : UnEdited);
        thres->setEditedState (pedited->locallab.thres ? Edited : UnEdited);
        proxi->setEditedState (pedited->locallab.proxi ? Edited : UnEdited);
        lightness->setEditedState (pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setEditedState (pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setEditedState (pedited->locallab.chroma ? Edited : UnEdited);
        sharradius->setEditedState (pedited->locallab.sharradius ? Edited : UnEdited);
        sharamount->setEditedState (pedited->locallab.sharamount ? Edited : UnEdited);
        shardamping->setEditedState (pedited->locallab.shardamping ? Edited : UnEdited);
        shariter->setEditedState (pedited->locallab.shariter ? Edited : UnEdited);
        sensisha->setEditedState (pedited->locallab.sensisha ? Edited : UnEdited);
        noiselumf->setEditedState (pedited->locallab.noiselumf ? Edited : UnEdited);
        noiselumc->setEditedState (pedited->locallab.noiselumc ? Edited : UnEdited);
        noisechrof->setEditedState (pedited->locallab.noisechrof ? Edited : UnEdited);
        noisechroc->setEditedState (pedited->locallab.noisechroc ? Edited : UnEdited);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setEditedState (pedited->locallab.mult[i] ? Edited : UnEdited);
        }

        threshold->setEditedState (pedited->locallab.threshold ? Edited : UnEdited);

        sensi->setEditedState (pedited->locallab.sensi ? Edited : UnEdited);
        sensih->setEditedState (pedited->locallab.sensih ? Edited : UnEdited);
        retrab->setEditedState (pedited->locallab.retrab ? Edited : UnEdited);
        sensicb->setEditedState (pedited->locallab.sensicb ? Edited : UnEdited);
        sensibn->setEditedState (pedited->locallab.sensibn ? Edited : UnEdited);
        sensitm->setEditedState (pedited->locallab.sensitm ? Edited : UnEdited);
        radius->setEditedState (pedited->locallab.radius ? Edited : UnEdited);
        strength->setEditedState (pedited->locallab.strength ? Edited : UnEdited);
        stren->setEditedState (pedited->locallab.stren ? Edited : UnEdited);
        gamma->setEditedState (pedited->locallab.gamma ? Edited : UnEdited);
        estop->setEditedState (pedited->locallab.estop ? Edited : UnEdited);
        scaltm->setEditedState (pedited->locallab.scaltm ? Edited : UnEdited);
        rewei->setEditedState (pedited->locallab.rewei ? Edited : UnEdited);

        nbspot->setEditedState (pedited->locallab.nbspot ? Edited : UnEdited);
        anbspot->setEditedState (pedited->locallab.anbspot ? Edited : UnEdited);
        hueref->setEditedState (pedited->locallab.hueref ? Edited : UnEdited);
        chromaref->setEditedState (pedited->locallab.chromaref ? Edited : UnEdited);
        lumaref->setEditedState (pedited->locallab.lumaref ? Edited : UnEdited);
        transit->setEditedState (pedited->locallab.transit ? Edited : UnEdited);
        str->setEditedState (pedited->locallab.str ? Edited : UnEdited);
        neigh->setEditedState (pedited->locallab.neigh ? Edited : UnEdited);
        vart->setEditedState (pedited->locallab.vart ? Edited : UnEdited);
        chrrt->setEditedState (pedited->locallab.chrrt ? Edited : UnEdited);
        set_inconsistent (multiImage && !pedited->locallab.enabled);
        avoid->set_inconsistent (multiImage && !pedited->locallab.avoid);
        activlum->set_inconsistent (multiImage && !pedited->locallab.activlum);
        invers->set_inconsistent (multiImage && !pedited->locallab.invers);
        curvactiv->set_inconsistent (multiImage && !pedited->locallab.curvactiv);
        inversrad->set_inconsistent (multiImage && !pedited->locallab.inversrad);
        inverssha->set_inconsistent (multiImage && !pedited->locallab.inverssha);
        cTgainshape->setUnChanged  (!pedited->locallab.localTgaincurve);
        llshape->setUnChanged  (!pedited->locallab.llcurve);
        ccshape->setUnChanged  (!pedited->locallab.cccurve);
        LHshape->setUnChanged  (!pedited->locallab.LHcurve);
        inversret->set_inconsistent (multiImage && !pedited->locallab.inversret);
        cTgainshaperab->setUnChanged  (!pedited->locallab.localTgaincurverab);
        expcolor->set_inconsistent   (!pedited->locallab.expcolor);
        expblur->set_inconsistent   (!pedited->locallab.expblur);
        exptonemap->set_inconsistent   (!pedited->locallab.exptonemap);
        expreti->set_inconsistent   (!pedited->locallab.expreti);
        expsharp->set_inconsistent   (!pedited->locallab.expsharp);
        expcbdl->set_inconsistent   (!pedited->locallab.expcbdl);
        expdenoi->set_inconsistent   (!pedited->locallab.expdenoi);

        if (!pedited->locallab.Smethod) {
            Smethod->set_active_text (M ("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.retinexMethod) {
            retinexMethod->set_active_text (M ("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.qualityMethod) {
            qualityMethod->set_active_text (M ("GENERAL_UNCHANGED"));
        }

        if (!pedited->locallab.qualitycurveMethod) {
            qualitycurveMethod->set_active_text (M ("GENERAL_UNCHANGED"));
        }

    }

    setEnabled (pp->locallab.enabled);

    Smethodconn.block (true);
    retinexMethodConn.block (true);
    qualityMethodConn.block (true);
    qualitycurveMethodConn.block (true);

    avoidConn.block (true);
    avoid->set_active (pp->locallab.avoid);
    avoidConn.block (false);
    activlumConn.block (true);
    activlum->set_active (pp->locallab.activlum);
    activlumConn.block (false);
    inversConn.block (true);
    invers->set_active (pp->locallab.invers);
    inversConn.block (false);
    curvactivConn.block (true);
    curvactiv->set_active (pp->locallab.curvactiv);
    curvactivConn.block (false);
    inversradConn.block (true);
    inversrad->set_active (pp->locallab.inversrad);
    inversradConn.block (false);
    inversretConn.block (true);
    inversret->set_active (pp->locallab.inversret);
    inversretConn.block (false);
    inversshaConn.block (true);
    inverssha->set_active (pp->locallab.inverssha);
    inversshaConn.block (false);

    degree->setValue (pp->locallab.degree);
    locY->setValue (pp->locallab.locY);
    locX->setValue (pp->locallab.locX);
    locYT->setValue (pp->locallab.locYT);
    locXL->setValue (pp->locallab.locXL);
    centerX->setValue (pp->locallab.centerX);
    centerY->setValue (pp->locallab.centerY);
    circrad->setValue (pp->locallab.circrad);
    thres->setValue (pp->locallab.thres);
    proxi->setValue (pp->locallab.proxi);
    lightness->setValue (pp->locallab.lightness);
    contrast->setValue (pp->locallab.contrast);
    chroma->setValue (pp->locallab.chroma);
    sharradius->setValue (pp->locallab.sharradius);
    sharamount->setValue (pp->locallab.sharamount);
    shardamping->setValue (pp->locallab.shardamping);
    shariter->setValue (pp->locallab.shariter);
    sensisha->setValue (pp->locallab.sensisha);
    sensi->setValue (pp->locallab.sensi);
    sensih->setValue (pp->locallab.sensih);
    retrab->setValue (pp->locallab.retrab);
    sensicb->setValue (pp->locallab.sensicb);
    sensibn->setValue (pp->locallab.sensibn);
    sensitm->setValue (pp->locallab.sensitm);
    transit->setValue (pp->locallab.transit);
    radius->setValue (pp->locallab.radius);
    strength->setValue (pp->locallab.strength);
    stren->setValue (pp->locallab.stren);
    gamma->setValue (pp->locallab.gamma);
    estop->setValue (pp->locallab.estop);
    scaltm->setValue (pp->locallab.scaltm);
    rewei->setValue (pp->locallab.rewei);
    str->setValue (pp->locallab.str);
    neigh->setValue (pp->locallab.neigh);
    nbspot->setValue (pp->locallab.nbspot);
    anbspot->setValue (pp->locallab.anbspot);
    hueref->setValue (pp->locallab.hueref);
    chromaref->setValue (pp->locallab.chromaref);
    lumaref->setValue (pp->locallab.lumaref);
    vart->setValue (pp->locallab.vart);
    chrrt->setValue (pp->locallab.chrrt);
    cTgainshape->setCurve (pp->locallab.localTgaincurve);
    cTgainshaperab->setCurve (pp->locallab.localTgaincurverab);
    llshape->setCurve (pp->locallab.llcurve);
    ccshape->setCurve (pp->locallab.cccurve);
    LHshape->setCurve (pp->locallab.LHcurve);
    lastactivlum = pp->locallab.activlum;
    lastanbspot = pp->locallab.anbspot;
    noiselumf->setValue (pp->locallab.noiselumf);
    noiselumc->setValue (pp->locallab.noiselumc);
    noisechrof->setValue (pp->locallab.noisechrof);
    noisechroc->setValue (pp->locallab.noisechroc);
    expcolor->setEnabled (pp->locallab.expcolor);
    expblur->setEnabled (pp->locallab.expblur);
    exptonemap->setEnabled (pp->locallab.exptonemap);
    expreti->setEnabled (pp->locallab.expreti);
    expsharp->setEnabled (pp->locallab.expsharp);
    expcbdl->setEnabled (pp->locallab.expcbdl);
    expdenoi->setEnabled (pp->locallab.expdenoi);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue (pp->locallab.mult[i]);
    }

    threshold->setValue (pp->locallab.threshold);

    lastavoid = pp->locallab.avoid;
    lastinvers = pp->locallab.invers;
    lastcurvactiv = pp->locallab.curvactiv;
    lastinversrad = pp->locallab.inversrad;
    lastinversret = pp->locallab.inversret;
    lastinverssha = pp->locallab.inverssha;
    activlumChanged();
    inversChanged();
    curvactivChanged();
    inversradChanged();
    inversretChanged();
    inversshaChanged();

    updateGeometry (pp->locallab.centerX, pp->locallab.centerY, pp->locallab.circrad, pp->locallab.locY, pp->locallab.degree,  pp->locallab.locX, pp->locallab.locYT, pp->locallab.locXL);

    if (pp->locallab.Smethod == "IND") {
        Smethod->set_active (0);
    } else if (pp->locallab.Smethod == "SYM") {
        Smethod->set_active (1);
    } else if (pp->locallab.Smethod == "INDSL") {
        Smethod->set_active (2);
    } else if (pp->locallab.Smethod == "SYMSL") {
        Smethod->set_active (3);
    }

    SmethodChanged();
    Smethodconn.block (false);

    if (pp->locallab.retinexMethod == "low") {
        retinexMethod->set_active (0);
    } else if (pp->locallab.retinexMethod == "uni") {
        retinexMethod->set_active (1);
    } else if (pp->locallab.retinexMethod == "high") {
        retinexMethod->set_active (2);
    }

    retinexMethodChanged ();
    retinexMethodConn.block (false);

    if (pp->locallab.qualityMethod == "std") {
        qualityMethod->set_active (0);
    } else if (pp->locallab.qualityMethod == "enh") {
        qualityMethod->set_active (1);
    } else if (pp->locallab.qualityMethod == "enhden") {
        qualityMethod->set_active (2);
    }

    qualityMethodChanged ();
    qualityMethodConn.block (false);

    if (pp->locallab.qualitycurveMethod == "none") {
        qualitycurveMethod->set_active (0);
    } else if (pp->locallab.qualitycurveMethod == "std") {
        qualitycurveMethod->set_active (1);
    } else if (pp->locallab.qualitycurveMethod == "enh") {
        qualitycurveMethod->set_active (2);
    }

    qualitycurveMethodChanged ();
    qualitycurveMethodConn.block (false);

    anbspot->hide();
    hueref->hide();
    chromaref->hide();
    lumaref->hide();

    if (pp->locallab.Smethod == "SYM" || pp->locallab.Smethod == "SYMSL") {
        locXL->setValue (locX->getValue());
        locYT->setValue (locY->getValue());
    } else if (pp->locallab.Smethod == "LOC") {
        locXL->setValue (locX->getValue());
        locYT->setValue (locX->getValue());
        locY->setValue (locX->getValue());
    } else if (pp->locallab.Smethod == "INDSL" || pp->locallab.Smethod == "IND") {
        locX->setValue (pp->locallab.locX);
        locY->setValue (pp->locallab.locY);
        locXL->setValue (pp->locallab.locXL);
        locYT->setValue (pp->locallab.locYT);

    }

    enablecolorConn.block (false);
    enableblurConn.block (false);
    enabletonemapConn.block (false);
    enableretiConn.block (false);
    enablesharpConn.block (false);
    enablecbdlConn.block (false);
    enabledenoiConn.block (false);

    enableListener ();
}

void Locallab::updateGeometry (const int centerX_, const int centerY_, const int circrad_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth, const int fullHeight)
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
        dataProvider->getImageSize (imW, imH);

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
    rtengine::Coord origin (imW / 2 + centerX_ * imW / 2000.f, imH / 2 + centerY_ * imH / 2000.f);
//   printf("deX=%f dexL=%f deY=%f deyT=%f locX=%i locY=%i\n", decayX, decayXL, decayY, decayYT, locX_, locY_);

    if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
        decayYT = decayY;
        decayXL = decayX;
    }

    Line *currLine;
    Circle *currCircle;
    //  Arcellipse *currArcellipse;
    Beziers *currBeziers;
    double decay;
    const auto updateLine = [&] (Geometry * geometry, const float radius, const float begin, const float end) {
        const auto line = static_cast<Line*> (geometry);
        line->begin = PolarCoord (radius, -degree_ + begin);
        line->begin += origin;
        line->end = PolarCoord (radius, -degree_ + end);
        line->end += origin;
    };

    const auto updateLineWithDecay = [&] (Geometry * geometry, const float radius, const float decal, const float offSetAngle) {
        const auto line = static_cast<Line*> (geometry); //180
        line->begin = PolarCoord (radius, -degree_ + decal) + PolarCoord (decay, -degree_ + offSetAngle);
        line->begin += origin;//0
        line->end = PolarCoord (radius, -degree_ + (decal - 180)) + PolarCoord (decay, -degree_ + offSetAngle);
        line->end += origin;
    };

    const auto updateCircle = [&] (Geometry * geometry) {
        const auto circle = static_cast<Circle*> (geometry);
        circle->center = origin;
        circle->radius = circrad_;
    };

    const auto updateBeziers = [&] (Geometry * geometry, const double dX_, const double dI_, const double dY_,  const float begi, const float inte, const float en) {
        const auto beziers = static_cast<Beziers*> (geometry);
        beziers->begin = PolarCoord (dX_, begi);
        beziers->begin += origin;//0
        beziers->inter = PolarCoord (dI_, inte);
        beziers->inter += origin;//0
        beziers->end = PolarCoord (dY_,  en);
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
    updateLineWithDecay (visibleGeometry.at (0), dimline, 90., 0.);
    updateLineWithDecay (mouseOverGeometry.at (0), dimline, 90., 0.);

    decay = decayXL;

    updateLineWithDecay (visibleGeometry.at (1), dimline, 90., 180.);
    updateLineWithDecay (mouseOverGeometry.at (1), dimline, 90., 180.);

    decay = decayYT;
    updateLineWithDecay (visibleGeometry.at (2), dimline, 180., 270.);
    updateLineWithDecay (mouseOverGeometry.at (2), dimline, 180., 270.);

    decay = decayY;

    updateLineWithDecay (visibleGeometry.at (3), dimline, 180, 90.);
    updateLineWithDecay (mouseOverGeometry.at (3), dimline, 180., 90.);


    updateCircle (visibleGeometry.at (4));
    updateCircle (mouseOverGeometry.at (4));

    if (options.showdelimspot) {
        //this decayww evaluate approximation of a point in the ellipse for an angle alpha
        double decay15 = 1.07854 * ((decayX * decayY) / sqrt (0.07179 * SQR (decayX) + SQR (decayY))); //0.07179 = SQR(sin(15)/cos(15))  1.07854 = 1 / cos(15)
        double decay30 = 1.15473 * ((decayX * decayY) / sqrt (0.33335 * SQR (decayX) + SQR (decayY)));
        double decay60 = 2. * ((decayX * decayY) / sqrt (3.0 * SQR (decayX) + SQR (decayY)));
        double decay75 = 3.86398 * ((decayX * decayY) / sqrt (13.929 * SQR (decayX) + SQR (decayY)));

        double decay15L = 1.07854 * ((decayXL * decayY) / sqrt (0.07179 * SQR (decayXL) + SQR (decayY)));
        double decay30L = 1.15473 * ((decayXL * decayY) / sqrt (0.33335 * SQR (decayXL) + SQR (decayY)));
        double decay60L = 2. * ((decayXL * decayY) / sqrt (3.0 * SQR (decayXL) + SQR (decayY)));
        double decay75L = 3.86398 * ((decayXL * decayY) / sqrt (13.929 * SQR (decayXL) + SQR (decayY)));

        double decay15LT = 1.07854 * ((decayXL * decayYT) / sqrt (0.07179 * SQR (decayXL) + SQR (decayYT)));
        double decay30LT = 1.15473 * ((decayXL * decayYT) / sqrt (0.33335 * SQR (decayXL) + SQR (decayYT)));
        double decay60LT = 2. * ((decayXL * decayYT) / sqrt (3.0 * SQR (decayXL) + SQR (decayYT)));
        double decay75LT = 3.86398 * ((decayXL * decayYT) / sqrt (13.929 * SQR (decayXL) + SQR (decayYT)));

        double decay15T = 1.07854 * ((decayX * decayYT) / sqrt (0.07179 * SQR (decayX) + SQR (decayYT)));
        double decay30T = 1.15473 * ((decayX * decayYT) / sqrt (0.33335 * SQR (decayX) + SQR (decayYT)));
        double decay60T = 2. * ((decayX * decayYT) / sqrt (3.0 * SQR (decayX) + SQR (decayYT)));
        double decay75T = 3.86398 * ((decayX * decayYT) / sqrt (13.929 * SQR (decayX) + SQR (decayYT)));

        double decay45 = (1.414 * decayX * decayY) / sqrt (SQR (decayX) + SQR (decayY));
        double decay45L = (1.414 * decayXL * decayY) / sqrt (SQR (decayXL) + SQR (decayY));
        double decay45LT = (1.414 * decayXL * decayYT) / sqrt (SQR (decayXL) + SQR (decayYT));
        double decay45T = (1.414 * decayX * decayYT) / sqrt (SQR (decayX) + SQR (decayYT));

        //printf("decayX=%f decayY=%f decay10=%f decay45=%f oriX=%i origY=%i\n", decayX, decayY, decay10, decay45, origin.x, origin.y);
        updateBeziers (visibleGeometry.at (5), decayX, decay15  , decay30, 0., 15., 30.);
        updateBeziers (mouseOverGeometry.at (5), decayX, decay15 , decay30, 0., 15., 30.);

        updateBeziers (visibleGeometry.at (6), decay30, decay45 , decay60, 30., 45., 60.);
        updateBeziers (mouseOverGeometry.at (6), decay30, decay45 , decay60, 30., 45., 60.);

        updateBeziers (visibleGeometry.at (7), decay60, decay75 , decayY, 60., 75., 90.);
        updateBeziers (mouseOverGeometry.at (7), decay60, decay75 , decayY, 60., 75., 90.);

        updateBeziers (visibleGeometry.at (8), decayY, decay75L  , decay60L, 90., 105., 120.);
        updateBeziers (mouseOverGeometry.at (8), decayY, decay75L , decay60L, 90., 105., 120.);

        updateBeziers (visibleGeometry.at (9), decay60L, decay45L  , decay30L, 120., 135., 150.);
        updateBeziers (mouseOverGeometry.at (9), decay60L, decay45L , decay30L, 120., 135., 150.);

        updateBeziers (visibleGeometry.at (10), decay30L, decay15L  , decayXL, 150., 165., 180.);
        updateBeziers (mouseOverGeometry.at (10), decay30L, decay15L , decayXL, 150., 165., 180.);

        updateBeziers (visibleGeometry.at (11), decayXL, decay15LT  , decay30LT, 180., 195., 210.);
        updateBeziers (mouseOverGeometry.at (11), decayXL, decay15LT , decay30LT, 180., 195., 210.);

        updateBeziers (visibleGeometry.at (12), decay30LT, decay45LT  , decay60LT, 210., 225., 240.);
        updateBeziers (mouseOverGeometry.at (12), decay30LT, decay45LT , decay60LT, 210., 225., 240.);

        updateBeziers (visibleGeometry.at (13), decay60LT, decay75LT  , decayYT, 240., 255., 270.);
        updateBeziers (mouseOverGeometry.at (13), decay60LT, decay75LT , decayYT, 240., 255., 270.);

        updateBeziers (visibleGeometry.at (14), decayYT, decay75T  , decay60T, 270., 285., 300.);
        updateBeziers (mouseOverGeometry.at (14), decayYT, decay75T , decay60T, 270., 285., 300.);

        updateBeziers (visibleGeometry.at (15), decay60T, decay45T  , decay30T, 300., 315., 330.);
        updateBeziers (mouseOverGeometry.at (15), decay60T, decay45T , decay30T, 300., 315., 330.);

        updateBeziers (visibleGeometry.at (16), decay30T, decay15T  , decayX, 330., 345., 360.);
        updateBeziers (mouseOverGeometry.at (16), decay30T, decay15T , decayX, 330., 345., 360.);

    }

    //  updateArcellipse (visibleGeometry.at (5), decayX, decayY, 0., 0.5);
    //  updateArcellipse (mouseOverGeometry.at (5), decayX, decayY, 0., 0.5);

}

void Locallab::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->locallab.degree = degree->getValue ();
    pp->locallab.locY = locY->getIntValue ();
    pp->locallab.locX = locX->getValue ();
    pp->locallab.locYT = locYT->getIntValue ();
    pp->locallab.locXL = locXL->getValue ();
    pp->locallab.centerX = centerX->getIntValue ();
    pp->locallab.centerY = centerY->getIntValue ();
    pp->locallab.circrad = circrad->getIntValue ();
    pp->locallab.proxi = proxi->getIntValue ();
    pp->locallab.thres = thres->getIntValue ();
    pp->locallab.lightness = lightness->getIntValue ();
    pp->locallab.contrast = contrast->getIntValue ();
    pp->locallab.chroma = chroma->getIntValue ();
    pp->locallab.noiselumc = noiselumc->getIntValue ();
    pp->locallab.noiselumf = noiselumf->getIntValue ();
    pp->locallab.noisechrof = noisechrof->getIntValue ();
    pp->locallab.noisechroc = noisechroc->getIntValue ();
    pp->locallab.sharradius = sharradius->getIntValue ();
    pp->locallab.sharamount = sharamount->getIntValue ();
    pp->locallab.shardamping = shardamping->getIntValue ();
    pp->locallab.shariter = shariter->getIntValue ();
    pp->locallab.sensisha = sensisha->getIntValue ();
    pp->locallab.sensi = sensi->getIntValue ();
    pp->locallab.sensih = sensih->getIntValue ();
    pp->locallab.retrab = retrab->getIntValue ();
    pp->locallab.sensicb = sensicb->getIntValue ();
    pp->locallab.sensibn = sensibn->getIntValue ();
    pp->locallab.sensitm = sensitm->getIntValue ();
    pp->locallab.radius = radius->getIntValue ();
    pp->locallab.strength = strength->getIntValue ();
    pp->locallab.stren = stren->getIntValue ();
    pp->locallab.gamma = gamma->getIntValue ();
    pp->locallab.estop = estop->getIntValue ();
    pp->locallab.scaltm = scaltm->getIntValue ();
    pp->locallab.rewei = rewei->getIntValue ();
    pp->locallab.enabled = getEnabled();
    pp->locallab.transit = transit->getIntValue ();
    pp->locallab.avoid = avoid->get_active();
    pp->locallab.activlum = activlum->get_active();
    pp->locallab.invers = invers->get_active();
    pp->locallab.curvactiv = curvactiv->get_active();
    pp->locallab.inversrad = inversrad->get_active();
    pp->locallab.inversret = inversret->get_active();
    pp->locallab.inverssha = inverssha->get_active();
    pp->locallab.str = str->getIntValue ();
    pp->locallab.neigh = neigh->getIntValue ();
    pp->locallab.nbspot = nbspot->getIntValue ();
    pp->locallab.anbspot = anbspot->getIntValue ();
    pp->locallab.hueref = hueref->getValue ();
    pp->locallab.chromaref = chromaref->getValue ();
    pp->locallab.lumaref = lumaref->getValue ();
    pp->locallab.vart = vart->getIntValue ();
    pp->locallab.chrrt = chrrt->getIntValue ();
    pp->locallab.localTgaincurve       = cTgainshape->getCurve ();
    pp->locallab.localTgaincurverab       = cTgainshaperab->getCurve ();
    pp->locallab.llcurve       = llshape->getCurve ();
    pp->locallab.cccurve       = ccshape->getCurve ();
    pp->locallab.LHcurve       = LHshape->getCurve ();
    pp->locallab.expcolor      = expcolor->getEnabled();
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

    if (pedited) {
        pedited->locallab.degree = degree->getEditedState ();
        pedited->locallab.Smethod  = Smethod->get_active_text() != M ("GENERAL_UNCHANGED");
        pedited->locallab.retinexMethod    = retinexMethod->get_active_text() != M ("GENERAL_UNCHANGED");
        pedited->locallab.qualityMethod    = qualityMethod->get_active_text() != M ("GENERAL_UNCHANGED");
        pedited->locallab.qualitycurveMethod    = qualitycurveMethod->get_active_text() != M ("GENERAL_UNCHANGED");
        pedited->locallab.locY = locY->getEditedState ();
        pedited->locallab.locX = locX->getEditedState ();
        pedited->locallab.locYT = locYT->getEditedState ();
        pedited->locallab.locXL = locXL->getEditedState ();
        pedited->locallab.centerX = centerX->getEditedState ();
        pedited->locallab.centerY = centerY->getEditedState ();
        pedited->locallab.circrad = circrad->getEditedState ();
        pedited->locallab.proxi = proxi->getEditedState ();
        pedited->locallab.thres = thres->getEditedState ();
        pedited->locallab.lightness = lightness->getEditedState ();
        pedited->locallab.contrast = contrast->getEditedState ();
        pedited->locallab.chroma = chroma->getEditedState ();
        pedited->locallab.noiselumf = noiselumf->getEditedState ();
        pedited->locallab.noiselumc = noiselumc->getEditedState ();
        pedited->locallab.noisechrof = noisechrof->getEditedState ();
        pedited->locallab.noisechroc = noisechroc->getEditedState ();
        pedited->locallab.sharradius = sharradius->getEditedState ();
        pedited->locallab.sharamount = sharamount->getEditedState ();
        pedited->locallab.shardamping = shardamping->getEditedState ();
        pedited->locallab.shariter = shariter->getEditedState ();
        pedited->locallab.sensisha = sensisha->getEditedState ();
        pedited->locallab.sensi = sensi->getEditedState ();
        pedited->locallab.sensih = sensih->getEditedState ();
        pedited->locallab.retrab = retrab->getEditedState ();
        pedited->locallab.sensicb = sensicb->getEditedState ();
        pedited->locallab.sensibn = sensibn->getEditedState ();
        pedited->locallab.sensitm = sensitm->getEditedState ();
        pedited->locallab.radius = radius->getEditedState ();
        pedited->locallab.strength = strength->getEditedState ();
        pedited->locallab.stren = stren->getEditedState ();
        pedited->locallab.gamma = gamma->getEditedState ();
        pedited->locallab.estop = estop->getEditedState ();
        pedited->locallab.scaltm = scaltm->getEditedState ();
        pedited->locallab.rewei = rewei->getEditedState ();
        pedited->locallab.transit = transit->getEditedState ();
        pedited->locallab.enabled = !get_inconsistent();
        pedited->locallab.avoid = !avoid->get_inconsistent();
        pedited->locallab.invers = !invers->get_inconsistent();
        pedited->locallab.curvactiv = !curvactiv->get_inconsistent();
        pedited->locallab.activlum = !activlum->get_inconsistent();
        pedited->locallab.inversret = !inversret->get_inconsistent();
        pedited->locallab.inversrad = !inversrad->get_inconsistent();
        pedited->locallab.inverssha = !inverssha->get_inconsistent();
        pedited->locallab.str = str->getEditedState ();
        pedited->locallab.neigh = neigh->getEditedState ();
        pedited->locallab.nbspot = nbspot->getEditedState ();
        pedited->locallab.anbspot = anbspot->getEditedState ();
        pedited->locallab.hueref = hueref->getEditedState ();
        pedited->locallab.chromaref = chromaref->getEditedState ();
        pedited->locallab.lumaref = lumaref->getEditedState ();

        pedited->locallab.vart = vart->getEditedState ();
        pedited->locallab.chrrt = chrrt->getEditedState ();
        pedited->locallab.localTgaincurve        = !cTgainshape->isUnChanged ();
        pedited->locallab.localTgaincurverab        = !cTgainshaperab->isUnChanged ();
        pedited->locallab.llcurve        = !llshape->isUnChanged ();
        pedited->locallab.cccurve        = !ccshape->isUnChanged ();
        pedited->locallab.LHcurve        = !LHshape->isUnChanged ();
        pedited->locallab.expcolor     = !expcolor->get_inconsistent();
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

    }

    if (retinexMethod->get_active_row_number() == 0) {
        pp->locallab.retinexMethod = "low";
    } else if (retinexMethod->get_active_row_number() == 1) {
        pp->locallab.retinexMethod = "uni";
    } else if (retinexMethod->get_active_row_number() == 2) {
        pp->locallab.retinexMethod = "high";
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

void Locallab::curveChanged (CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == cTgainshape) {
            listener->panelChanged (EvlocallabCTgainCurve, M ("HISTORY_CUSTOMCURVE"));//HISTORY_CUSTOMCURVE
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue (strval + 1);
            adjusterChanged (retrab, strval + 1);
            usleep (10000); //to test
            retrab->setValue (strval);

            adjusterChanged (retrab, strval);
        }

        else if (ce == cTgainshaperab) {
            listener->panelChanged (EvlocallabCTgainCurverab, M (""));
        } else if (ce == LHshape) {
            listener->panelChanged (EvlocallabLHshape, M (""));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue (strval + 1);
            adjusterChanged (retrab, strval + 1);
            usleep (10000); //to test
            retrab->setValue (strval);

            adjusterChanged (retrab, strval);


        } else if (ce == llshape) {
            listener->panelChanged (Evlocallabllshape, M ("HISTORY_CUSTOMCURVE"));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue (strval + 1);
            adjusterChanged (retrab, strval + 1);
            usleep (10000); //to test
            retrab->setValue (strval);

            adjusterChanged (retrab, strval);

        } else if (ce == ccshape) {
            listener->panelChanged (Evlocallabccshape, M ("HISTORY_CUSTOMCURVE"));
            int strval = retrab->getValue();
            //update MIP
            retrab->setValue (strval + 1);
            adjusterChanged (retrab, strval + 1);
            usleep (10000); //to test
            retrab->setValue (strval);

            adjusterChanged (retrab, strval);

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
        listener->panelChanged (EvlocallabretinexMethod, retinexMethod->get_active_text ());
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
        listener->panelChanged (EvlocallabqualityMethod, qualityMethod->get_active_text ());
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
        listener->panelChanged (EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text ());
    }
}

void Locallab::SmethodChanged ()
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
            listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
            locXL->setValue (locX->getValue());
            locYT->setValue (locY->getValue());
        }
        //   else if(Smethod->get_active_row_number()==2) {
        //          listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
        //           locXL->setValue (locX->getValue());
        //           locYT->setValue (locX->getValue());
        //          locY->setValue (locX->getValue());
        //     }
        else

        {
            listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());

        }
    }
}
void Locallab::inversChanged ()
{

    if (batchMode) {
        if (invers->get_inconsistent()) {
            invers->set_inconsistent (false);
            inversConn.block (true);
            invers->set_active (false);
            inversConn.block (false);
        } else if (lastinvers) {
            invers->set_inconsistent (true);
        }

        lastinvers = invers->get_active ();
    }

    if (invers->get_active ()) {
        sensi->hide();
        llCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        artifFrame->hide();
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
            listener->panelChanged (Evlocallabinvers, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabinvers, M ("GENERAL_DISABLED"));
        }
    }
}

void Locallab::curvactivChanged ()
{

    if (batchMode) {
        if (curvactiv->get_inconsistent()) {
            curvactiv->set_inconsistent (false);
            curvactivConn.block (true);
            curvactiv->set_active (false);
            curvactivConn.block (false);
        } else if (lastcurvactiv) {
            curvactiv->set_inconsistent (true);
        }

        lastcurvactiv = curvactiv->get_active ();
    }


    if (listener) {
        if (getEnabled()) {
            listener->panelChanged (Evlocallabcurvactiv, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabcurvactiv, M ("GENERAL_DISABLED"));
        }
    }
}

void Locallab::activlumChanged ()
{

    if (batchMode) {
        if (activlum->get_inconsistent()) {
            activlum->set_inconsistent (false);
            activlumConn.block (true);
            activlum->set_active (false);
            activlumConn.block (false);
        } else if (lastactivlum) {
            activlum->set_inconsistent (true);
        }

        lastactivlum = activlum->get_active ();
    }


    if (listener) {

        if (getEnabled()) {
            listener->panelChanged (Evlocallabactivlum, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabactivlum, M ("GENERAL_DISABLED"));
        }
    }
}

void Locallab::inversradChanged ()
{

    if (batchMode) {
        if (inversrad->get_inconsistent()) {
            inversrad->set_inconsistent (false);
            inversradConn.block (true);
            inversrad->set_active (false);
            inversradConn.block (false);
        } else if (lastinversrad) {
            inversrad->set_inconsistent (true);
        }

        lastinversrad = inversrad->get_active ();
    }

    if (inversrad->get_active ()) {
        sensibn->hide();
    } else {
        sensibn->show();
    }


    if (listener) {
        if (getEnabled()) {
            listener->panelChanged (Evlocallabinversrad, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabinversrad, M ("GENERAL_DISABLED"));
        }
    }
}

void Locallab::inversshaChanged ()
{

    if (batchMode) {
        if (inverssha->get_inconsistent()) {
            inverssha->set_inconsistent (false);
            inversshaConn.block (true);
            inverssha->set_active (false);
            inversshaConn.block (false);
        } else if (lastinverssha) {
            inverssha->set_inconsistent (true);
        }

        lastinverssha = inverssha->get_active ();
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
            listener->panelChanged (Evlocallabinverssha, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabinverssha, M ("GENERAL_DISABLED"));
        }
    }
}

void Locallab::inversretChanged ()
{

    if (batchMode) {
        if (inversret->get_inconsistent()) {
            inversret->set_inconsistent (false);
            inversretConn.block (true);
            inversret->set_active (false);
            inversretConn.block (false);
        } else if (lastinversret) {
            inversret->set_inconsistent (true);
        }

        lastinversret = inversret->get_active ();
    }

    if (inversret->get_active ()) {
        sensih->hide();
    } else {
        sensih->show();
    }


    if (listener) {
        if (getEnabled()) {
            listener->panelChanged (Evlocallabinversret, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabinversret, M ("GENERAL_DISABLED"));
        }
    }
}


void Locallab::setDefaults (const ProcParams * defParams, const ParamsEdited * pedited)
{
    degree->setDefault (defParams->locallab.degree);
    locY->setDefault (defParams->locallab.locY);
    locX->setDefault (defParams->locallab.locX);
    locYT->setDefault (defParams->locallab.locYT);
    locXL->setDefault (defParams->locallab.locXL);
    centerX->setDefault (defParams->locallab.centerX);
    centerY->setDefault (defParams->locallab.centerY);
    circrad->setDefault (defParams->locallab.circrad);
    thres->setDefault (defParams->locallab.thres);
    proxi->setDefault (defParams->locallab.proxi);
    lightness->setDefault (defParams->locallab.lightness);
    contrast->setDefault (defParams->locallab.contrast);
    chroma->setDefault (defParams->locallab.chroma);
    noiselumf->setDefault (defParams->locallab.noiselumf);
    noiselumc->setDefault (defParams->locallab.noiselumc);
    noisechrof->setDefault (defParams->locallab.noisechrof);
    noisechroc->setDefault (defParams->locallab.noisechroc);
    sharradius->setDefault (defParams->locallab.sharradius);
    sharamount->setDefault (defParams->locallab.sharamount);
    shardamping->setDefault (defParams->locallab.shardamping);
    shariter->setDefault (defParams->locallab.shariter);
    sensisha->setDefault (defParams->locallab.sensisha);
    sensi->setDefault (defParams->locallab.sensi);
    sensih->setDefault (defParams->locallab.sensih);
    retrab->setDefault (defParams->locallab.retrab);
    sensicb->setDefault (defParams->locallab.sensicb);
    sensibn->setDefault (defParams->locallab.sensibn);
    sensitm->setDefault (defParams->locallab.sensitm);
    transit->setDefault (defParams->locallab.transit);
    radius->setDefault (defParams->locallab.radius);
    strength->setDefault (defParams->locallab.strength);
    stren->setDefault (defParams->locallab.stren);
    gamma->setDefault (defParams->locallab.gamma);
    estop->setDefault (defParams->locallab.estop);
    gamma->setDefault (defParams->locallab.gamma);
    scaltm->setDefault (defParams->locallab.scaltm);
    rewei->setDefault (defParams->locallab.rewei);
    neigh->setDefault (defParams->locallab.neigh);
    nbspot->setDefault (defParams->locallab.nbspot);
    anbspot->setDefault (defParams->locallab.anbspot);
    hueref->setDefault (defParams->locallab.hueref);
    chromaref->setDefault (defParams->locallab.chromaref);
    lumaref->setDefault (defParams->locallab.lumaref);

    vart->setDefault (defParams->locallab.vart);
    chrrt->setDefault (defParams->locallab.chrrt);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setDefault (defParams->locallab.mult[i]);
    }

    threshold->setDefault (defParams->locallab.threshold);


    if (pedited) {
        degree->setDefaultEditedState (pedited->locallab.degree ? Edited : UnEdited);
        locY->setDefaultEditedState (pedited->locallab.locY ? Edited : UnEdited);
        locX->setDefaultEditedState (pedited->locallab.locX ? Edited : UnEdited);
        locYT->setDefaultEditedState (pedited->locallab.locYT ? Edited : UnEdited);
        locXL->setDefaultEditedState (pedited->locallab.locXL ? Edited : UnEdited);
        centerX->setDefaultEditedState (pedited->locallab.centerX ? Edited : UnEdited);
        centerY->setDefaultEditedState (pedited->locallab.centerY ? Edited : UnEdited);
        circrad->setDefaultEditedState (pedited->locallab.circrad ? Edited : UnEdited);
        thres->setDefaultEditedState (pedited->locallab.thres ? Edited : UnEdited);
        proxi->setDefaultEditedState (pedited->locallab.proxi ? Edited : UnEdited);
        lightness->setDefaultEditedState (pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setDefaultEditedState (pedited->locallab.chroma ? Edited : UnEdited);
        noiselumf->setDefaultEditedState (pedited->locallab.noiselumf ? Edited : UnEdited);
        noiselumc->setDefaultEditedState (pedited->locallab.noiselumc ? Edited : UnEdited);
        noisechrof->setDefaultEditedState (pedited->locallab.noisechrof ? Edited : UnEdited);
        noisechroc->setDefaultEditedState (pedited->locallab.noisechroc ? Edited : UnEdited);
        sharradius->setDefaultEditedState (pedited->locallab.sharradius ? Edited : UnEdited);
        sharamount->setDefaultEditedState (pedited->locallab.sharamount ? Edited : UnEdited);
        shardamping->setDefaultEditedState (pedited->locallab.shardamping ? Edited : UnEdited);
        shariter->setDefaultEditedState (pedited->locallab.shariter ? Edited : UnEdited);
        sensisha->setDefaultEditedState (pedited->locallab.sensisha ? Edited : UnEdited);
        sensi->setDefaultEditedState (pedited->locallab.sensi ? Edited : UnEdited);
        sensih->setDefaultEditedState (pedited->locallab.sensih ? Edited : UnEdited);
        retrab->setDefaultEditedState (pedited->locallab.retrab ? Edited : UnEdited);
        sensicb->setDefaultEditedState (pedited->locallab.sensicb ? Edited : UnEdited);
        sensibn->setDefaultEditedState (pedited->locallab.sensibn ? Edited : UnEdited);
        sensitm->setDefaultEditedState (pedited->locallab.sensitm ? Edited : UnEdited);
        radius->setDefaultEditedState (pedited->locallab.radius ? Edited : UnEdited);
        strength->setDefaultEditedState (pedited->locallab.strength ? Edited : UnEdited);
        stren->setDefaultEditedState (pedited->locallab.stren ? Edited : UnEdited);
        gamma->setDefaultEditedState (pedited->locallab.gamma ? Edited : UnEdited);
        estop->setDefaultEditedState (pedited->locallab.estop ? Edited : UnEdited);
        scaltm->setDefaultEditedState (pedited->locallab.scaltm ? Edited : UnEdited);
        rewei->setDefaultEditedState (pedited->locallab.rewei ? Edited : UnEdited);
        transit->setDefaultEditedState (pedited->locallab.transit ? Edited : UnEdited);
        str->setDefaultEditedState (pedited->locallab.str ? Edited : UnEdited);
        neigh->setDefaultEditedState (pedited->locallab.neigh ? Edited : UnEdited);
        nbspot->setDefaultEditedState (pedited->locallab.nbspot ? Edited : UnEdited);
        anbspot->setDefaultEditedState (pedited->locallab.anbspot ? Edited : UnEdited);
        hueref->setDefaultEditedState (pedited->locallab.hueref ? Edited : UnEdited);
        chromaref->setDefaultEditedState (pedited->locallab.chromaref ? Edited : UnEdited);
        lumaref->setDefaultEditedState (pedited->locallab.lumaref ? Edited : UnEdited);
        vart->setDefaultEditedState (pedited->locallab.vart ? Edited : UnEdited);
        chrrt->setDefaultEditedState (pedited->locallab.chrrt ? Edited : UnEdited);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState (pedited->locallab.mult[i] ? Edited : UnEdited);
        }

        threshold->setDefaultEditedState (pedited->locallab.threshold ? Edited : UnEdited);

    } else {
        degree->setDefaultEditedState (Irrelevant);
        locY->setDefaultEditedState (Irrelevant);
        locX->setDefaultEditedState (Irrelevant);
        locYT->setDefaultEditedState (Irrelevant);
        locXL->setDefaultEditedState (Irrelevant);
        centerX->setDefaultEditedState (Irrelevant);
        centerY->setDefaultEditedState (Irrelevant);
        circrad->setDefaultEditedState (Irrelevant);
        thres->setDefaultEditedState (Irrelevant);
        proxi->setDefaultEditedState (Irrelevant);
        lightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
        chroma->setDefaultEditedState (Irrelevant);
        noiselumf->setDefaultEditedState (Irrelevant);
        noiselumc->setDefaultEditedState (Irrelevant);
        noisechrof->setDefaultEditedState (Irrelevant);
        noisechroc->setDefaultEditedState (Irrelevant);
        sharradius->setDefaultEditedState (Irrelevant);
        sharamount->setDefaultEditedState (Irrelevant);
        shardamping->setDefaultEditedState (Irrelevant);
        shariter->setDefaultEditedState (Irrelevant);
        sensisha->setDefaultEditedState (Irrelevant);
        sensi->setDefaultEditedState (Irrelevant);
        sensih->setDefaultEditedState (Irrelevant);
        retrab->setDefaultEditedState (Irrelevant);
        sensicb->setDefaultEditedState (Irrelevant);
        sensibn->setDefaultEditedState (Irrelevant);
        sensitm->setDefaultEditedState (Irrelevant);
        radius->setDefaultEditedState (Irrelevant);
        strength->setDefaultEditedState (Irrelevant);
        stren->setDefaultEditedState (Irrelevant);
        gamma->setDefaultEditedState (Irrelevant);
        estop->setDefaultEditedState (Irrelevant);
        scaltm->setDefaultEditedState (Irrelevant);
        rewei->setDefaultEditedState (Irrelevant);
        transit->setDefaultEditedState (Irrelevant);
        str->setDefaultEditedState (Irrelevant);
        neigh->setDefaultEditedState (Irrelevant);
        nbspot->setDefaultEditedState (Irrelevant);
        anbspot->setDefaultEditedState (Irrelevant);
        hueref->setDefaultEditedState (Irrelevant);
        chromaref->setDefaultEditedState (Irrelevant);
        lumaref->setDefaultEditedState (Irrelevant);
        vart->setDefaultEditedState (Irrelevant);
        chrrt->setDefaultEditedState (Irrelevant);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState (Irrelevant);
        }

        threshold->setDefaultEditedState (Irrelevant);


    }
}

void Locallab::adjusterChanged (Adjuster * a, double newval)
{

    updateGeometry (int (centerX->getValue()), int (centerY->getValue()), int (circrad->getValue()), (int)locY->getValue(), degree->getValue(), (int)locX->getValue(), (int)locYT->getValue(), (int)locXL->getValue());
    anbspot->hide();
    hueref->hide();
    chromaref->hide();
    lumaref->hide();

    if (listener && getEnabled()) {
        if (a == degree) {
            listener->panelChanged (EvlocallabDegree, degree->getTextValue());
        } else if (a == locY) {
            if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) { // 0 2
                listener->panelChanged (EvlocallablocY, locY->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocY, locY->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged (EvlocallablocY, locY->getTextValue());
                locYT->setValue (locY->getValue());
            }
        } else if (a == locX) {
            //listener->panelChanged (EvlocallablocX, locX->getTextValue());
            if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocallablocX, locX->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged (EvlocallablocX, locX->getTextValue());
                locXL->setValue (locX->getValue());
            }
        } else if (a == locYT) {
            if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocallablocYT, locYT->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocYT, locYT->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            else if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
                listener->panelChanged (EvlocallablocYT, locYT->getTextValue());
                locYT->setValue (locY->getValue());
            }
        } else if (a == locXL) {
            if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocallablocXL, locXL->getTextValue());
                listener->panelChanged (EvlocallablocXL, locXL->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocXL, locXL->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged (EvlocallablocXL, locXL->getTextValue());
                locXL->setValue (locX->getValue());
            }
        } else if (a == lightness) {
            listener->panelChanged (Evlocallablightness, lightness->getTextValue());
        } else if (a == contrast) {
            listener->panelChanged (Evlocallabcontrast, contrast->getTextValue());
        } else if (a == chroma) {
            listener->panelChanged (Evlocallabchroma, chroma->getTextValue());
        } else if (a == noiselumf) {
            listener->panelChanged (Evlocallabnoiselumf, noiselumf->getTextValue());
        } else if (a == noiselumc) {
            listener->panelChanged (Evlocallabnoiselumc, noiselumc->getTextValue());
        } else if (a == noisechrof) {
            listener->panelChanged (Evlocallabnoisechrof, noisechrof->getTextValue());
        } else if (a == noisechroc) {
            listener->panelChanged (Evlocallabnoisechroc, noisechroc->getTextValue());
        } else if (a == sharradius) {
            listener->panelChanged (Evlocallabsharradius, sharradius->getTextValue());
        } else if (a == sharamount) {
            listener->panelChanged (Evlocallabsharamount, sharamount->getTextValue());
        } else if (a == shardamping) {
            listener->panelChanged (Evlocallabshardamping, shardamping->getTextValue());
        } else if (a == shariter) {
            listener->panelChanged (Evlocallabshariter, shariter->getTextValue());
        } else if (a == sensisha) {
            listener->panelChanged (Evlocallabsensis, sensisha->getTextValue());
        } else if (a == sensi) {
            listener->panelChanged (Evlocallabsensi, sensi->getTextValue());
        } else if (a == sensih) {
            listener->panelChanged (Evlocallabsensih, sensih->getTextValue());
        } else if (a == retrab) {
            listener->panelChanged (Evlocallabretrab, "");//retrab->getTextValue());
        } else if (a == radius) {
            listener->panelChanged (Evlocallabradius, radius->getTextValue());
        } else if (a == strength) {
            listener->panelChanged (Evlocallabstrength, strength->getTextValue());
        } else if (a == stren) {
            listener->panelChanged (Evlocallabstren, stren->getTextValue());
        } else if (a == gamma) {
            listener->panelChanged (Evlocallabgamma, gamma->getTextValue());
        } else if (a == estop) {
            listener->panelChanged (Evlocallabestop, estop->getTextValue());
        } else if (a == scaltm) {
            listener->panelChanged (Evlocallabscaltm, scaltm->getTextValue());
        } else if (a == rewei) {
            listener->panelChanged (Evlocallabrewei, rewei->getTextValue());
        } else if (a == sensitm) {
            listener->panelChanged (Evlocallabsensitm, sensitm->getTextValue());
        } else if (a == transit) {
            listener->panelChanged (Evlocallabtransit, transit->getTextValue());
        } else if (a == str) {
            listener->panelChanged (Evlocallabstr, str->getTextValue());
        } else if (a == neigh) {
            listener->panelChanged (Evlocallabneigh, neigh->getTextValue());
        } else if (a == nbspot) {
            listener->panelChanged (Evlocallabnbspot, nbspot->getTextValue());
        } else if (a == anbspot) {
            listener->panelChanged (Evlocallabanbspot, "");//anbspot->getTextValue());
        } else if (a == hueref) {
            listener->panelChanged (Evlocallabhueref, "");//anbspot->getTextValue());
        } else if (a == chromaref) {
            listener->panelChanged (Evlocallabchromaref, "");//anbspot->getTextValue());
        } else if (a == lumaref) {
            listener->panelChanged (Evlocallablumaref, "");//anbspot->getTextValue());
        } else if (a == vart) {
            listener->panelChanged (Evlocallabvart, vart->getTextValue());
        } else if (a == chrrt) {
            listener->panelChanged (Evlocallabchrrt, chrrt->getTextValue());
        } else if (a == circrad) {
            listener->panelChanged (Evlocallabcircrad, circrad->getTextValue());
        } else if (a == thres) {
            listener->panelChanged (Evlocallabthres, thres->getTextValue());
        } else if (a == threshold) {
            listener->panelChanged (EvlocallabThresho, threshold->getTextValue());
        } else if (a == sensicb) {
            listener->panelChanged (Evlocallabsensicb, sensicb->getTextValue());
        } else if (a == sensibn) {
            listener->panelChanged (Evlocallabsensibn, sensibn->getTextValue());
        } else if (a == proxi) {
            listener->panelChanged (Evlocallabproxi, proxi->getTextValue());
        } else if (a == centerX || a == centerY) {
            listener->panelChanged (EvlocallabCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
        } else {
            listener->panelChanged (EvlocallabEqualizer,
                                    Glib::ustring::compose ("%1, %2, %3, %4, %5",
                                            Glib::ustring::format (std::fixed, std::setprecision (0), multiplier[0]->getValue()),
                                            Glib::ustring::format (std::fixed, std::setprecision (0), multiplier[1]->getValue()),
                                            Glib::ustring::format (std::fixed, std::setprecision (0), multiplier[2]->getValue()),
                                            Glib::ustring::format (std::fixed, std::setprecision (0), multiplier[3]->getValue()),
                                            Glib::ustring::format (std::fixed, std::setprecision (0), multiplier[4]->getValue()))
                                   );

        }
    }
}

void Locallab::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvlocallabEnabled, M ("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvlocallabEnabled, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvlocallabEnabled, M ("GENERAL_DISABLED"));
        }
    }
}

void Locallab::avoidChanged ()
{

    if (batchMode) {
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

    if (listener) {
        if (getEnabled()) {
            listener->panelChanged (Evlocallabavoid, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabavoid, M ("GENERAL_DISABLED"));
        }
    }
}

void Locallab::setAdjusterBehavior (bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd,  bool strengthadd)
{
    degree->setAddMode (degreeadd);
    locY->setAddMode (locYadd);
    locX->setAddMode (locXadd);
    locYT->setAddMode (locYTadd);
    locXL->setAddMode (locXLadd);
    centerX->setAddMode (centeradd);
    centerY->setAddMode (centeradd);
    lightness->setAddMode (lightnessadd);
    contrast->setAddMode (contrastadd);
    chroma->setAddMode (chromaadd);
    sensi->setAddMode (sensiadd);
    transit->setAddMode (transitadd);
    radius->setAddMode (radiusadd);
    strength->setAddMode (strengthadd);

}

void Locallab::trimValues (rtengine::procparams::ProcParams * pp)
{
    degree->trimValue (pp->locallab.degree);
    locY->trimValue (pp->locallab.locY);
    locX->trimValue (pp->locallab.locX);
    locYT->trimValue (pp->locallab.locYT);
    locXL->trimValue (pp->locallab.locXL);
    centerX->trimValue (pp->locallab.centerX);
    centerY->trimValue (pp->locallab.centerY);
    circrad->trimValue (pp->locallab.circrad);
    thres->trimValue (pp->locallab.thres);
    proxi->trimValue (pp->locallab.proxi);
    lightness->trimValue (pp->locallab.lightness);
    contrast->trimValue (pp->locallab.contrast);
    chroma->trimValue (pp->locallab.chroma);
    noiselumf->trimValue (pp->locallab.noiselumf);
    noiselumc->trimValue (pp->locallab.noiselumc);
    noisechrof->trimValue (pp->locallab.noisechrof);
    noisechroc->trimValue (pp->locallab.noisechroc);
    sharradius->trimValue (pp->locallab.sharradius);
    sharamount->trimValue (pp->locallab.sharamount);
    shardamping->trimValue (pp->locallab.shardamping);
    shariter->trimValue (pp->locallab.shariter);
    sensisha->trimValue (pp->locallab.sensisha);
    sensi->trimValue (pp->locallab.sensi);
    sensih->trimValue (pp->locallab.sensih);
    retrab->trimValue (pp->locallab.retrab);
    sensicb->trimValue (pp->locallab.sensicb);
    sensibn->trimValue (pp->locallab.sensibn);
    sensitm->trimValue (pp->locallab.sensitm);
    radius->trimValue (pp->locallab.radius);
    strength->trimValue (pp->locallab.strength);
    stren->trimValue (pp->locallab.stren);
    gamma->trimValue (pp->locallab.gamma);
    estop->trimValue (pp->locallab.estop);
    scaltm->trimValue (pp->locallab.scaltm);
    rewei->trimValue (pp->locallab.rewei);
    transit->trimValue (pp->locallab.transit);
    str->trimValue (pp->locallab.str);
    neigh->trimValue (pp->locallab.neigh);
    nbspot->trimValue (pp->locallab.nbspot);
    anbspot->trimValue (pp->locallab.anbspot);
    hueref->trimValue (pp->locallab.hueref);
    chromaref->trimValue (pp->locallab.chromaref);
    lumaref->trimValue (pp->locallab.lumaref);

    vart->trimValue (pp->locallab.vart);
    chrrt->trimValue (pp->locallab.chrrt);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->trimValue (pp->locallab.mult[i]);
    }

    threshold->trimValue (pp->locallab.threshold);

}

void Locallab::setBatchMode (bool batchMode)
{
    removeIfThere (this, edit, false);
    ToolPanel::setBatchMode (batchMode);
    degree->showEditedCB ();
    locY->showEditedCB ();
    locX->showEditedCB ();
    locYT->showEditedCB ();
    locXL->showEditedCB ();
    centerX->showEditedCB ();
    centerY->showEditedCB ();
    circrad->showEditedCB ();
    thres->showEditedCB ();
    proxi->showEditedCB ();
    lightness->showEditedCB ();
    contrast->showEditedCB ();
    chroma->showEditedCB ();
    noiselumf->showEditedCB ();
    noiselumc->showEditedCB ();
    noisechroc->showEditedCB ();
    noiselumf->showEditedCB ();
    sharradius->showEditedCB ();
    sharamount->showEditedCB ();
    shardamping->showEditedCB ();
    shariter->showEditedCB ();
    sensisha->showEditedCB ();
    sensi->showEditedCB ();
    sensih->showEditedCB ();
    retrab->showEditedCB ();
    sensicb->showEditedCB ();
    sensibn->showEditedCB ();
    sensitm->showEditedCB ();
    radius->showEditedCB ();
    strength->showEditedCB ();
    stren->showEditedCB ();
    gamma->showEditedCB ();
    estop->showEditedCB ();
    scaltm->showEditedCB ();
    rewei->showEditedCB ();
    transit->showEditedCB ();
    Smethod->append (M ("GENERAL_UNCHANGED"));
    str->showEditedCB ();
    neigh->showEditedCB ();
    nbspot->showEditedCB ();
    anbspot->showEditedCB ();
    hueref->showEditedCB ();
    chromaref->showEditedCB ();
    lumaref->showEditedCB ();
    vart->showEditedCB ();
    LocalcurveEditorgainT->setBatchMode (batchMode);
    LocalcurveEditorgainTrab->setBatchMode (batchMode);
    llCurveEditorG->setBatchMode (batchMode);
//    llCurveEditorG2->setBatchMode (batchMode);
    chrrt->showEditedCB ();

    for (int i = 0; i < 5; i++) {
        multiplier[i]->showEditedCB();
    }

    threshold->showEditedCB();

}

void Locallab::setEditProvider (EditDataProvider * provider)
{
    EditSubscriber::setEditProvider (provider);
    cTgainshape->setEditProvider (provider);
    cTgainshaperab->setEditProvider (provider);

}

void Locallab::editToggled ()
{
    if (edit->get_active()) {
        subscribe();
    } else {
        unsubscribe();
    }
}

void Locallab::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R, G, B;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve

        Color::hsv2rgb01 (float (valX), float (valY), 0.5f, R, G, B);
    } else if (callerId == 2) {  // cc - bottom bar

        float value = (1.f - 0.7f) * float (valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01 (float (valY), float (valX), value, R, G, B);
    } else if (callerId == 3) {  // lc - bottom bar

        float value = (1.f - 0.7f) * float (valX) + 0.7f;
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01 (float (valY), float (valX), value, R, G, B);
    } else if (callerId == 4) {  // LH - bottom bar
        Color::hsv2rgb01 (float (valX), 0.5f, float (valY), R, G, B);
    } else if (callerId == 5) {  // HH - bottom bar
        float h = float ((valY - 0.5) * 0.3 + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01 (h, 0.5f, 0.5f, R, G, B);
    }

    caller->ccRed = double (R);
    caller->ccGreen = double (G);
    caller->ccBlue = double (B);
}



CursorShape Locallab::getCursor (int objectID)
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

bool Locallab::mouseOver (int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (editProvider && editProvider->object != lastObject) {
        if (lastObject > -1) {
            if (lastObject == 2 || lastObject == 3) {
                EditSubscriber::visibleGeometry.at (2)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at (3)->state = Geometry::NORMAL;

            } else if (lastObject == 0 || lastObject == 1) {
                EditSubscriber::visibleGeometry.at (0)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at (1)->state = Geometry::NORMAL;

            }

            else {
                EditSubscriber::visibleGeometry.at (4)->state = Geometry::NORMAL;
//               EditSubscriber::visibleGeometry.at (lastObject)->state = Geometry::NORMAL;
            }
        }

        if (editProvider->object > -1) {
            if (editProvider->object == 2 || editProvider->object == 3) {
                EditSubscriber::visibleGeometry.at (2)->state = Geometry::PRELIGHT;
                EditSubscriber::visibleGeometry.at (3)->state = Geometry::PRELIGHT;

            } else if (editProvider->object == 0 || editProvider->object == 1) {
                EditSubscriber::visibleGeometry.at (0)->state = Geometry::PRELIGHT;
                EditSubscriber::visibleGeometry.at (1)->state = Geometry::PRELIGHT;

            }

            else {
                EditSubscriber::visibleGeometry.at (4)->state = Geometry::PRELIGHT;
                //              EditSubscriber::visibleGeometry.at (editProvider->object)->state = Geometry::PRELIGHT;
            }
        }

        lastObject = editProvider->object;
        return true;
    }

    return false;
}

bool Locallab::button1Pressed (int modifierKey)
{
    if (lastObject < 0) {
        return false;
    }

    EditDataProvider *provider = getEditProvider();

    if (! (modifierKey & GDK_CONTROL_MASK)) {
        // button press is valid (no modifier key)
        PolarCoord pCoord;
        //  EditDataProvider *provider = getEditProvider();
        int imW, imH;
        provider->getImageSize (imW, imH);
        double halfSizeW = imW / 2.;
        double halfSizeH = imH / 2.;
        draggedCenter.set (int (halfSizeW + halfSizeW * (centerX->getValue() / 1000.)), int (halfSizeH + halfSizeH * (centerY->getValue() / 1000.)));

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
                draggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

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
                draggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

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
                draggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

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
                draggedlocXOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);
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
                draggedlocXOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

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
                draggedlocXOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

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
                EditSubscriber::visibleGeometry.at (2)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at (3)->state = Geometry::NORMAL;
            }

            if (lastObject == 0 || lastObject == 1) {
                EditSubscriber::visibleGeometry.at (0)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at (1)->state = Geometry::NORMAL;

            } else {
                EditSubscriber::visibleGeometry.at (4)->state = Geometry::NORMAL;
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

bool Locallab::drag1 (int modifierKey)
{
    // compute the polar coordinate of the mouse position
    EditDataProvider *provider = getEditProvider();
    int imW, imH;
    provider->getImageSize (imW, imH);
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
            double currDraggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

            if (lastObject == 2) {
                currDraggedlocYOffset -= draggedlocYOffset;
            }

            //else if (lastObject==3)
            // Dragging the lower locY bar
            //  currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            currDraggedlocYOffset = currDraggedlocYOffset * 2000. / verti;

            if (int (currDraggedlocYOffset) != locYT->getIntValue()) {
                locYT->setValue ((int (currDraggedlocYOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue() );

                if (listener) {
                    listener->panelChanged (EvlocallablocY, locYT->getTextValue());
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
            double currDraggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

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

                locY->setValue ((int (currDraggedlocYOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();

                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocY, locY->getTextValue());
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
            double currDraggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

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
                locY->setValue ((int (currDraggedlocYOffset)));
                //Smethod->get_active_row_number()==2
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();

                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(),  locYT->getValue(), locXL->getValue());

                if (listener) {
                    if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
                        listener->panelChanged (EvlocallablocY, locY->getTextValue());
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
            double currDraggedStrOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

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
                locX->setValue ((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
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
            double currDraggedStrOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

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
                locXL->setValue ((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
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
            double currDraggedStrOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

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
                locX->setValue ((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
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
        currPos.clip (imW, imH);
        int newCenterX = int ((double (currPos.x) - halfSizeW) / halfSizeW * 1000.);
        int newCenterY = int ((double (currPos.y) - halfSizeH) / halfSizeH * 1000.);

        if (newCenterX != centerX->getIntValue() || newCenterY != centerY->getIntValue()) {
            centerX->setValue (newCenterX);
            centerY->setValue (newCenterY);
            updateGeometry (newCenterX, newCenterY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

            if (listener) {
                listener->panelChanged (EvlocallabCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
            }

            return true;
        }
    }

    return false;
}

void Locallab::switchOffEditMode ()
{
    if (edit->get_active()) {
        // switching off the toggle button
        bool wasBlocked = editConn.block (true);
        edit->set_active (false);

        if (!wasBlocked) {
            editConn.block (false);
        }
    }

    EditSubscriber::switchOffEditMode();  // disconnect
}

