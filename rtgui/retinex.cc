/*
 *  This file is part of RawTherapee.
 */
#include "retinex.h"
#include "mycurve.h"

using namespace rtengine;
using namespace rtengine::procparams;

Retinex::Retinex () : FoldableToolPanel(this, "retinex", M("TP_RETINEX_LABEL"), false, true)
{
    CurveListener::setMulti(true);
    std::vector<GradientMilestone> milestones;
    nextmin = 0.;
    nextmax = 0.;
    nextminiT = 0.;
    nextmaxiT = 0.;
    nextmeanT = 0.;
    nextsigma = 0.;
    nextminT = 0.;
    nextmaxT = 0.;

    Gtk::VBox * retinexVBox = Gtk::manage ( new Gtk::VBox());
    retinexVBox->set_border_width(4);
    retinexVBox->set_spacing(2);

    Gtk::VBox * settingsVBox = Gtk::manage ( new Gtk::VBox());
    settingsVBox->set_border_width(4);
    settingsVBox->set_spacing(2);

    dhbox = Gtk::manage (new Gtk::HBox ());
    labmdh = Gtk::manage (new Gtk::Label (M("TP_RETINEX_METHOD") + ":"));
    dhbox->pack_start (*labmdh, Gtk::PACK_SHRINK, 1);

    retinexMethod = Gtk::manage (new MyComboBoxText ());
    retinexMethod->append_text (M("TP_RETINEX_LOW"));
    retinexMethod->append_text (M("TP_RETINEX_UNIFORM"));
    retinexMethod->append_text (M("TP_RETINEX_HIGH"));
    retinexMethod->set_active(0);
    retinexMethodConn = retinexMethod->signal_changed().connect ( sigc::mem_fun(*this, &Retinex::retinexMethodChanged) );
    retinexMethod->set_tooltip_markup (M("TP_RETINEX_METHOD_TOOLTIP"));

    retinexcolorspace = Gtk::manage (new MyComboBoxText ());
    retinexcolorspace->append_text (M("TP_RETINEX_LABSPACE"));
    retinexcolorspace->append_text (M("TP_RETINEX_HSLSPACE_LOG"));
    retinexcolorspace->append_text (M("TP_RETINEX_HSLSPACE_LIN"));
    retinexcolorspace->set_active(0);
    retinexColorSpaceConn = retinexcolorspace->signal_changed().connect ( sigc::mem_fun(*this, &Retinex::retinexColorSpaceChanged) );
    

    dhbox->pack_start(*retinexMethod);
    dhbox->pack_start(*retinexcolorspace);
    retinexVBox->pack_start(*dhbox);
    std::vector<double> defaultCurve;

    curveEditorGD = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_CONTEDIT_LAB"));
    curveEditorGD->setCurveListener (this);
    rtengine::RetinexParams::getDefaultCDCurve(defaultCurve);
    cdshape = static_cast<DiagonalCurveEditor*>(curveEditorGD->addCurve(CT_Diagonal, M("TP_RETINEX_CURVEEDITOR_CD")));
    cdshape->setResetCurve(DiagonalCurveType(defaultCurve.at(0)), defaultCurve);
    cdshape->setTooltip(M("TP_RETINEX_CURVEEDITOR_CD_TOOLTIP"));
    std::vector<GradientMilestone> milestones22;

    milestones22.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones22.push_back( GradientMilestone(1., 1., 1., 1.) );
    cdshape->setBottomBarBgGradient(milestones22);
    cdshape->setLeftBarBgGradient(milestones22);

    curveEditorGD->curveListComplete();

    curveEditorGDH = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_CONTEDIT_HSL"));
    curveEditorGDH->setCurveListener (this);
    rtengine::RetinexParams::getDefaultCDHCurve(defaultCurve);
    cdshapeH = static_cast<DiagonalCurveEditor*>(curveEditorGDH->addCurve(CT_Diagonal, M("TP_RETINEX_CURVEEDITOR_CD")));
    cdshapeH->setResetCurve(DiagonalCurveType(defaultCurve.at(0)), defaultCurve);
    cdshapeH->setTooltip(M("TP_RETINEX_CURVEEDITOR_CD_TOOLTIP"));
    std::vector<GradientMilestone> milestones22H;

    milestones22H.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones22H.push_back( GradientMilestone(1., 1., 1., 1.) );
    cdshapeH->setBottomBarBgGradient(milestones22H);
    cdshapeH->setLeftBarBgGradient(milestones22H);

    curveEditorGDH->curveListComplete();


    transmissionCurveEditorG = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_TRANSMISSION"));
    transmissionCurveEditorG->setCurveListener (this);

    rtengine::RetinexParams::getDefaulttransmissionCurve(defaultCurve);
    transmissionShape = static_cast<FlatCurveEditor*>(transmissionCurveEditorG->addCurve(CT_Flat, "", NULL, false));
    transmissionShape->setIdentityValue(0.);
    transmissionShape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
    transmissionShape->setBottomBarBgGradient(milestones);
    transmissionCurveEditorG->set_tooltip_markup (M("TP_RETINEX_TRANSMISSION_TOOLTIP"));

    transmissionCurveEditorG->curveListComplete();

    gambox = Gtk::manage (new Gtk::HBox ());
    labgam = Gtk::manage (new Gtk::Label (M("TP_RETINEX_GAM") + ":"));
    gambox->pack_start (*labgam, Gtk::PACK_SHRINK, 1);
    
    gammaretinex = Gtk::manage (new MyComboBoxText ());
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_NONE"));
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_LOW"));
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_MID"));
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_HIGH"));
    gammaretinex->set_active(0);
    gammaretinexConn = gammaretinex->signal_changed().connect ( sigc::mem_fun(*this, &Retinex::gammaretinexChanged) );
    gammaretinex->set_tooltip_markup (M("TP_RETINEX_GAMMA_TOOLTIP"));


    str = Gtk::manage (new Adjuster (M("TP_RETINEX_STRENGTH"), 0, 100., 1., 20.));
    neigh = Gtk::manage (new Adjuster (M("TP_RETINEX_NEIGHBOR"), 6, 100., 1., 80.));

    expsettings = new MyExpander (false, M("TP_RETINEX_SETTINGS"));
    expsettings->signal_button_release_event().connect_notify( sigc::bind( sigc::mem_fun(this, &Retinex::foldAllButMe), expsettings) );
    
    retinexVBox->pack_start (*str);
    str->show ();

    retinexVBox->pack_start (*neigh);
    neigh->show ();

    mMLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
    mMLabels->set_tooltip_markup (M("TP_RETINEX_MLABEL_TOOLTIP"));

    transLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
    transLabels->set_tooltip_markup (M("TP_RETINEX_TLABEL_TOOLTIP"));

    scal   = Gtk::manage (new Adjuster (M("TP_RETINEX_SCALES"), 1, 8., 1., 3.));
    gain   = Gtk::manage (new Adjuster (M("TP_RETINEX_GAIN"), 20, 200, 1, 50));
    offs   = Gtk::manage (new Adjuster (M("TP_RETINEX_OFFSET"), -10000, 10000, 1, 0));
    vart   = Gtk::manage (new Adjuster (M("TP_RETINEX_VARIANCE"), 50, 500, 1, 125));
    limd   = Gtk::manage (new Adjuster (M("TP_RETINEX_THRESHOLD"), 2, 100, 1, 8));
    gain->set_tooltip_markup (M("TP_RETINEX_GAIN_TOOLTIP"));
    scal->set_tooltip_markup (M("TP_RETINEX_SCALES_TOOLTIP"));
    vart->set_tooltip_markup (M("TP_RETINEX_VARIANCE_TOOLTIP"));
    limd->set_tooltip_markup (M("TP_RETINEX_THRESHOLD_TOOLTIP"));

    medianmap = Gtk::manage (new Gtk::CheckButton (M("TP_RETINEX_MEDIAN")));
    medianmap->set_active (true);
    medianmapConn  = medianmap->signal_toggled().connect( sigc::mem_fun(*this, &Retinex::medianmapChanged) );

    settingsVBox->pack_start (*mMLabels);
    mMLabels->show ();

    settingsVBox->pack_start (*transLabels);
    transLabels->show ();

    settingsVBox->pack_start (*curveEditorGD, Gtk::PACK_SHRINK, 4);
    curveEditorGD->show();

    settingsVBox->pack_start (*curveEditorGDH, Gtk::PACK_SHRINK, 4);
    curveEditorGDH->show();

    gambox->pack_start(*gammaretinex);
    
    settingsVBox->pack_start(*gambox);
    gammaretinex->show();
    
    settingsVBox->pack_start (*scal);
    scal->show ();

    settingsVBox->pack_start (*gain);
    gain->show ();

    settingsVBox->pack_start (*offs);
    offs->show ();

    settingsVBox->pack_start (*vart);
    vart->show ();

    settingsVBox->pack_start (*limd);
    limd->show ();

    settingsVBox->pack_start( *transmissionCurveEditorG, Gtk::PACK_SHRINK, 2);
    transmissionCurveEditorG->show();

    settingsVBox->pack_start (*medianmap);
    medianmap->show ();

    expsettings->add(*settingsVBox);

    str->setAdjusterListener (this);

    if (str->delay < 200) {
        str->delay = 200;
    }

    scal->setAdjusterListener (this);

    if (scal->delay < 200) {
        scal->delay = 200;
    }

    neigh->setAdjusterListener (this);

    if (neigh->delay < 200) {
        neigh->delay = 200;
    }

    gain->setAdjusterListener (this);

    if (gain->delay < 200) {
        gain->delay = 200;
    }

    offs->setAdjusterListener (this);

    if (offs->delay < 200) {
        offs->delay = 200;
    }

    vart->setAdjusterListener (this);

    if (vart->delay < 200) {
        vart->delay = 200;
    }

    limd->setAdjusterListener (this);

    if (limd->delay < 200) {
        limd->delay = 200;
    }

    pack_start (*retinexVBox);
    pack_start (*expsettings);

    disableListener();
    retinexColorSpaceChanged();
    gammaretinexChanged();
    medianmapChanged();
    enableListener();

}

Retinex::~Retinex()
{
    delete curveEditorGD;
    delete curveEditorGDH;
    delete transmissionCurveEditorG;

}

void Retinex::foldAllButMe (GdkEventButton* event, MyExpander *expander)
{
    if (event->button == 3) {
        expsettings->set_expanded(expsettings == expander);
    }
}

void Retinex::writeOptions(std::vector<int> &tpOpen)
{
    tpOpen.push_back (expsettings->get_expanded ());
}

void Retinex::updateToolState(std::vector<int> &tpOpen)
{
    if(tpOpen.size() == 9) {
        expsettings->set_expanded(tpOpen.at(0));
    }
}





int minmaxChangedUI (void* data)
{
    GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected
    (static_cast<Retinex*>(data))->minmaxComputed_ ();
    return 0;
}

void Retinex::minmaxChanged (double cdma, double cdmin, double mini, double maxi, double Tmean, double Tsigma, double Tmin, double Tmax)
{
    nextmin = cdmin;
    nextmax = cdma;
    nextminiT = mini;
    nextmaxiT = maxi;
    nextmeanT = Tmean;
    nextsigma = Tsigma;
    nextminT = Tmin;
    nextmaxT = Tmax;
    g_idle_add (minmaxChangedUI, this);

}

bool Retinex::minmaxComputed_ ()
{

    disableListener ();
    enableListener ();
    updateLabel ();
    updateTrans ();
    return false;

}
void Retinex::updateLabel ()
{
    if (!batchMode) {
        float nX, nY;
        nX = nextmin;
        nY = nextmax;
        {
            mMLabels->set_text(
                Glib::ustring::compose(M("TP_RETINEX_MLABEL"),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), nX),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), nY))
            );
        }
    }
}

void Retinex::updateTrans ()
{
    if (!batchMode) {
        float nm, nM, nZ, nA, nB, nS;
        nm = nextminiT;
        nM = nextmaxiT;
        nZ = nextmeanT;
        nA = nextminT;
        nB = nextmaxT;
        nS = nextsigma;
        {
            transLabels->set_text(
                Glib::ustring::compose(M("TP_RETINEX_TLABEL"),
                                       Glib::ustring::format(std::fixed, std::setprecision(1), nm),
                                       Glib::ustring::format(std::fixed, std::setprecision(1), nM),
                                       Glib::ustring::format(std::fixed, std::setprecision(1), nZ),
                                       Glib::ustring::format(std::fixed, std::setprecision(1), nS),
                                       Glib::ustring::format(std::fixed, std::setprecision(1), nA),
                                       Glib::ustring::format(std::fixed, std::setprecision(1), nB))
            );
        }
    }
}



void Retinex::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    retinexMethodConn.block(true);
    retinexColorSpaceConn.block(true);
    gammaretinexConn.block(true);


    if (pedited) {
        scal->setEditedState (pedited->retinex.scal ? Edited : UnEdited);
        neigh->setEditedState (pedited->retinex.neigh ? Edited : UnEdited);
        gain->setEditedState (pedited->retinex.gain ? Edited : UnEdited);
        offs->setEditedState (pedited->retinex.offs ? Edited : UnEdited);
        vart->setEditedState (pedited->retinex.vart ? Edited : UnEdited);
        limd->setEditedState (pedited->retinex.limd ? Edited : UnEdited);
        set_inconsistent (multiImage && !pedited->retinex.enabled);
        medianmap->set_inconsistent (!pedited->retinex.medianmap);


        if (!pedited->retinex.retinexMethod) {
            retinexMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->retinex.retinexcolorspace) {
            retinexcolorspace->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->retinex.gammaretinex) {
            gammaretinex->set_active_text(M("GENERAL_UNCHANGED"));
        }
        
        cdshape->setUnChanged  (!pedited->retinex.cdcurve);
        cdshapeH->setUnChanged  (!pedited->retinex.cdHcurve);
        transmissionShape->setUnChanged (!pedited->retinex.transmissionCurve);

    }

    neigh->setValue    (pp->retinex.neigh);
    gain->setValue      (pp->retinex.gain);
    offs->setValue  (pp->retinex.offs);
    str->setValue    (pp->retinex.str);
    scal->setValue      (pp->retinex.scal);
    vart->setValue  (pp->retinex.vart);
    limd->setValue  (pp->retinex.limd);

    setEnabled (pp->retinex.enabled);

    medianmapConn.block (true);
    medianmap->set_active (pp->retinex.medianmap);
    medianmapConn.block (false);
    lastmedianmap = pp->retinex.medianmap;

    if (pp->retinex.retinexMethod == "low") {
        retinexMethod->set_active (0);
    } else if (pp->retinex.retinexMethod == "uni") {
        retinexMethod->set_active (1);
    } else if (pp->retinex.retinexMethod == "high") {
        retinexMethod->set_active (2);
    }


    if (pp->retinex.retinexcolorspace == "Lab") {
        retinexcolorspace->set_active (0);
    } else if (pp->retinex.retinexcolorspace == "HSLLOG") {
        retinexcolorspace->set_active (1);
    } else if (pp->retinex.retinexcolorspace == "HSLLIN") {
        retinexcolorspace->set_active (2);
    }

    if (pp->retinex.gammaretinex == "none") {
        gammaretinex->set_active (0);
    } else if (pp->retinex.gammaretinex == "low") {
        gammaretinex->set_active (1);
    } else if (pp->retinex.gammaretinex == "mid") {
        gammaretinex->set_active (2);
    } else if (pp->retinex.gammaretinex == "hig") {
        gammaretinex->set_active (3);
    }
    
    retinexMethodChanged ();
    retinexColorSpaceChanged();
    gammaretinexChanged();

    medianmapConn.block(true);
    medianmapChanged ();
    medianmapConn.block(false);

    cdshape->setCurve  (pp->retinex.cdcurve);
    cdshapeH->setCurve  (pp->retinex.cdHcurve);

    retinexMethodConn.block(false);
    retinexColorSpaceConn.block(false);
    gammaretinexConn.block(false);
    transmissionShape->setCurve (pp->retinex.transmissionCurve);


    enableListener ();
}



void Retinex::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->retinex.str    = str->getValue ();
    pp->retinex.scal      = (int)scal->getValue ();
    pp->retinex.neigh    = neigh->getValue ();
    pp->retinex.gain      = (int)gain->getValue ();
    pp->retinex.offs  = (int)offs->getValue ();
    pp->retinex.vart  = (int)vart->getValue ();
    pp->retinex.limd  = (int)limd->getValue ();
    pp->retinex.cdcurve = cdshape->getCurve ();
    pp->retinex.cdHcurve = cdshapeH->getCurve ();
    pp->retinex.transmissionCurve = transmissionShape->getCurve ();
    pp->retinex.enabled      = getEnabled();
    pp->retinex.medianmap                = medianmap->get_active();

    if (pedited) {
        pedited->retinex.retinexMethod    = retinexMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->retinex.retinexcolorspace    = retinexcolorspace->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->retinex.gammaretinex    = gammaretinex->get_active_text() != M("GENERAL_UNCHANGED");

        //%%%%%%%%%%%%%%%%%%%%%%
        pedited->retinex.str   = str->getEditedState ();
        pedited->retinex.scal     = scal->getEditedState ();
        pedited->retinex.neigh   = neigh->getEditedState ();
        pedited->retinex.gain     = gain->getEditedState ();
        pedited->retinex.offs = offs->getEditedState ();
        pedited->retinex.vart = vart->getEditedState ();
        pedited->retinex.limd = limd->getEditedState ();
        pedited->retinex.cdcurve   = !cdshape->isUnChanged ();
        pedited->retinex.cdHcurve   = !cdshapeH->isUnChanged ();
        pedited->retinex.transmissionCurve  = !transmissionShape->isUnChanged ();
        pedited->retinex.enabled       = !get_inconsistent();
        pedited->retinex.medianmap       = !medianmap->get_inconsistent();

    }

    if (retinexMethod->get_active_row_number() == 0) {
        pp->retinex.retinexMethod = "low";
    } else if (retinexMethod->get_active_row_number() == 1) {
        pp->retinex.retinexMethod = "uni";
    } else if (retinexMethod->get_active_row_number() == 2) {
        pp->retinex.retinexMethod = "high";
    }

    if (retinexcolorspace->get_active_row_number() == 0) {
        pp->retinex.retinexcolorspace = "Lab";
    } else if (retinexcolorspace->get_active_row_number() == 1) {
        pp->retinex.retinexcolorspace = "HSLLOG";
    } else if (retinexcolorspace->get_active_row_number() == 2) {
        pp->retinex.retinexcolorspace = "HSLLIN";
    }
    
    if (gammaretinex->get_active_row_number() == 0) {
        pp->retinex.gammaretinex = "none";
    } else if (gammaretinex->get_active_row_number() == 1) {
        pp->retinex.gammaretinex = "low";
    } else if (gammaretinex->get_active_row_number() == 2) {
        pp->retinex.gammaretinex = "mid";
    } else if (gammaretinex->get_active_row_number() == 3) {
        pp->retinex.gammaretinex = "hig";
    }
    
}

void Retinex::retinexMethodChanged()
{
    if (listener) {
        listener->panelChanged (EvretinexMethod, retinexMethod->get_active_text ());
    }
}

void Retinex::ColorSpaceUpdateUI ()
{
    if (!batchMode) {
        if(retinexcolorspace->get_active_row_number() == 0) {
            curveEditorGD->show();
            curveEditorGDH->hide();
        } else if(retinexcolorspace->get_active_row_number() == 1) {
            curveEditorGD->hide();
            curveEditorGDH->show();
        } else if(retinexcolorspace->get_active_row_number() == 2) {
            curveEditorGD->hide();
            curveEditorGDH->show();
        }
    }
}


void Retinex::retinexColorSpaceChanged()
{

    ColorSpaceUpdateUI();

    if (listener) {
        listener->panelChanged (EvretinexColorSpace, retinexcolorspace->get_active_text ());
    }
}

void Retinex::gammaretinexChanged()
{

    ColorSpaceUpdateUI();

    if (listener) {
        listener->panelChanged (Evretinexgamma, gammaretinex->get_active_text ());
    }
}

void Retinex::medianmapChanged ()
{
    if (batchMode) {
        if (medianmap->get_inconsistent()) {
            medianmap->set_inconsistent (false);
            medianmapConn.block (true);
            medianmap->set_active (false);
            medianmapConn.block (false);
        } else if (lastmedianmap) {
            medianmap->set_inconsistent (true);
        }

        lastmedianmap = medianmap->get_active ();
    }

    if (listener) {
        if (medianmap->get_active()) {
            if (getEnabled()) {
                listener->panelChanged (EvRetinexmedianmap, M("GENERAL_ENABLED"));
            }
        } else {
            if (getEnabled()) {
                listener->panelChanged (EvRetinexmedianmap, M("GENERAL_DISABLED"));
            }
        }

    }
}


void Retinex::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    neigh->setDefault (defParams->retinex.neigh);
    gain->setDefault (defParams->retinex.gain);
    offs->setDefault (defParams->retinex.offs);
    str->setDefault (defParams->retinex.str);
    scal->setDefault (defParams->retinex.scal);
    vart->setDefault (defParams->retinex.vart);
    limd->setDefault (defParams->retinex.limd);

    if (pedited) {
        neigh->setDefaultEditedState (pedited->retinex.neigh ? Edited : UnEdited);
        gain->setDefaultEditedState (pedited->retinex.gain ? Edited : UnEdited);
        offs->setDefaultEditedState (pedited->retinex.offs ? Edited : UnEdited);
        str->setDefaultEditedState (pedited->retinex.str ? Edited : UnEdited);
        scal->setDefaultEditedState (pedited->retinex.scal ? Edited : UnEdited);
        vart->setDefaultEditedState (pedited->retinex.vart ? Edited : UnEdited);
        limd->setDefaultEditedState (pedited->retinex.limd ? Edited : UnEdited);

    } else {
        neigh->setDefaultEditedState (Irrelevant);
        gain->setDefaultEditedState (Irrelevant);
        offs->setDefaultEditedState (Irrelevant);
        vart->setDefaultEditedState (Irrelevant);
        limd->setDefaultEditedState (Irrelevant);
        str->setDefaultEditedState (Irrelevant);
        scal->setDefaultEditedState (Irrelevant);
    }
}

void Retinex::setAdjusterBehavior (bool strAdd, bool neighAdd, bool scalAdd, bool limdAdd, bool gainAdd, bool offsAdd, bool vartAdd)
{

    str->setAddMode(strAdd);
    neigh->setAddMode(neighAdd);
    scal->setAddMode(scalAdd);
    limd->setAddMode(limdAdd);
    gain->setAddMode(gainAdd);
    offs->setAddMode(offsAdd);
    vart->setAddMode(vartAdd);
}


void Retinex::adjusterChanged (Adjuster* a, double newval)
{

    if (!listener || !getEnabled()) {
        return;
    }

    if (a == neigh) {
        listener->panelChanged (EvLneigh, neigh->getTextValue());
    } else if (a == str) {
        listener->panelChanged (EvLstr, str->getTextValue());
    } else if (a == scal) {
        listener->panelChanged (EvLscal, scal->getTextValue());
    } else if (a == gain) {
        listener->panelChanged (EvLgain, gain->getTextValue());
    } else if (a == offs) {
        listener->panelChanged (EvLoffs, offs->getTextValue());
    } else if (a == vart) {
        listener->panelChanged (EvLvart, vart->getTextValue());
    } else if (a == limd) {
        listener->panelChanged (EvLlimd, limd->getTextValue());
    }

}



void Retinex::autoOpenCurve  ()
{
    cdshape->openIfNonlinear();
    cdshapeH->openIfNonlinear();
    transmissionShape->openIfNonlinear();

}


void Retinex::curveChanged (CurveEditor* ce)
{
    if (listener && getEnabled()) {
        if (ce == cdshape) {
            listener->panelChanged (EvLCDCurve, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == cdshapeH) {
            listener->panelChanged (EvLCDHCurve, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == transmissionShape) {
            listener->panelChanged (EvRetinextransmission, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void Retinex::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvRetinexEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvRetinexEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvRetinexEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void Retinex::trimValues (rtengine::procparams::ProcParams* pp)
{
    str->trimValue(pp->retinex.str);
    scal->trimValue(pp->retinex.scal);
    neigh->trimValue(pp->retinex.neigh);
    gain->trimValue(pp->retinex.gain);
    offs->trimValue(pp->retinex.offs);
    vart->trimValue(pp->retinex.vart);
    limd->trimValue(pp->retinex.limd);

}
void Retinex::updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve,/* LUTu & histCLurve, LUTu & histLLCurve,*/ LUTu & histLCAM,  LUTu & histCCAM, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histLRETI)
{

    cdshape->updateBackgroundHistogram (histLRETI);
    cdshapeH->updateBackgroundHistogram (histLRETI);
}

void Retinex::setBatchMode (bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    neigh->showEditedCB ();
    gain->showEditedCB ();
    offs->showEditedCB ();
    str->showEditedCB ();
    scal->showEditedCB ();
    vart->showEditedCB ();
    limd->showEditedCB ();
    curveEditorGD->setBatchMode (batchMode);
    curveEditorGDH->setBatchMode (batchMode);
    transmissionCurveEditorG->setBatchMode (batchMode);


}
