/*
 *  This file is part of RawTherapee.
 */
#include "retinex.h"
#include "mycurve.h"
#include "rtimage.h"

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
    retinexMethod->append_text (M("TP_RETINEX_HIGHLIG"));
//    retinexMethod->append_text (M("TP_RETINEX_HIGHLIGPLUS"));
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

    Gtk::VBox *equalVBox = Gtk::manage (new Gtk::VBox());

    equalFrame = Gtk::manage (new Gtk::Frame(M("TP_RETINEX_EQUAL")));

    // Histogram equalizer Lab curve
    curveEditorGD = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_CONTEDIT_LAB"));
    curveEditorGD->setCurveListener (this);
    cdshape = static_cast<DiagonalCurveEditor*>(curveEditorGD->addCurve(CT_Diagonal, M("TP_RETINEX_CURVEEDITOR_CD")));
    cdshape->setTooltip(M("TP_RETINEX_CURVEEDITOR_CD_TOOLTIP"));
    std::vector<GradientMilestone> milestones22;

    milestones22.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones22.push_back( GradientMilestone(1., 1., 1., 1.) );
    cdshape->setBottomBarBgGradient(milestones22);
    cdshape->setLeftBarBgGradient(milestones22);

    curveEditorGD->curveListComplete();


    // Histogram equalizer HSL curve
    curveEditorGDH = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_CONTEDIT_HSL"));
    curveEditorGDH->setCurveListener (this);
    cdshapeH = static_cast<DiagonalCurveEditor*>(curveEditorGDH->addCurve(CT_Diagonal, M("TP_RETINEX_CURVEEDITOR_CD")));
    cdshapeH->setTooltip(M("TP_RETINEX_CURVEEDITOR_CD_TOOLTIP"));
    std::vector<GradientMilestone> milestones22H;

    milestones22H.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones22H.push_back( GradientMilestone(1., 1., 1., 1.) );
    cdshapeH->setBottomBarBgGradient(milestones22H);
    cdshapeH->setLeftBarBgGradient(milestones22H);

    curveEditorGDH->curveListComplete();

    curveEditormap = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_CONTEDIT_MAP"));
    curveEditormap->setCurveListener (this);
    mapshape = static_cast<DiagonalCurveEditor*>(curveEditormap->addCurve(CT_Diagonal, M("TP_RETINEX_CURVEEDITOR_MAP")));
    mapshape->setTooltip(M("TP_RETINEX_CURVEEDITOR_MAP_TOOLTIP"));
    std::vector<GradientMilestone> milestones222;
    milestones222.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones222.push_back( GradientMilestone(1., 1., 1., 1.) );
    mapshape->setBottomBarBgGradient(milestones222);
    mapshape->setLeftBarBgGradient(milestones222);

    curveEditormap->curveListComplete();

    // Transmission map curve
    transmissionCurveEditorG = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_TRANSMISSION"));
    transmissionCurveEditorG->setCurveListener (this);

    std::vector<double> defaultCurve;
    rtengine::RetinexParams::getDefaulttransmissionCurve(defaultCurve);
    transmissionShape = static_cast<FlatCurveEditor*>(transmissionCurveEditorG->addCurve(CT_Flat, "", NULL, false));
    transmissionShape->setIdentityValue(0.);
    transmissionShape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
    transmissionShape->setBottomBarBgGradient(milestones);
    transmissionCurveEditorG->set_tooltip_markup (M("TP_RETINEX_TRANSMISSION_TOOLTIP"));

    transmissionCurveEditorG->curveListComplete();

    gambox = Gtk::manage (new Gtk::HBox ());
    labgam = Gtk::manage (new Gtk::Label (M("TP_RETINEX_GAMMA") + ":"));
    gambox->pack_start (*labgam, Gtk::PACK_SHRINK, 1);

    gammaretinex = Gtk::manage (new MyComboBoxText ());
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_NONE"));
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_LOW"));
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_MID"));
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_HIGH"));
    gammaretinex->append_text (M("TP_RETINEX_GAMMA_FREE"));
    gammaretinex->set_active(0);
    gammaretinexConn = gammaretinex->signal_changed().connect ( sigc::mem_fun(*this, &Retinex::gammaretinexChanged) );
    gammaretinex->set_tooltip_markup (M("TP_RETINEX_GAMMA_TOOLTIP"));

    gam = Gtk::manage (new Adjuster (M("TP_RETINEX_FREEGAMMA"), 0.6, 3.0, 0.01, 1.30));
    slope = Gtk::manage (new Adjuster (M("TP_RETINEX_SLOPE"), 1., 20., 0.1, 3.));

    str = Gtk::manage (new Adjuster (M("TP_RETINEX_STRENGTH"), 0, 100., 1., 20.));
    neigh = Gtk::manage (new Adjuster (M("TP_RETINEX_NEIGHBOR"), 6, 100., 1., 80.));
    highl   = Gtk::manage (new Adjuster (M("TP_RETINEX_HIGHLIGHT"), 1, 20, 1, 4));
    highl->set_tooltip_markup (M("TP_RETINEX_HIGHLIGHT_TOOLTIP"));
    vart   = Gtk::manage (new Adjuster (M("TP_RETINEX_VARIANCE"), 50, 500, 1, 200));
    vart->set_tooltip_markup (M("TP_RETINEX_VARIANCE_TOOLTIP"));

    expsettings = new MyExpander (false, M("TP_RETINEX_SETTINGS"));
    expsettings->signal_button_release_event().connect_notify( sigc::bind( sigc::mem_fun(this, &Retinex::foldAllButMe), expsettings) );

    retinexVBox->pack_start (*str);
    str->show ();

    retinexVBox->pack_start (*neigh);
    neigh->show ();

    retinexVBox->pack_start (*vart);
    vart->show ();

    retinexVBox->pack_start (*highl);
    highl->show ();

    mMLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
    mMLabels->set_tooltip_markup (M("TP_RETINEX_MLABEL_TOOLTIP"));

    transLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
    transLabels->set_tooltip_markup (M("TP_RETINEX_TLABEL_TOOLTIP"));
    transLabels2 = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));

    scal   = Gtk::manage (new Adjuster (M("TP_RETINEX_SCALES"), -1, 6., 1., 3.));
    iter   = Gtk::manage (new Adjuster (M("TP_RETINEX_ITER"), 1, 5., 1., 1.));
    grad   = Gtk::manage (new Adjuster (M("TP_RETINEX_GRAD"), -2., 2., 1., 1.));
    grads   = Gtk::manage (new Adjuster (M("TP_RETINEX_GRADS"), -2., 2., 1., 1.));
    gain   = Gtk::manage (new Adjuster (M("TP_RETINEX_GAIN"), 20, 200, 1, 50));
    offs   = Gtk::manage (new Adjuster (M("TP_RETINEX_OFFSET"), -1000, 5000, 1, 0));
//    vart   = Gtk::manage (new Adjuster (M("TP_RETINEX_VARIANCE"), 50, 500, 1, 125));
    limd   = Gtk::manage (new Adjuster (M("TP_RETINEX_THRESHOLD"), 2, 100, 1, 8));
    baselog   = Gtk::manage (new Adjuster (M("TP_RETINEX_BASELOG"), 1., 10., 1., 3.));
    skal   = Gtk::manage (new Adjuster (M("TP_RETINEX_SKAL"), 1, 8, 1, 3));
    gain->set_tooltip_markup (M("TP_RETINEX_GAIN_TOOLTIP"));
    scal->set_tooltip_markup (M("TP_RETINEX_SCALES_TOOLTIP"));
    iter->set_tooltip_markup (M("TP_RETINEX_ITER_TOOLTIP"));
    grad->set_tooltip_markup (M("TP_RETINEX_GRAD_TOOLTIP"));
    grads->set_tooltip_markup (M("TP_RETINEX_GRADS_TOOLTIP"));
//    vart->set_tooltip_markup (M("TP_RETINEX_VARIANCE_TOOLTIP"));
    limd->set_tooltip_markup (M("TP_RETINEX_THRESHOLD_TOOLTIP"));
    baselog->set_tooltip_markup (M("TP_RETINEX_BASELOG_TOOLTIP"));

    // Gain Transmission map curve
    gaintransmissionCurve = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_GAINTRANSMISSION"));
    gaintransmissionCurve->setCurveListener (this);

//    std::vector<double> defaultCurve;
    rtengine::RetinexParams::getDefaultgaintransmissionCurve(defaultCurve);
    gaintransmissionShape = static_cast<FlatCurveEditor*>(gaintransmissionCurve->addCurve(CT_Flat, "", NULL, false));
    gaintransmissionShape->setIdentityValue(0.);
    gaintransmissionShape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
    gaintransmissionShape->setBottomBarBgGradient(milestones);
    gaintransmissionCurve->set_tooltip_markup (M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));

    gaintransmissionCurve->curveListComplete();


    Gtk::Frame *p1Frame;
    p1Frame = Gtk::manage (new Gtk::Frame(M("TP_RETINEX_LABEL_MASK")) );
    p1Frame->set_border_width(0);
    p1Frame->set_label_align(0.025, 0.5);

    Gtk::VBox *p1VBox;
    p1VBox = Gtk::manage ( new Gtk::VBox());
    p1VBox->set_border_width(4);
    p1VBox->set_spacing(2);

    mapbox = Gtk::manage (new Gtk::HBox ());
    labmap = Gtk::manage (new Gtk::Label (M("TP_RETINEX_MAP") + ":"));
    mapbox->pack_start (*labmap, Gtk::PACK_SHRINK, 1);

    mapMethod = Gtk::manage (new MyComboBoxText ());
    mapMethod->append_text (M("TP_RETINEX_MAP_NONE"));
//    mapMethod->append_text (M("TP_RETINEX_MAP_CURV"));
    mapMethod->append_text (M("TP_RETINEX_MAP_GAUS"));
    mapMethod->append_text (M("TP_RETINEX_MAP_MAPP"));
    mapMethod->append_text (M("TP_RETINEX_MAP_MAPT"));
    mapMethod->set_active(0);
    mapMethodConn = mapMethod->signal_changed().connect ( sigc::mem_fun(*this, &Retinex::mapMethodChanged) );
    mapMethod->set_tooltip_markup (M("TP_RETINEX_MAP_METHOD_TOOLTIP"));

    highlights   = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), 0, 100, 1, 0));
    h_tonalwidth = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_HLTONALW"), 10, 100, 1, 80));
    shadows      = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_SHADOWS"), 0, 100, 1, 0));
    s_tonalwidth = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_SHTONALW"), 10, 100, 1, 80));
    radius = Gtk::manage (new Adjuster (M("TP_SHADOWSHLIGHTS_RADIUS"), 5, 100, 1, 40));

    viewbox = Gtk::manage (new Gtk::HBox ());
    labview = Gtk::manage (new Gtk::Label (M("TP_RETINEX_VIEW") + ":"));
    viewbox->pack_start (*labview, Gtk::PACK_SHRINK, 1);

    viewMethod = Gtk::manage (new MyComboBoxText ());
    viewMethod->append_text (M("TP_RETINEX_VIEW_NONE"));
    viewMethod->append_text (M("TP_RETINEX_VIEW_UNSHARP"));
    viewMethod->append_text (M("TP_RETINEX_VIEW_MASK"));
    viewMethod->append_text (M("TP_RETINEX_VIEW_TRAN"));
    viewMethod->append_text (M("TP_RETINEX_VIEW_TRAN2"));
    viewMethod->set_active(0);
    viewMethodConn = viewMethod->signal_changed().connect ( sigc::mem_fun(*this, &Retinex::viewMethodChanged) );
    viewMethod->set_tooltip_markup (M("TP_RETINEX_VIEW_METHOD_TOOLTIP"));

    curveEditorGH = new CurveEditorGroup (options.lastRetinexDir, M("TP_RETINEX_CONTEDIT_LH"));
    curveEditorGH->setCurveListener (this);

    lhshape = static_cast<FlatCurveEditor*>(curveEditorGH->addCurve(CT_Flat, M("TP_RETINEX_CURVEEDITOR_LH")));
    lhshape->setTooltip(M("TP_RETINEX_CURVEEDITOR_LH_TOOLTIP"));
    lhshape->setCurveColorProvider(this, 4);

    milestones.clear();

    for (int i = 0; i < 7; i++) {
        float R, G, B;
        float x = float(i) * (1.0f / 6.0);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        milestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
    }

    lhshape->setBottomBarBgGradient(milestones);

    curveEditorGH->curveListComplete();

    medianmap = Gtk::manage (new Gtk::CheckButton (M("TP_RETINEX_MEDIAN")));
    medianmap->set_active (true);
    medianmapConn  = medianmap->signal_toggled().connect( sigc::mem_fun(*this, &Retinex::medianmapChanged) );

    settingsVBox->pack_start (*mMLabels);
    mMLabels->show ();

    settingsVBox->pack_start (*transLabels);
    transLabels->show ();

    settingsVBox->pack_start (*transLabels2);
    transLabels2->show ();

    equalVBox->pack_start (*curveEditorGD, Gtk::PACK_SHRINK, 4);
    curveEditorGD->show();

    equalVBox->pack_start (*curveEditorGDH, Gtk::PACK_SHRINK, 4);
    curveEditorGDH->show();

    equalVBox->pack_start (*curveEditorGH, Gtk::PACK_SHRINK, 4);
    curveEditorGH->show();

    gambox->pack_start(*gammaretinex);

    equalVBox->pack_start(*gambox);
    gammaretinex->show();

    equalVBox->pack_start (*gam);
    gam->show ();

    equalVBox->pack_start (*slope);
    slope->show ();

    equalFrame->add(*equalVBox);
    settingsVBox->pack_start (*equalFrame);


    Gtk::VBox *iterVBox = Gtk::manage (new Gtk::VBox());

    iterFrame = Gtk::manage (new Gtk::Frame(M("TP_RETINEX_ITERF")));

    iterVBox->pack_start (*iter);
    iter->show ();

    iterVBox->pack_start (*scal);
    scal->show ();

    iterVBox->pack_start (*grad);
    grad->show ();

    iterVBox->pack_start (*grads);
    grads->show ();

    iterFrame->add(*iterVBox);
    settingsVBox->pack_start (*iterFrame);

    Gtk::VBox *tranVBox = Gtk::manage (new Gtk::VBox());

    tranFrame = Gtk::manage (new Gtk::Frame(M("TP_RETINEX_TRANF")));

    tranVBox->pack_start( *transmissionCurveEditorG, Gtk::PACK_SHRINK, 2);
    transmissionCurveEditorG->show();

    tranVBox->pack_start (*skal);
    skal->show ();

    tranVBox->pack_start (*limd);
    limd->show ();

    tranVBox->pack_start (*medianmap);
    medianmap->show ();

    tranFrame->add(*tranVBox);
    settingsVBox->pack_start (*tranFrame);

    Gtk::VBox *gainBox = Gtk::manage (new Gtk::VBox());

    Gtk::HSeparator *separator = Gtk::manage (new  Gtk::HSeparator());
    settingsVBox->pack_start(*separator, Gtk::PACK_SHRINK, 2);
    gainFrame = Gtk::manage (new Gtk::Frame(M("TP_RETINEX_GAINOFFS")));

    gainBox->pack_start( *gaintransmissionCurve, Gtk::PACK_SHRINK, 2);
    gaintransmissionCurve->show();

    gainBox->pack_start (*offs);
    offs->show ();
    gainFrame->add(*gainBox);
    settingsVBox->pack_start (*gainFrame);

    viewbox->pack_start(*viewMethod);
//   settingsVBox->pack_start(*viewbox);
    retinexVBox->pack_start(*viewbox);
    //settingsVBox->pack_start (*viewMethod);

    mapbox->pack_start(*mapMethod);
    // settingsVBox->pack_start(*mapbox);
    p1VBox->pack_start(*mapbox);

    p1VBox->pack_start (*curveEditormap, Gtk::PACK_SHRINK, 4);
    curveEditormap->show();

    p1VBox->pack_start (*highlights);
    highlights->show();
    p1VBox->pack_start (*h_tonalwidth);
    h_tonalwidth->show();
    p1VBox->pack_start (*shadows);
    shadows->show();
    p1VBox->pack_start (*s_tonalwidth);
    s_tonalwidth->show();
    p1VBox->pack_start (*radius);
    radius->show();



//    settingsVBox->pack_start (*highl);
//    highl->show ();

//   settingsVBox->pack_start (*baselog);
//   baselog->show ();

    //  settingsVBox->pack_start (*Gtk::manage (new  Gtk::HSeparator()));


    expsettings->add(*settingsVBox);

    neutrHBox = Gtk::manage (new Gtk::HBox ());
    neutrHBox->set_border_width (2);

    neutral = Gtk::manage (new Gtk::Button (M("TP_RETINEX_NEUTRAL")));
    RTImage *resetImg = Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
    neutral->set_image(*resetImg);
    neutral->set_tooltip_text (M("TP_RETINEX_NEUTRAL_TIP"));
    neutralconn = neutral->signal_pressed().connect( sigc::mem_fun(*this, &Retinex::neutral_pressed) );
    neutral->show();
    neutrHBox->pack_start (*neutral);



    str->setAdjusterListener (this);

    if (str->delay < 200) {
        str->delay = 200;
    }

    scal->setAdjusterListener (this);

    if (scal->delay < 200) {
        scal->delay = 200;
    }

    iter->setAdjusterListener (this);

    if (iter->delay < 200) {
        iter->delay = 200;
    }

    grad->setAdjusterListener (this);

    if (grad->delay < 200) {
        grad->delay = 200;
    }

    grads->setAdjusterListener (this);

    if (grads->delay < 200) {
        grads->delay = 200;
    }

    gam->setAdjusterListener (this);

    if (gam->delay < 500) {
        gam->delay = 500;
    }

    slope->setAdjusterListener (this);

    if (slope->delay < 500) {
        slope->delay = 500;
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

    highl->setAdjusterListener (this);

    if (highl->delay < 200) {
        highl->delay = 200;
    }

    baselog->setAdjusterListener (this);

    if (baselog->delay < 200) {
        baselog->delay = 200;
    }


    radius->setAdjusterListener (this);

    if (radius->delay < 200) {
        radius->delay = 200;
    }

    highlights->setAdjusterListener (this);

    if (highlights->delay < 200) {
        highlights->delay = 200;
    }

    h_tonalwidth->setAdjusterListener (this);

    if (h_tonalwidth->delay < 200) {
        h_tonalwidth->delay = 200;
    }

    shadows->setAdjusterListener (this);

    if (shadows->delay < 200) {
        shadows->delay = 200;
    }

    s_tonalwidth->setAdjusterListener (this);

    if (s_tonalwidth->delay < 200) {
        s_tonalwidth->delay = 200;
    }

    skal->setAdjusterListener (this);

    if (skal->delay < 200) {
        skal->delay = 200;
    }

    pack_start (*retinexVBox);
    p1Frame->add(*p1VBox);
    pack_start (*p1Frame, Gtk::PACK_EXPAND_WIDGET, 4);

    pack_start (*expsettings);
    pack_start (*neutrHBox);

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
    delete gaintransmissionCurve;
    delete curveEditorGH;
    delete curveEditormap;

}
void Retinex::neutral_pressed ()
{
    neigh->resetValue(false);
    gain->resetValue(false);
    offs->resetValue(false);
    str->resetValue(false);
    scal->resetValue(false);
    iter->resetValue(false);
    grad->resetValue(false);
    grads->resetValue(false);
    vart->resetValue(false);
    limd->resetValue(false);
    highl->resetValue(false);
    baselog->resetValue(false);
    gam->resetValue(false);
    slope->resetValue(false);
    highlights->resetValue(false);
    h_tonalwidth->resetValue(false);
    shadows->resetValue(false);
    s_tonalwidth->resetValue(false);
    radius->resetValue(false);
    mapMethod->set_active(0);
    viewMethod->set_active(0);
    retinexMethod->set_active(2);
    retinexcolorspace->set_active(0);
    gammaretinex->set_active(0);
    transmissionShape->reset();
    gaintransmissionShape->reset();
    cdshape->reset();
    cdshapeH->reset();
    lhshape->reset();
    mapshape->reset();
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
    if(tpOpen.size() == 10) {
        expsettings->set_expanded(tpOpen.at(9));
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
                                       Glib::ustring::format(std::fixed, std::setprecision(1), nS))
            );
            transLabels2->set_text(
                Glib::ustring::compose(M("TP_RETINEX_TLABEL2"),
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
    mapMethodConn.block(true);
    viewMethodConn.block(true);


    if (pedited) {
        scal->setEditedState (pedited->retinex.scal ? Edited : UnEdited);
        iter->setEditedState (pedited->retinex.iter ? Edited : UnEdited);
        grad->setEditedState (pedited->retinex.grad ? Edited : UnEdited);
        grads->setEditedState (pedited->retinex.grads ? Edited : UnEdited);
        neigh->setEditedState (pedited->retinex.neigh ? Edited : UnEdited);
        gam->setEditedState (pedited->retinex.gam ? Edited : UnEdited);
        slope->setEditedState (pedited->retinex.slope ? Edited : UnEdited);
        gain->setEditedState (pedited->retinex.gain ? Edited : UnEdited);
        offs->setEditedState (pedited->retinex.offs ? Edited : UnEdited);
        vart->setEditedState (pedited->retinex.vart ? Edited : UnEdited);
        limd->setEditedState (pedited->retinex.limd ? Edited : UnEdited);
        highl->setEditedState (pedited->retinex.highl ? Edited : UnEdited);
        baselog->setEditedState (pedited->retinex.baselog ? Edited : UnEdited);
        skal->setEditedState (pedited->retinex.skal ? Edited : UnEdited);
        set_inconsistent (multiImage && !pedited->retinex.enabled);
        medianmap->set_inconsistent (!pedited->retinex.medianmap);
        radius->setEditedState       (pedited->retinex.radius ? Edited : UnEdited);
        highlights->setEditedState   (pedited->retinex.highlights ? Edited : UnEdited);
        h_tonalwidth->setEditedState (pedited->retinex.htonalwidth ? Edited : UnEdited);
        shadows->setEditedState      (pedited->retinex.shadows ? Edited : UnEdited);
        s_tonalwidth->setEditedState (pedited->retinex.stonalwidth ? Edited : UnEdited);


        if (!pedited->retinex.retinexMethod) {
            retinexMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->retinex.mapMethod) {
            mapMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        if (!pedited->retinex.viewMethod) {
            viewMethod->set_active_text(M("GENERAL_UNCHANGED"));
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
        gaintransmissionShape->setUnChanged (!pedited->retinex.gaintransmissionCurve);
        lhshape->setUnChanged  (!pedited->retinex.lhcurve);
        mapshape->setUnChanged  (!pedited->retinex.mapcurve);

    }

    neigh->setValue    (pp->retinex.neigh);
    gain->setValue      (pp->retinex.gain);
    offs->setValue  (pp->retinex.offs);
    str->setValue    (pp->retinex.str);
    scal->setValue      (pp->retinex.scal);
    iter->setValue      (pp->retinex.iter);
    grad->setValue      (pp->retinex.grad);
    grads->setValue      (pp->retinex.grads);
    vart->setValue  (pp->retinex.vart);
    limd->setValue  (pp->retinex.limd);
    gam->setValue      (pp->retinex.gam);
    slope->setValue      (pp->retinex.slope);
    highl->setValue  (pp->retinex.highl);
    baselog->setValue  (pp->retinex.baselog);

    radius->setValue        (pp->retinex.radius);
    highlights->setValue    (pp->retinex.highlights);
    h_tonalwidth->setValue  (pp->retinex.htonalwidth);
    shadows->setValue       (pp->retinex.shadows);
    s_tonalwidth->setValue  (pp->retinex.stonalwidth);

    skal->setValue  (pp->retinex.skal);

    if(pp->retinex.iter == 1)   {
        grad->set_sensitive(false);
        scal->set_sensitive(false);
        grads->set_sensitive(false);
    } else {
        grad->set_sensitive(true);
        scal->set_sensitive(true);
        grads->set_sensitive(true);
    }

    setEnabled (pp->retinex.enabled);

    medianmap->set_active (pp->retinex.medianmap);
    medianmapConn.block (false);
    lastmedianmap = pp->retinex.medianmap;

    if (pp->retinex.retinexMethod == "low") {
        retinexMethod->set_active (0);
    } else if (pp->retinex.retinexMethod == "uni") {
        retinexMethod->set_active (1);
    } else if (pp->retinex.retinexMethod == "high") {
        retinexMethod->set_active (2);
    } else if (pp->retinex.retinexMethod == "highli") {
        retinexMethod->set_active (3);
//    } else if (pp->retinex.retinexMethod == "highliplus") {
//        retinexMethod->set_active (4);
    }

    if (pp->retinex.mapMethod == "none") {
        mapMethod->set_active (0);
//    } else if (pp->retinex.mapMethod == "curv") {
//        mapMethod->set_active (1);
    } else if (pp->retinex.mapMethod == "gaus") {
        mapMethod->set_active (1);
    } else if (pp->retinex.mapMethod == "map") {
        mapMethod->set_active (2);
    } else if (pp->retinex.mapMethod == "mapT") {
        mapMethod->set_active (3);
    }

    if (pp->retinex.viewMethod == "none") {
        viewMethod->set_active (0);
    } else if (pp->retinex.viewMethod == "unsharp") {
        viewMethod->set_active (1);
    } else if (pp->retinex.viewMethod == "mask") {
        viewMethod->set_active (2);
    } else if (pp->retinex.viewMethod == "tran") {
        viewMethod->set_active (3);
    } else if (pp->retinex.viewMethod == "tran2") {
        viewMethod->set_active (4);
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
    } else if (pp->retinex.gammaretinex == "fre") {
        gammaretinex->set_active (4);
    }

    retinexMethodChanged ();
    retinexColorSpaceChanged();
    gammaretinexChanged();
    mapMethodChanged ();
    viewMethodChanged ();

    medianmapConn.block(true);
    medianmapChanged ();
    medianmapConn.block(false);

    cdshape->setCurve  (pp->retinex.cdcurve);
    cdshapeH->setCurve  (pp->retinex.cdHcurve);
    lhshape->setCurve  (pp->retinex.lhcurve);
    mapshape->setCurve  (pp->retinex.mapcurve);

    retinexMethodConn.block(false);
    retinexColorSpaceConn.block(false);
    gammaretinexConn.block(false);
    mapMethodConn.block(false);
    viewMethodConn.block(false);
    transmissionShape->setCurve (pp->retinex.transmissionCurve);
    gaintransmissionShape->setCurve (pp->retinex.gaintransmissionCurve);


    enableListener ();
}



void Retinex::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->retinex.str    = str->getValue ();
    pp->retinex.scal      = (int)scal->getValue ();
    pp->retinex.iter      = (int) iter->getValue ();
    pp->retinex.grad      = (int) grad->getValue ();
    pp->retinex.grads      = (int) grads->getValue ();
    pp->retinex.gam      = gam->getValue ();
    pp->retinex.slope      = slope->getValue ();
    pp->retinex.neigh    = neigh->getValue ();
    pp->retinex.gain      = (int)gain->getValue ();
    pp->retinex.offs  = (int)offs->getValue ();
    pp->retinex.vart  = (int)vart->getValue ();
    pp->retinex.limd  = (int)limd->getValue ();
    pp->retinex.highl  = (int)highl->getValue ();
    pp->retinex.baselog  = baselog->getValue ();
    pp->retinex.skal  = (int)skal->getValue ();
    pp->retinex.cdcurve = cdshape->getCurve ();
    pp->retinex.lhcurve = lhshape->getCurve ();
    pp->retinex.cdHcurve = cdshapeH->getCurve ();
    pp->retinex.mapcurve = mapshape->getCurve ();
    pp->retinex.transmissionCurve = transmissionShape->getCurve ();
    pp->retinex.gaintransmissionCurve = gaintransmissionShape->getCurve ();
    pp->retinex.enabled      = getEnabled();
    pp->retinex.medianmap                = medianmap->get_active();

    pp->retinex.radius        = (int)radius->getValue ();
    pp->retinex.highlights    = (int)highlights->getValue ();
    pp->retinex.htonalwidth   = (int)h_tonalwidth->getValue ();
    pp->retinex.shadows       = (int)shadows->getValue ();
    pp->retinex.stonalwidth   = (int)s_tonalwidth->getValue ();

    if (pedited) {
        pedited->retinex.retinexMethod    = retinexMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->retinex.retinexcolorspace    = retinexcolorspace->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->retinex.gammaretinex    = gammaretinex->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->retinex.mapMethod    = mapMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->retinex.viewMethod    = viewMethod->get_active_text() != M("GENERAL_UNCHANGED");

        //%%%%%%%%%%%%%%%%%%%%%%
        pedited->retinex.str   = str->getEditedState ();
        pedited->retinex.scal     = scal->getEditedState ();
        pedited->retinex.iter     = iter->getEditedState ();
        pedited->retinex.grad     = grad->getEditedState ();
        pedited->retinex.grads     = grads->getEditedState ();
        pedited->retinex.gam     = gam->getEditedState ();
        pedited->retinex.slope     = slope->getEditedState ();
        pedited->retinex.neigh   = neigh->getEditedState ();
        pedited->retinex.gain     = gain->getEditedState ();
        pedited->retinex.offs = offs->getEditedState ();
        pedited->retinex.vart = vart->getEditedState ();
        pedited->retinex.limd = limd->getEditedState ();
        pedited->retinex.highl = highl->getEditedState ();
        pedited->retinex.baselog = baselog->getEditedState ();
        pedited->retinex.skal = skal->getEditedState ();
        pedited->retinex.cdcurve   = !cdshape->isUnChanged ();
        pedited->retinex.cdHcurve   = !cdshapeH->isUnChanged ();
        pedited->retinex.transmissionCurve  = !transmissionShape->isUnChanged ();
        pedited->retinex.gaintransmissionCurve  = !gaintransmissionShape->isUnChanged ();
        pedited->retinex.mapcurve   = !mapshape->isUnChanged ();
        pedited->retinex.enabled       = !get_inconsistent();
        pedited->retinex.medianmap       = !medianmap->get_inconsistent();
        pedited->retinex.lhcurve   = !lhshape->isUnChanged ();

        pedited->retinex.radius          = radius->getEditedState ();
        pedited->retinex.highlights      = highlights->getEditedState ();
        pedited->retinex.htonalwidth     = h_tonalwidth->getEditedState ();
        pedited->retinex.shadows         = shadows->getEditedState ();
        pedited->retinex.stonalwidth     = s_tonalwidth->getEditedState ();

    }

    if (retinexMethod->get_active_row_number() == 0) {
        pp->retinex.retinexMethod = "low";
    } else if (retinexMethod->get_active_row_number() == 1) {
        pp->retinex.retinexMethod = "uni";
    } else if (retinexMethod->get_active_row_number() == 2) {
        pp->retinex.retinexMethod = "high";
    } else if (retinexMethod->get_active_row_number() == 3) {
        pp->retinex.retinexMethod = "highli";
//    } else if (retinexMethod->get_active_row_number() == 4) {
//        pp->retinex.retinexMethod = "highliplus";
    }

    if (mapMethod->get_active_row_number() == 0) {
        pp->retinex.mapMethod = "none";
//   } else if (mapMethod->get_active_row_number() == 1) {
//       pp->retinex.mapMethod = "curv";
    } else if (mapMethod->get_active_row_number() == 1) {
        pp->retinex.mapMethod = "gaus";
    } else if (mapMethod->get_active_row_number() == 2) {
        pp->retinex.mapMethod = "map";
    } else if (mapMethod->get_active_row_number() == 3) {
        pp->retinex.mapMethod = "mapT";
    }

    if (viewMethod->get_active_row_number() == 0) {
        pp->retinex.viewMethod = "none";
    } else if (viewMethod->get_active_row_number() == 1) {
        pp->retinex.viewMethod = "unsharp";
    } else if (viewMethod->get_active_row_number() == 2) {
        pp->retinex.viewMethod = "mask";
    } else if (viewMethod->get_active_row_number() == 3) {
        pp->retinex.viewMethod = "tran";
    } else if (viewMethod->get_active_row_number() == 4) {
        pp->retinex.viewMethod = "tran2";
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
    } else if (gammaretinex->get_active_row_number() == 4) {
        pp->retinex.gammaretinex = "fre";
    }

}

void Retinex::retinexMethodChanged()
{

    if(retinexMethod->get_active_row_number() == 3) {
        highl->show();
    } else {
        highl->hide();
    }

    if (listener) {
        listener->panelChanged (EvretinexMethod, retinexMethod->get_active_text ());
    }
}



void Retinex::mapMethodChanged()
{

    if(mapMethod->get_active_row_number() == 1  /*|| mapMethod->get_active_row_number() == 2*/) {
        curveEditormap->show();
        highlights->show();
        h_tonalwidth->show();
        shadows->show();
        s_tonalwidth->show();
        radius->show();
    } else if(mapMethod->get_active_row_number() == 2  || mapMethod->get_active_row_number() == 3) {
        curveEditormap->show();
        highlights->show();
        h_tonalwidth->show();
        shadows->show();
        s_tonalwidth->show();
        radius->hide();
    } else {
        curveEditormap->hide();
        highlights->hide();
        h_tonalwidth->hide();
        shadows->hide();
        s_tonalwidth->hide();
        radius->hide();

    }

    if (listener) {
        listener->panelChanged (EvmapMethod, mapMethod->get_active_text ());
    }
}

void Retinex::viewMethodChanged()
{
    if(viewMethod->get_active_row_number() == 1 || viewMethod->get_active_row_number() == 2) {
        //    vart->hide();
        gain->hide();
        offs->hide();
        limd->hide();
        transmissionCurveEditorG->hide();
        medianmap->hide();
        iter->hide();
        scal->hide();
        grad->hide();
        grads->hide();
        curveEditorGH->hide();
    } else if(viewMethod->get_active_row_number() == 3 || viewMethod->get_active_row_number() == 4) {
        gain->hide();
        offs->hide();
        transmissionCurveEditorG->show();

        //    vart->hide();
        curveEditorGH->hide();
    } else {
        vart->show();
        neigh->show();
        gain->show();
        offs->show();
        limd->show();
        transmissionCurveEditorG->show();
        medianmap->show();
        iter->show();
        scal->show();
        grad->show();
        grads->show();
        curveEditorGH->show();
    }

    if (listener) {
        listener->panelChanged (EvviewMethod, viewMethod->get_active_text ());
    }
}



void Retinex::ColorSpaceUpdateUI ()
{
    if (!batchMode) {
        curveEditorGH->show();

        if(retinexcolorspace->get_active_row_number() == 0) {
            curveEditorGD->show();
            curveEditorGDH->hide();
            baselog->show();
        } else if(retinexcolorspace->get_active_row_number() == 1) {
            curveEditorGD->hide();
            curveEditorGDH->show();
            baselog->show();
        } else if(retinexcolorspace->get_active_row_number() == 2) {
            curveEditorGD->hide();
            curveEditorGDH->show();
            baselog->hide();
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
    if (!batchMode) {
        if(gammaretinex->get_active_row_number() == 4) {
            gam->show();
            slope->show();
        } else if(gammaretinex->get_active_row_number() != 4) {
            gam->hide();
            slope->hide();
        }
    }

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
    iter->setDefault (defParams->retinex.iter);
    grad->setDefault (defParams->retinex.grad);
    grads->setDefault (defParams->retinex.grads);
    vart->setDefault (defParams->retinex.vart);
    limd->setDefault (defParams->retinex.limd);
    highl->setDefault (defParams->retinex.highl);
    baselog->setDefault (defParams->retinex.baselog);
    skal->setDefault (defParams->retinex.skal);
    gam->setDefault (defParams->retinex.gam);
    slope->setDefault (defParams->retinex.slope);

    radius->setDefault (defParams->retinex.radius);
    highlights->setDefault (defParams->retinex.highlights);
    h_tonalwidth->setDefault (defParams->retinex.htonalwidth);
    shadows->setDefault (defParams->retinex.shadows);
    s_tonalwidth->setDefault (defParams->retinex.stonalwidth);

    if (pedited) {
        neigh->setDefaultEditedState (pedited->retinex.neigh ? Edited : UnEdited);
        gain->setDefaultEditedState (pedited->retinex.gain ? Edited : UnEdited);
        offs->setDefaultEditedState (pedited->retinex.offs ? Edited : UnEdited);
        str->setDefaultEditedState (pedited->retinex.str ? Edited : UnEdited);
        scal->setDefaultEditedState (pedited->retinex.scal ? Edited : UnEdited);
        iter->setDefaultEditedState (pedited->retinex.iter ? Edited : UnEdited);
        grad->setDefaultEditedState (pedited->retinex.grad ? Edited : UnEdited);
        grads->setDefaultEditedState (pedited->retinex.grads ? Edited : UnEdited);
        vart->setDefaultEditedState (pedited->retinex.vart ? Edited : UnEdited);
        limd->setDefaultEditedState (pedited->retinex.limd ? Edited : UnEdited);
        highl->setDefaultEditedState (pedited->retinex.highl ? Edited : UnEdited);
        baselog->setDefaultEditedState (pedited->retinex.baselog ? Edited : UnEdited);
        skal->setDefaultEditedState (pedited->retinex.skal ? Edited : UnEdited);
        gam->setDefaultEditedState (pedited->retinex.gam ? Edited : UnEdited);
        slope->setDefaultEditedState (pedited->retinex.slope ? Edited : UnEdited);

        radius->setDefaultEditedState       (pedited->retinex.radius ? Edited : UnEdited);
        highlights->setDefaultEditedState   (pedited->retinex.highlights ? Edited : UnEdited);
        h_tonalwidth->setDefaultEditedState (pedited->retinex.htonalwidth ? Edited : UnEdited);
        shadows->setDefaultEditedState      (pedited->retinex.shadows ? Edited : UnEdited);
        s_tonalwidth->setDefaultEditedState (pedited->retinex.stonalwidth ? Edited : UnEdited);

    } else {
        neigh->setDefaultEditedState (Irrelevant);
        gain->setDefaultEditedState (Irrelevant);
        offs->setDefaultEditedState (Irrelevant);
        vart->setDefaultEditedState (Irrelevant);
        limd->setDefaultEditedState (Irrelevant);
        highl->setDefaultEditedState (Irrelevant);
        baselog->setDefaultEditedState (Irrelevant);
        skal->setDefaultEditedState (Irrelevant);
        str->setDefaultEditedState (Irrelevant);
        scal->setDefaultEditedState (Irrelevant);
        iter->setDefaultEditedState (Irrelevant);
        grad->setDefaultEditedState (Irrelevant);
        grads->setDefaultEditedState (Irrelevant);
        gam->setDefaultEditedState (Irrelevant);
        slope->setDefaultEditedState (Irrelevant);

        radius->setDefaultEditedState       (Irrelevant);
        highlights->setDefaultEditedState   (Irrelevant);
        h_tonalwidth->setDefaultEditedState (Irrelevant);
        shadows->setDefaultEditedState      (Irrelevant);
        s_tonalwidth->setDefaultEditedState (Irrelevant);

    }
}

void Retinex::setAdjusterBehavior (bool strAdd, bool neighAdd, bool limdAdd, bool gainAdd, bool offsAdd, bool vartAdd, bool gamAdd, bool slopeAdd)
{
    str->setAddMode(strAdd);
    neigh->setAddMode(neighAdd);
    limd->setAddMode(limdAdd);
    gain->setAddMode(gainAdd);
    offs->setAddMode(offsAdd);
    vart->setAddMode(vartAdd);
    gam->setAddMode(gamAdd);
    slope->setAddMode(slopeAdd);
}


void Retinex::adjusterChanged (Adjuster* a, double newval)
{

    if (!listener || !getEnabled()) {
        return;
    }

    if(iter->getTextValue() > "1") {
        scal->set_sensitive(true);
        grad->set_sensitive(true);
        grads->set_sensitive(true);
    } else {
        scal->set_sensitive(false);
        grad->set_sensitive(false);
        grads->set_sensitive(false);
    }


    if (a == neigh) {
        listener->panelChanged (EvLneigh, neigh->getTextValue());
    } else if (a == str) {
        listener->panelChanged (EvLstr, str->getTextValue());
    } else if (a == scal) {
        listener->panelChanged (EvLscal, scal->getTextValue());
    } else if (a == iter) {
        listener->panelChanged (EvLiter, iter->getTextValue());
    } else if (a == grad) {
        listener->panelChanged (EvLgrad, grad->getTextValue());
    } else if (a == grads) {
        listener->panelChanged (EvLgrads, grads->getTextValue());
    } else if (a == gain) {
        listener->panelChanged (EvLgain, gain->getTextValue());
    } else if (a == offs) {
        listener->panelChanged (EvLoffs, offs->getTextValue());
    } else if (a == vart) {
        listener->panelChanged (EvLvart, vart->getTextValue());
    } else if (a == limd) {
        listener->panelChanged (EvLlimd, limd->getTextValue());
    } else if (a == highl) {
        listener->panelChanged (EvLhighl, highl->getTextValue());
    } else if (a == baselog) {
        listener->panelChanged (EvLbaselog, baselog->getTextValue());
    } else if (a == skal) {
        listener->panelChanged (EvLskal, skal->getTextValue());
    } else if (a == gam) {
        listener->panelChanged (EvLgam, gam->getTextValue());
    } else if (a == slope) {
        listener->panelChanged (EvLslope, slope->getTextValue());
    } else if (a == highlights) {
        listener->panelChanged (EvLhighlights, highlights->getTextValue());
    } else if (a == h_tonalwidth) {
        listener->panelChanged (EvLh_tonalwidth, h_tonalwidth->getTextValue());
    } else if (a == shadows) {
        listener->panelChanged (EvLshadows, shadows->getTextValue());
    } else if (a ==  s_tonalwidth) {
        listener->panelChanged (EvLs_tonalwidth,  s_tonalwidth->getTextValue());
    } else if (a ==  radius) {
        listener->panelChanged (EvLradius,  radius->getTextValue());

    }


}



void Retinex::autoOpenCurve  ()
{
    cdshape->openIfNonlinear();
    cdshapeH->openIfNonlinear();
    transmissionShape->openIfNonlinear();
    gaintransmissionShape->openIfNonlinear();
    lhshape->openIfNonlinear();
    mapshape->openIfNonlinear();

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
        } else if (ce == gaintransmissionShape) {
            listener->panelChanged (EvRetinexgaintransmission, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == lhshape) {
            listener->panelChanged (EvRetinexlhcurve, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == mapshape) {
            listener->panelChanged (EvRetinexmapcurve, M("HISTORY_CUSTOMCURVE"));
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
    iter->trimValue(pp->retinex.iter);
    grad->trimValue(pp->retinex.grad);
    grads->trimValue(pp->retinex.grads);
    neigh->trimValue(pp->retinex.neigh);
    gain->trimValue(pp->retinex.gain);
    offs->trimValue(pp->retinex.offs);
    vart->trimValue(pp->retinex.vart);
    limd->trimValue(pp->retinex.limd);
    highl->trimValue(pp->retinex.highl);
    baselog->trimValue(pp->retinex.baselog);
    gam->trimValue(pp->retinex.gam);
    slope->trimValue(pp->retinex.slope);
    highlights->trimValue(pp->retinex.highlights);
    shadows->trimValue(pp->retinex.shadows);


}
void Retinex::updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve,/* LUTu & histCLurve, LUTu & histLLCurve,*/ LUTu & histLCAM,  LUTu & histCCAM, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histLRETI)
{

    cdshape->updateBackgroundHistogram (histLRETI);
    cdshapeH->updateBackgroundHistogram (histLRETI);
}

void Retinex::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R, G, B;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve

        Color::hsv2rgb01(float(valX), float(valY), 0.5f, R, G, B);
    } else if (callerId == 2) {  // cc - bottom bar

        float value = (1.f - 0.7f) * float(valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
    } else if (callerId == 3) {  // lc - bottom bar

        float value = (1.f - 0.7f) * float(valX) + 0.7f;

        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
    } else if (callerId == 4) {  // LH - bottom bar
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



void Retinex::setBatchMode (bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    neigh->showEditedCB ();
    gain->showEditedCB ();
    offs->showEditedCB ();
    str->showEditedCB ();
    scal->showEditedCB ();
    iter->showEditedCB ();
    grad->showEditedCB ();
    grads->showEditedCB ();
    gam->showEditedCB ();
    slope->showEditedCB ();
    vart->showEditedCB ();
    limd->showEditedCB ();
    highl->showEditedCB ();
    baselog->showEditedCB ();

    radius->showEditedCB ();
    highlights->showEditedCB ();
    h_tonalwidth->showEditedCB ();
    shadows->showEditedCB ();
    s_tonalwidth->showEditedCB ();

    skal->showEditedCB ();
    curveEditorGD->setBatchMode (batchMode);
    curveEditorGDH->setBatchMode (batchMode);
    transmissionCurveEditorG->setBatchMode (batchMode);
    gaintransmissionCurve->setBatchMode (batchMode);
    curveEditorGH->setBatchMode (batchMode);
    curveEditormap->setBatchMode (batchMode);


}
