/*
 *  This file is part of RawTherapee.
 */
#include "colortoning.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "mycurve.h"
#include "rtimage.h"
#include "eventmapper.h"
#include "labgrid.h"
#include "options.h"
#include "rtengine/color.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring ColorToning::TOOL_NAME = "colortoning";

namespace {

constexpr int ID_LABREGION_HUE = 5;

inline bool hasMask(const std::vector<double> &dflt, const std::vector<double> &mask)
{
    return !(mask.empty() || mask[0] == FCT_Linear || mask == dflt);
}


inline float round_ab(float v)
{
    return int(v * 1000) / 1000.f;
}

} // namespace


ColorToning::ColorToning () : FoldableToolPanel(this, TOOL_NAME, M("TP_COLORTONING_LABEL"), false, true)
{
    nextbw = 0;
    CurveListener::setMulti(true);

    //---------------method

    method = Gtk::manage (new MyComboBoxText ());
    method->append (M("TP_COLORTONING_LAB"));
    method->append (M("TP_COLORTONING_RGBSLIDERS"));
    method->append (M("TP_COLORTONING_RGBCURVES"));
    method->append (M("TP_COLORTONING_SPLITCOCO"));
    method->append (M("TP_COLORTONING_SPLITLR"));
    method->append(M("TP_COLORTONING_LABGRID"));
    method->append(M("TP_COLORTONING_LABREGIONS"));
    method->set_active (0);
    method->set_tooltip_text (M("TP_COLORTONING_METHOD_TOOLTIP"));

    ctbox = Gtk::manage (new Gtk::Box ());
    Gtk::Label* lab = Gtk::manage (new Gtk::Label (M("TP_COLORTONING_METHOD")));
    ctbox->pack_start (*lab, Gtk::PACK_SHRINK, 4);
    ctbox->pack_start (*method);
    pack_start (*ctbox);

    methodconn = method->signal_changed().connect ( sigc::mem_fun(*this, &ColorToning::methodChanged) );

    //----------- Color curve ------------------------------

    colorSep = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    pack_start (*colorSep);

    colorCurveEditorG = new CurveEditorGroup (options.lastColorToningCurvesDir, M("TP_COLORTONING_COLOR"));
    colorCurveEditorG->setCurveListener (this);

    colorShape = static_cast<FlatCurveEditor*>(colorCurveEditorG->addCurve(CT_Flat, "", nullptr, false, false));
    colorShape->setCurveColorProvider(this, 4);
    std::vector<GradientMilestone> milestones;

    // whole hue range
    for (int i = 0; i < 7; i++) {
        float R, G, B;
        float x = float(i) * (1.0f / 6.f);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        milestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
    }

    colorShape->setLeftBarBgGradient(milestones);

    const ColorToningParams default_params;

    // luminance gradient
    milestones.clear();
    milestones.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones.push_back( GradientMilestone(1., 1., 1., 1.) );
    colorShape->setBottomBarBgGradient(milestones);
    colorShape->setResetCurve(FCT_MinMaxCPoints, default_params.colorCurve);

    // This will add the reset button at the end of the curveType buttons
    colorCurveEditorG->curveListComplete();
    colorCurveEditorG->show();

    pack_start( *colorCurveEditorG, Gtk::PACK_SHRINK, 2);

    //----------------------red green  blue yellow colours

    twocolor = Gtk::manage (new MyComboBoxText ());
    twocolor->append (M("TP_COLORTONING_TWOSTD"));
    twocolor->append (M("TP_COLORTONING_TWOALL"));
    twocolor->append (M("TP_COLORTONING_TWOBY"));
    twocolor->append (M("TP_COLORTONING_TWO2"));
    twocolor->set_tooltip_text (M("TP_COLORTONING_TWOCOLOR_TOOLTIP"));
    twocolor->set_active (0);

    twocconn = twocolor->signal_changed().connect( sigc::mem_fun(*this, &ColorToning::twoColorChangedByGui) );

    pack_start (*twocolor, Gtk::PACK_SHRINK, 4);

    //----------- Opacity curve ------------------------------

    opacityCurveEditorG = new CurveEditorGroup (options.lastColorToningCurvesDir, M("TP_COLORTONING_OPACITY"));
    opacityCurveEditorG->setCurveListener (this);

    opacityShape = static_cast<FlatCurveEditor*>(opacityCurveEditorG->addCurve(CT_Flat, "", nullptr, false, false));
    opacityShape->setIdentityValue(0.);
    opacityShape->setResetCurve(FlatCurveType(default_params.opacityCurve.at(0)), default_params.opacityCurve);
    opacityShape->setBottomBarBgGradient(milestones);

    // This will add the reset button at the end of the curveType buttons
    opacityCurveEditorG->curveListComplete();
    opacityCurveEditorG->show();

    pack_start( *opacityCurveEditorG, Gtk::PACK_SHRINK, 2);

    //---------Chroma curve 1 --------------------
    iby   = Gtk::manage (new RTImage ("circle-yellow-blue-small"));
    irg   = Gtk::manage (new RTImage ("circle-green-red-small"));

    clCurveEditorG = new CurveEditorGroup (options.lastColorToningCurvesDir, M("TP_COLORTONING_CHROMAC"));
    clCurveEditorG->setCurveListener (this);

    clshape = static_cast<DiagonalCurveEditor*>(clCurveEditorG->addCurve(CT_Diagonal, M("TP_COLORTONING_AB"), irg, false));
    clshape->setResetCurve(DiagonalCurveType(default_params.clcurve.at(0)), default_params.clcurve);
    clshape->setTooltip(M("TP_COLORTONING_CURVEEDITOR_CL_TOOLTIP"));

    clshape->setLeftBarColorProvider(this, 1);
    clshape->setRangeDefaultMilestones(0.25, 0.5, 0.75);
    milestones.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones.push_back( GradientMilestone(1., 1., 1., 1.) );

    clshape->setBottomBarBgGradient(milestones);
    clCurveEditorG->curveListComplete();

    pack_start( *clCurveEditorG, Gtk::PACK_SHRINK, 2);

    //---------Chroma curve 2 --------------------

    cl2CurveEditorG = new CurveEditorGroup (options.lastColorToningCurvesDir, M("TP_COLORTONING_CHROMAC"));
    cl2CurveEditorG->setCurveListener (this);

    cl2shape = static_cast<DiagonalCurveEditor*>(cl2CurveEditorG->addCurve(CT_Diagonal, M("TP_COLORTONING_BY"), iby, false));
    cl2shape->setResetCurve(DiagonalCurveType(default_params.cl2curve.at(0)), default_params.cl2curve);
    cl2shape->setTooltip(M("TP_COLORTONING_CURVEEDITOR_CL_TOOLTIP"));

    cl2shape->setLeftBarColorProvider(this, 1);
    cl2shape->setRangeDefaultMilestones(0.25, 0.5, 0.75);
    milestones.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones.push_back( GradientMilestone(1., 1., 1., 1.) );

    cl2shape->setBottomBarBgGradient(milestones);
    cl2CurveEditorG->curveListComplete();

    pack_start( *cl2CurveEditorG, Gtk::PACK_SHRINK, 2);

    //--------------------- Reset curves -----------------------------
    /*   Each curve can reset to a different curve, so this button only save one click now... so we remove it.
    neutralCurves = Gtk::manage (new Gtk::Button (M("TP_COLORTONING_NEUTRALCUR")));
    RTImage *resetImgc = Gtk::manage (new RTImage ("undo-small"));
    neutralCurves->set_image(*resetImgc);
    neutralCurves->set_tooltip_text (M("TP_COLORTONING_NEUTRALCUR_TIP"));
    neutralcurvesconn = neutralCurves->signal_pressed().connect( sigc::mem_fun(*this, &ColorToning::neutralCurves_pressed) );
    neutralCurves->show();

    pack_start (*neutralCurves);
    */

    //----------- Sliders + balance ------------------------------

    hlColSat = Gtk::manage (new ThresholdAdjuster (M("TP_COLORTONING_HIGHLIGHT"), 0., 100., 60., M("TP_COLORTONING_STRENGTH"), 1., 0., 360., 80., M("TP_COLORTONING_HUE"), 1., nullptr, false));
    hlColSat->setAdjusterListener (this);
    hlColSat->setBgColorProvider(this, 2);
    hlColSat->setUpdatePolicy(RTUP_DYNAMIC);

    pack_start( *hlColSat, Gtk::PACK_SHRINK, 0);

    shadowsColSat = Gtk::manage (new ThresholdAdjuster (M("TP_COLORTONING_SHADOWS"), 0., 100., 80., M("TP_COLORTONING_STRENGTH"), 1., 0., 360., 208., M("TP_COLORTONING_HUE"), 1., nullptr, false));
    shadowsColSat->setAdjusterListener (this);
    shadowsColSat->setBgColorProvider(this, 3);
    shadowsColSat->setUpdatePolicy(RTUP_DYNAMIC);

    pack_start( *shadowsColSat, Gtk::PACK_SHRINK, 0);


    balance = Gtk::manage( new Adjuster(M("TP_COLORTONING_BALANCE"), -100., 100., 1., 0.) );
    balance->setAdjusterListener(this);

    pack_start (*balance, Gtk::PACK_SHRINK, 2);

    //----------- Saturation and strength ------------------------------

//  satLimiterSep = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));


//  pack_start (*satLimiterSep, Gtk::PACK_SHRINK);

//  Gtk::Frame *p1Frame;
    // Vertical box container for the content of the Process 1 frame
    Gtk::Box* p1VBox;
    p1Frame = Gtk::manage (new Gtk::Frame(M("TP_COLORTONING_SA")) );
    p1Frame->set_label_align(0.025, 0.5);

    p1VBox = Gtk::manage ( new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    p1VBox->set_spacing(2);

    autosat = Gtk::manage (new Gtk::CheckButton (M("TP_COLORTONING_AUTOSAT")));
    autosat->set_active (true);
    autosatConn  = autosat->signal_toggled().connect( sigc::mem_fun(*this, &ColorToning::autosatChanged) );
    //satFrame->set_label_widget(*autosat);

    p1VBox->pack_start (*autosat, Gtk::PACK_SHRINK, 2);

    satProtectionThreshold = Gtk::manage( new Adjuster(M("TP_COLORTONING_SATURATIONTHRESHOLD"), 0., 100., 1., 80.) );
    satProtectionThreshold->setAdjusterListener(this);
    satProtectionThreshold->set_sensitive(false);

    p1VBox->pack_start( *satProtectionThreshold, Gtk::PACK_SHRINK, 2);

    saturatedOpacity = Gtk::manage( new Adjuster(M("TP_COLORTONING_SATURATEDOPACITY"), 0., 100., 1., 30.) );;
    saturatedOpacity->setAdjusterListener(this);
    saturatedOpacity->set_sensitive(false);

    p1VBox->pack_start( *saturatedOpacity, Gtk::PACK_SHRINK, 2); //I have moved after Chanmixer
    p1Frame->add(*p1VBox);
    pack_start (*p1Frame, Gtk::PACK_EXPAND_WIDGET, 4);

    strength = Gtk::manage( new Adjuster(M("TP_COLORTONING_STR"), 0., 100., 1., 50.) );;
    strength->setAdjusterListener(this);



    //  --------------------Sliders BW Colortoning -------------------

    chanMixerBox = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    Gtk::Box* chanMixerHLBox = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    Gtk::Box* chanMixerMidBox = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    Gtk::Box* chanMixerShadowsBox = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    Gtk::Image* iblueR   = Gtk::manage (new RTImage ("circle-blue-small"));
    Gtk::Image* iyelL    = Gtk::manage (new RTImage ("circle-yellow-small"));
    Gtk::Image* imagL    = Gtk::manage (new RTImage ("circle-magenta-small"));
    Gtk::Image* igreenR  = Gtk::manage (new RTImage ("circle-green-small"));
    Gtk::Image* icyanL   = Gtk::manage (new RTImage ("circle-blue-small"));
    Gtk::Image* iredR    = Gtk::manage (new RTImage ("circle-red-small"));

    Gtk::Image* iblueRm  = Gtk::manage (new RTImage ("circle-blue-small"));
    Gtk::Image* iyelLm   = Gtk::manage (new RTImage ("circle-yellow-small"));
    Gtk::Image* imagLm   = Gtk::manage (new RTImage ("circle-magenta-small"));
    Gtk::Image* igreenRm = Gtk::manage (new RTImage ("circle-green-small"));
    Gtk::Image* icyanLm  = Gtk::manage (new RTImage ("circle-blue-small"));
    Gtk::Image* iredRm   = Gtk::manage (new RTImage ("circle-red-small"));

    Gtk::Image* iblueRh  = Gtk::manage (new RTImage ("circle-blue-small"));
    Gtk::Image* iyelLh   = Gtk::manage (new RTImage ("circle-yellow-small"));
    Gtk::Image* imagLh   = Gtk::manage (new RTImage ("circle-magenta-small"));
    Gtk::Image* igreenRh = Gtk::manage (new RTImage ("circle-green-small"));
    Gtk::Image* icyanLh  = Gtk::manage (new RTImage ("circle-blue-small"));
    Gtk::Image* iredRh   = Gtk::manage (new RTImage ("circle-red-small"));

    redhigh   = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., icyanLh, iredRh  ));
    greenhigh = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., imagLh , igreenRh));
    bluehigh  = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., iyelLh , iblueRh ));

    redmed    = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., icyanLm, iredRm  ));
    greenmed  = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., imagLm , igreenRm));
    bluemed   = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., iyelLm , iblueRm ));

    redlow    = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., icyanL, iredR  ));
    greenlow  = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., imagL , igreenR));
    bluelow   = Gtk::manage (new Adjuster ("", -100., 100., 1., 0., iyelL , iblueR ));

    chanMixerHLBox->pack_start (*redhigh);
    chanMixerHLBox->pack_start (*greenhigh);
    chanMixerHLBox->pack_start (*bluehigh);
    chanMixerMidBox->pack_start (*redmed);
    chanMixerMidBox->pack_start (*greenmed);
    chanMixerMidBox->pack_start (*bluemed);
    chanMixerShadowsBox->pack_start (*redlow);
    chanMixerShadowsBox->pack_start (*greenlow);
    chanMixerShadowsBox->pack_start (*bluelow);

    Gtk::Frame *chanMixerHLFrame = Gtk::manage (new Gtk::Frame(M("TP_COLORTONING_HIGHLIGHT")));
    Gtk::Frame *chanMixerMidFrame = Gtk::manage (new Gtk::Frame(M("TP_COLORTONING_MIDTONES")));
    Gtk::Frame *chanMixerShadowsFrame = Gtk::manage (new Gtk::Frame(M("TP_COLORTONING_SHADOWS")));

    chanMixerHLFrame->set_label_align (0.025, 0.5);
    chanMixerMidFrame->set_label_align (0.025, 0.5);
    chanMixerShadowsFrame->set_label_align (0.025, 0.5);

    chanMixerHLFrame->add(*chanMixerHLBox);
    chanMixerMidFrame->add(*chanMixerMidBox);
    chanMixerShadowsFrame->add(*chanMixerShadowsBox);

    chanMixerBox->pack_start(*chanMixerHLFrame, Gtk::PACK_SHRINK);
    chanMixerBox->pack_start(*chanMixerMidFrame, Gtk::PACK_SHRINK);
    chanMixerBox->pack_start(*chanMixerShadowsFrame, Gtk::PACK_SHRINK);

    pack_start(*chanMixerBox, Gtk::PACK_SHRINK);
    pack_start( *strength, Gtk::PACK_SHRINK, 2); //I have moved after Chanmixer

    //--------------------- Reset sliders  ---------------------------
    neutrHBox = Gtk::manage (new Gtk::Box());

    neutral = Gtk::manage (new Gtk::Button (M("TP_COLORTONING_NEUTRAL")));
    neutral->set_tooltip_text (M("TP_COLORTONING_NEUTRAL_TOOLTIP"));
    neutralconn = neutral->signal_pressed().connect( sigc::mem_fun(*this, &ColorToning::neutral_pressed) );
    neutral->show();
    neutrHBox->pack_start (*neutral);

    pack_start (*neutrHBox);

    //--------------------- Keep luminance checkbox -------------------
    lumamode = Gtk::manage (new Gtk::CheckButton (M("TP_COLORTONING_LUMAMODE")));
    lumamode->set_tooltip_markup (M("TP_COLORTONING_LUMAMODE_TOOLTIP"));
    lumamode->set_active (false);
    lumamode->show ();
    lumamodeConn = lumamode->signal_toggled().connect( sigc::mem_fun(*this, &ColorToning::lumamodeChanged) );

    pack_start (*lumamode);


    redlow->setAdjusterListener (this);
    greenlow->setAdjusterListener (this);
    bluelow->setAdjusterListener (this);
    balance->setAdjusterListener (this);
    redmed->setAdjusterListener (this);
    greenmed->setAdjusterListener (this);
    bluemed->setAdjusterListener (this);
    redhigh->setAdjusterListener (this);
    greenhigh->setAdjusterListener (this);
    bluehigh->setAdjusterListener (this);

    //------------------------------------------------------------------------
    // LAB grid
    auto m = ProcEventMapper::getInstance();
//    EvColorToningLabGridValue = m->newEvent(RGBCURVE, "HISTORY_MSG_COLORTONING_LABGRID_VALUE");
    EvColorToningLabGridValue = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABGRID_VALUE");
    labgrid = Gtk::manage(new LabGrid(EvColorToningLabGridValue, M("TP_COLORTONING_LABGRID_VALUES")));
    pack_start(*labgrid, Gtk::PACK_EXPAND_WIDGET, 4);
    //------------------------------------------------------------------------

    //------------------------------------------------------------------------
    // LAB regions

    const auto add_button =
        [&](Gtk::Button *btn, Gtk::Box *box) -> void
        {
            setExpandAlignProperties(btn, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);
            btn->set_relief(Gtk::RELIEF_NONE);
            btn->get_style_context()->add_class(GTK_STYLE_CLASS_FLAT);
            btn->set_can_focus(false);
            btn->set_size_request(-1, 20);
            box->pack_start(*btn, false, false);
        };

    EvLabRegionList = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_LIST");
    EvLabRegionAB = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_AB");
    EvLabRegionSaturation = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_SATURATION");
    EvLabRegionLightness = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_LIGHTNESS");
    EvLabRegionSlope = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_SLOPE");
    EvLabRegionOffset = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_OFFSET");
    EvLabRegionPower = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_POWER");
    EvLabRegionHueMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_HUEMASK");
    EvLabRegionChromaticityMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_CHROMATICITYMASK");
    EvLabRegionLightnessMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_LIGHTNESSMASK");
    EvLabRegionMaskBlur = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_MASKBLUR");
    EvLabRegionShowMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_SHOWMASK");
    EvLabRegionChannel = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORTONING_LABREGION_CHANNEL");
    labRegionBox = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    labRegionList = Gtk::manage(new Gtk::ListViewText(3));
    labRegionList->set_size_request(-1, 150);
    labRegionList->set_can_focus(false);
    labRegionList->set_column_title(0, "#");
    labRegionList->set_column_title(1, M("TP_COLORTONING_LABREGION_LIST_TITLE"));
    labRegionList->set_column_title(2, M("TP_COLORTONING_LABREGION_MASK"));
    labRegionList->set_activate_on_single_click(true);
    labRegionSelectionConn = labRegionList->get_selection()->signal_changed().connect(sigc::mem_fun(this, &ColorToning::onLabRegionSelectionChanged));
    Gtk::Box *hb = Gtk::manage(new Gtk::Box());
    hb->pack_start(*labRegionList, Gtk::PACK_EXPAND_WIDGET);
    Gtk::Box* vb = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    labRegionAdd = Gtk::manage(new Gtk::Button());
    labRegionAdd->add(*Gtk::manage(new RTImage("add-small", Gtk::ICON_SIZE_BUTTON)));
    labRegionAdd->signal_clicked().connect(sigc::mem_fun(*this, &ColorToning::labRegionAddPressed));
    add_button(labRegionAdd, vb);
    labRegionRemove = Gtk::manage(new Gtk::Button());
    labRegionRemove->add(*Gtk::manage(new RTImage("remove-small", Gtk::ICON_SIZE_BUTTON)));
    labRegionRemove->signal_clicked().connect(sigc::mem_fun(*this, &ColorToning::labRegionRemovePressed));
    add_button(labRegionRemove, vb);
    labRegionUp = Gtk::manage(new Gtk::Button());
    labRegionUp->add(*Gtk::manage(new RTImage("arrow-up-small", Gtk::ICON_SIZE_BUTTON)));
    labRegionUp->signal_clicked().connect(sigc::mem_fun(*this, &ColorToning::labRegionUpPressed));
    add_button(labRegionUp, vb);
    labRegionDown = Gtk::manage(new Gtk::Button());
    labRegionDown->add(*Gtk::manage(new RTImage("arrow-down-small", Gtk::ICON_SIZE_BUTTON)));
    labRegionDown->signal_clicked().connect(sigc::mem_fun(*this, &ColorToning::labRegionDownPressed));
    add_button(labRegionDown, vb);
    labRegionCopy = Gtk::manage(new Gtk::Button());
    labRegionCopy->add(*Gtk::manage(new RTImage("arrow-right-small", Gtk::ICON_SIZE_BUTTON)));
    labRegionCopy->signal_clicked().connect(sigc::mem_fun(*this, &ColorToning::labRegionCopyPressed));
    add_button(labRegionCopy, vb);
    hb->pack_start(*vb, Gtk::PACK_SHRINK);
    labRegionBox->pack_start(*hb, true, true);

    labRegionAB = Gtk::manage(new LabGrid(EvLabRegionAB, M("TP_COLORTONING_LABREGION_ABVALUES"), false));
    labRegionBox->pack_start(*labRegionAB);

    labRegionSaturation = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_SATURATION"), -100, 100, 1, 0));
    labRegionSaturation->setAdjusterListener(this);
    labRegionBox->pack_start(*labRegionSaturation);

    labRegionSlope = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_SLOPE"), 0.1, 4.0, 0.001, 1));
    labRegionSlope->setLogScale(4, 0.1);
    labRegionSlope->setAdjusterListener(this);
    labRegionBox->pack_start(*labRegionSlope);
    labRegionOffset = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_OFFSET"), -0.1, 0.1, 0.001, 0));
    labRegionOffset->setAdjusterListener(this);
    labRegionBox->pack_start(*labRegionOffset);
    labRegionPower = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_POWER"), 0.1, 4.0, 0.001, 1));
    labRegionPower->setAdjusterListener(this);
    labRegionPower->setLogScale(4, 0.1);
    labRegionBox->pack_start(*labRegionPower);

    hb = Gtk::manage(new Gtk::Box());
    labRegionChannel = Gtk::manage(new MyComboBoxText());
    labRegionChannel->append(M("TP_COLORTONING_LABREGION_CHANNEL_ALL"));
    labRegionChannel->append(M("TP_COLORTONING_LABREGION_CHANNEL_R"));
    labRegionChannel->append(M("TP_COLORTONING_LABREGION_CHANNEL_G"));
    labRegionChannel->append(M("TP_COLORTONING_LABREGION_CHANNEL_B"));
    labRegionChannel->set_active(0);
    labRegionChannel->signal_changed().connect(sigc::mem_fun(*this, &ColorToning::labRegionChannelChanged));

    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_COLORTONING_LABREGION_CHANNEL") + ": ")), Gtk::PACK_SHRINK);
    hb->pack_start(*labRegionChannel);
    labRegionBox->pack_start(*hb);

    labRegionBox->pack_start(*Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL)));

    CurveEditorGroup *labRegionEditorG = Gtk::manage(new CurveEditorGroup(options.lastColorToningCurvesDir, M("TP_COLORTONING_LABREGION_MASK")));
    labRegionEditorG->setCurveListener(this);

    labRegionHueMask = static_cast<FlatCurveEditor *>(labRegionEditorG->addCurve(CT_Flat, M("TP_COLORTONING_LABREGION_HUEMASK"), nullptr, false, true));
    labRegionHueMask->setIdentityValue(0.);
    labRegionHueMask->setResetCurve(FlatCurveType(default_params.labregions[0].hueMask[0]), default_params.labregions[0].hueMask);
    labRegionHueMask->setCurveColorProvider(this, ID_LABREGION_HUE);
    labRegionHueMask->setBottomBarColorProvider(this, ID_LABREGION_HUE);
    labRegionHueMask->setEditID(EUID_Lab_HHCurve, BT_SINGLEPLANE_FLOAT);

    labRegionChromaticityMask = static_cast<FlatCurveEditor *>(labRegionEditorG->addCurve(CT_Flat, M("TP_COLORTONING_LABREGION_CHROMATICITYMASK"), nullptr, false, false));
    labRegionChromaticityMask->setIdentityValue(0.);
    labRegionChromaticityMask->setResetCurve(FlatCurveType(default_params.labregions[0].chromaticityMask[0]), default_params.labregions[0].chromaticityMask);
    labRegionChromaticityMask->setBottomBarColorProvider(this, ID_LABREGION_HUE+1);
    labRegionChromaticityMask->setEditID(EUID_Lab_CCurve, BT_SINGLEPLANE_FLOAT);

    labRegionLightnessMask = static_cast<FlatCurveEditor *>(labRegionEditorG->addCurve(CT_Flat, M("TP_COLORTONING_LABREGION_LIGHTNESSMASK"), nullptr, false, false));
    labRegionLightnessMask->setIdentityValue(0.);
    labRegionLightnessMask->setResetCurve(FlatCurveType(default_params.labregions[0].lightnessMask[0]), default_params.labregions[0].lightnessMask);
    labRegionLightnessMask->setBottomBarBgGradient(milestones);
    labRegionLightnessMask->setEditID(EUID_Lab_LCurve, BT_SINGLEPLANE_FLOAT);

    labRegionData = default_params.labregions;
    labRegionSelected = 0;
    {
        auto n = labRegionList->append("1");
        labRegionList->set_text(n, 1, "a=0 b=0 s=0 l=0");
    }

    labRegionEditorG->curveListComplete();
    labRegionEditorG->show();
    labRegionBox->pack_start(*labRegionEditorG, Gtk::PACK_SHRINK, 2);

    labRegionMaskBlur = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_MASKBLUR"), -10, 100, 0.1, 0));
    labRegionMaskBlur->setLogScale(10, 0);
    labRegionMaskBlur->setAdjusterListener(this);
    labRegionBox->pack_start(*labRegionMaskBlur);

    labRegionShowMask = Gtk::manage(new Gtk::CheckButton(M("TP_COLORTONING_LABREGION_SHOWMASK")));
    labRegionShowMask->signal_toggled().connect(sigc::mem_fun(*this, &ColorToning::labRegionShowMaskChanged));
    labRegionBox->pack_start(*labRegionShowMask, Gtk::PACK_SHRINK, 4);

    pack_start(*labRegionBox, Gtk::PACK_EXPAND_WIDGET, 4);

    labRegionSaturation->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));
    labRegionSlope->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));
    labRegionOffset->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));
    labRegionPower->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));
    labRegionMaskBlur->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));
    //------------------------------------------------------------------------

    show_all();

    disableListener();
    methodChanged();
    enableListener();
}

ColorToning::~ColorToning()
{
    idle_register.destroy();

    delete colorCurveEditorG;
    delete opacityCurveEditorG;
    delete clCurveEditorG;
    delete cl2CurveEditorG;
}


void ColorToning::setListener(ToolPanelListener *tpl)
{
    ToolPanel::setListener(tpl);
    labgrid->setListener(tpl);
    labRegionAB->setListener(tpl);
}

/*
void ColorToning::neutralCurves_pressed () {
    disableListener();

    bool changed = false;
    changed |= colorShape->reset();
    changed |= opacityShape->reset();
    changed |= clshape->reset();
    changed |= cl2shape->reset();

    enableListener();

    if (listener && enabled->get_active() && changed)
        listener->panelChanged (EvColorToningNeutralcur, M("GENERAL_RESET"));
}
*/

// Will only reset the channel mixer
void ColorToning::neutral_pressed ()
{
    disableListener();
    redlow->resetValue(false);
    greenlow->resetValue(false);
    bluelow->resetValue(false);
    redmed->resetValue(false);
    greenmed->resetValue(false);
    bluemed->resetValue(false);
    redhigh->resetValue(false);
    greenhigh->resetValue(false);
    bluehigh->resetValue(false);
    //balance->resetValue(false);

    enableListener();

    if (listener && getEnabled()) {
        listener->panelChanged (EvColorToningNeutral, M("GENERAL_RESET"));
    }
}

void ColorToning::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    methodconn.block(true);
    twocconn.block(true);
    colorShape->setCurve (pp->colorToning.colorCurve);
    opacityShape->setCurve (pp->colorToning.opacityCurve);
    clshape->setCurve (pp->colorToning.clcurve);
    cl2shape->setCurve (pp->colorToning.cl2curve);

    labRegionData = pp->colorToning.labregions;
    if (labRegionData.empty()) {
        labRegionData.emplace_back(rtengine::procparams::ColorToningParams::LabCorrectionRegion());
    }
    if (pp->colorToning.labregionsShowMask >= 0) {
        labRegionSelected = pp->colorToning.labregionsShowMask;
        labRegionShowMask->set_active(true);
    } else {
        labRegionSelected = 0;
        labRegionShowMask->set_active(false);
    }
    labRegionPopulateList();
    labRegionShow(labRegionSelected);

    if (pedited) {
        redlow->setEditedState (pedited->colorToning.redlow ? Edited : UnEdited);
        greenlow->setEditedState (pedited->colorToning.greenlow ? Edited : UnEdited);
        bluelow->setEditedState (pedited->colorToning.bluelow ? Edited : UnEdited);
        balance->setEditedState (pedited->colorToning.balance ? Edited : UnEdited);
        redmed->setEditedState (pedited->colorToning.redmed ? Edited : UnEdited);
        greenmed->setEditedState (pedited->colorToning.greenmed ? Edited : UnEdited);
        bluemed->setEditedState (pedited->colorToning.bluemed ? Edited : UnEdited);
        redhigh->setEditedState (pedited->colorToning.redhigh ? Edited : UnEdited);
        greenhigh->setEditedState (pedited->colorToning.greenhigh ? Edited : UnEdited);
        bluehigh->setEditedState (pedited->colorToning.bluehigh ? Edited : UnEdited);

        hlColSat->setEditedState (pedited->colorToning.hlColSat ? Edited : UnEdited);
        shadowsColSat->setEditedState (pedited->colorToning.shadowsColSat ? Edited : UnEdited);

        set_inconsistent (multiImage && !pedited->colorToning.enabled);
        colorShape->setUnChanged (!pedited->colorToning.colorCurve);
        opacityShape->setUnChanged (!pedited->colorToning.opacityCurve);
        autosat->set_inconsistent (!pedited->colorToning.autosat);
        clshape->setUnChanged  (!pedited->colorToning.clcurve);
        cl2shape->setUnChanged  (!pedited->colorToning.cl2curve);
        lumamode->set_inconsistent (!pedited->colorToning.lumamode);

        labgrid->setEdited(pedited->colorToning.labgridALow || pedited->colorToning.labgridBLow || pedited->colorToning.labgridAHigh || pedited->colorToning.labgridBHigh);

        labRegionAB->setEdited(pedited->colorToning.labregions);
        labRegionShowMask->set_inconsistent(!pedited->colorToning.labregionsShowMask);
    }

    redlow->setValue    (pp->colorToning.redlow);
    greenlow->setValue    (pp->colorToning.greenlow);
    bluelow->setValue    (pp->colorToning.bluelow);
    balance->setValue    (pp->colorToning.balance);
    redmed->setValue    (pp->colorToning.redmed);
    greenmed->setValue    (pp->colorToning.greenmed);
    bluemed->setValue    (pp->colorToning.bluemed);
    redhigh->setValue    (pp->colorToning.redhigh);
    greenhigh->setValue    (pp->colorToning.greenhigh);
    bluehigh->setValue    (pp->colorToning.bluehigh);

    setEnabled (pp->colorToning.enabled);

    autosatConn.block (true);
    autosat->set_active (pp->colorToning.autosat);
    autosatConn.block (false);
    lastautosat = pp->colorToning.autosat;

    satProtectionThreshold->setValue (pp->colorToning.satProtectionThreshold);
    saturatedOpacity->setValue (pp->colorToning.saturatedOpacity);
    hlColSat->setValue<int> (pp->colorToning.hlColSat);
    shadowsColSat->setValue<int> (pp->colorToning.shadowsColSat);
    strength->setValue (pp->colorToning.strength);
    lumamodeConn.block (true);
    lumamode->set_active (pp->colorToning.lumamode);
    lumamodeConn.block (false);

    lastLumamode = pp->colorToning.lumamode;

    labgrid->setParams(pp->colorToning.labgridALow / ColorToningParams::LABGRID_CORR_MAX, pp->colorToning.labgridBLow / ColorToningParams::LABGRID_CORR_MAX, pp->colorToning.labgridAHigh / ColorToningParams::LABGRID_CORR_MAX, pp->colorToning.labgridBHigh / ColorToningParams::LABGRID_CORR_MAX, 0, 0, 0, 0, 0, 0, false);

    if (pedited && !pedited->colorToning.method) {
        method->set_active (7);
    } else if (pp->colorToning.method == "Lab") {
        method->set_active (0);
    } else if (pp->colorToning.method == "RGBSliders") {
        method->set_active (1);
    } else if (pp->colorToning.method == "RGBCurves") {
        method->set_active (2);
    } else if (pp->colorToning.method == "Splitco") {
        method->set_active (3);
    } else if (pp->colorToning.method == "Splitlr") {
        method->set_active (4);
    } else if (pp->colorToning.method == "LabGrid") {
        method->set_active(5);
    } else if (pp->colorToning.method == "LabRegions") {
        method->set_active(6);
    }

    methodChanged();
    methodconn.block(false);


    if (pedited && !pedited->colorToning.twocolor) {
        twocolor->set_active (4);
    } else if (pp->colorToning.twocolor == "Std") {
        twocolor->set_active (0);
    } else if (pp->colorToning.twocolor == "All") {
        twocolor->set_active (1);
    } else if (pp->colorToning.twocolor == "Separ") {
        twocolor->set_active (2);
    } else if (pp->colorToning.twocolor == "Two") {
        twocolor->set_active (3);
    }

    twocolorChanged(true);

    twocconn.block(false);

    enableListener ();
}

void ColorToning::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->colorToning.redlow    = redlow->getValue ();
    pp->colorToning.greenlow  = greenlow->getValue ();
    pp->colorToning.bluelow   = bluelow->getValue ();
    pp->colorToning.balance   = balance->getIntValue ();
    pp->colorToning.redmed    = redmed->getValue ();
    pp->colorToning.greenmed  = greenmed->getValue ();
    pp->colorToning.bluemed   = bluemed->getValue ();
    pp->colorToning.redhigh   = redhigh->getValue ();
    pp->colorToning.greenhigh = greenhigh->getValue ();
    pp->colorToning.bluehigh  = bluehigh->getValue ();

    pp->colorToning.enabled      = getEnabled();
    pp->colorToning.colorCurve   = colorShape->getCurve ();
    pp->colorToning.opacityCurve = opacityShape->getCurve ();
    pp->colorToning.clcurve      = clshape->getCurve ();
    pp->colorToning.cl2curve     = cl2shape->getCurve ();
    pp->colorToning.lumamode     = lumamode->get_active();

    pp->colorToning.hlColSat               = hlColSat->getValue<int> ();
    pp->colorToning.shadowsColSat          = shadowsColSat->getValue<int> ();
    pp->colorToning.autosat                = autosat->get_active();
    pp->colorToning.satProtectionThreshold = satProtectionThreshold->getIntValue();
    pp->colorToning.saturatedOpacity       = saturatedOpacity->getIntValue();
    pp->colorToning.strength               = strength->getIntValue();
    double zerox = 0.;
    double zeroy = 0.;

    labgrid->getParams(pp->colorToning.labgridALow, pp->colorToning.labgridBLow, pp->colorToning.labgridAHigh, pp->colorToning.labgridBHigh, zerox, zeroy, zerox, zeroy, zerox, zeroy);
    pp->colorToning.labgridALow *= ColorToningParams::LABGRID_CORR_MAX;
    pp->colorToning.labgridAHigh *= ColorToningParams::LABGRID_CORR_MAX;
    pp->colorToning.labgridBLow *= ColorToningParams::LABGRID_CORR_MAX;
    pp->colorToning.labgridBHigh *= ColorToningParams::LABGRID_CORR_MAX;

    labRegionGet(labRegionSelected);
    labRegionShow(labRegionSelected, true);
    pp->colorToning.labregions = labRegionData;
    if (labRegionShowMask->get_active()) {
        pp->colorToning.labregionsShowMask = labRegionSelected;
    } else {
        pp->colorToning.labregionsShowMask = -1;
    }

    if (pedited) {
        pedited->colorToning.redlow     = redlow->getEditedState ();
        pedited->colorToning.greenlow   = greenlow->getEditedState ();
        pedited->colorToning.bluelow    = bluelow->getEditedState ();
        pedited->colorToning.balance    = balance->getEditedState ();
        pedited->colorToning.redmed     = redmed->getEditedState ();
        pedited->colorToning.greenmed   = greenmed->getEditedState ();
        pedited->colorToning.bluemed    = bluemed->getEditedState ();
        pedited->colorToning.redhigh    = redhigh->getEditedState ();
        pedited->colorToning.greenhigh  = greenhigh->getEditedState ();
        pedited->colorToning.bluehigh   = bluehigh->getEditedState ();
        pedited->colorToning.method     = method->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->colorToning.twocolor   = twocolor->get_active_text() != M("GENERAL_UNCHANGED");

        pedited->colorToning.enabled       = !get_inconsistent();
        pedited->colorToning.autosat       = !autosat->get_inconsistent();
        pedited->colorToning.colorCurve    = !colorShape->isUnChanged ();
        pedited->colorToning.opacityCurve  = !opacityShape->isUnChanged ();
        pedited->colorToning.clcurve       = !clshape->isUnChanged ();
        pedited->colorToning.cl2curve      = !cl2shape->isUnChanged ();
        pedited->colorToning.lumamode      = !lumamode->get_inconsistent();

        pedited->colorToning.hlColSat      = hlColSat->getEditedState ();
        pedited->colorToning.shadowsColSat = shadowsColSat->getEditedState ();

        pedited->colorToning.labgridALow = pedited->colorToning.labgridBLow = pedited->colorToning.labgridAHigh = pedited->colorToning.labgridBHigh = labgrid->getEdited();

        pedited->colorToning.labregions = labRegionAB->getEdited();
        pedited->colorToning.labregionsShowMask = !labRegionShowMask->get_inconsistent();
    }

    if (method->get_active_row_number() == 0) {
        pp->colorToning.method = "Lab";
    } else if (method->get_active_row_number() == 1) {
        pp->colorToning.method = "RGBSliders";
    } else if (method->get_active_row_number() == 2) {
        pp->colorToning.method = "RGBCurves";
    } else if (method->get_active_row_number() == 3) {
        pp->colorToning.method = "Splitco";
    } else if (method->get_active_row_number() == 4) {
        pp->colorToning.method = "Splitlr";
    } else if (method->get_active_row_number() == 5) {
        pp->colorToning.method = "LabGrid";
    } else if (method->get_active_row_number() == 6) {
        pp->colorToning.method = "LabRegions";
    }

    if (twocolor->get_active_row_number() == 0) {
        pp->colorToning.twocolor = "Std";
    } else if (twocolor->get_active_row_number() == 1) {
        pp->colorToning.twocolor = "All";
    } else if (twocolor->get_active_row_number() == 2) {
        pp->colorToning.twocolor = "Separ";
    } else if (twocolor->get_active_row_number() == 3) {
        pp->colorToning.twocolor = "Two";
    }
}

void ColorToning::lumamodeChanged ()
{

    if (batchMode) {
        if (lumamode->get_inconsistent()) {
            lumamode->set_inconsistent (false);
            lumamodeConn.block (true);
            lumamode->set_active (false);
            lumamodeConn.block (false);
        } else if (lastLumamode) {
            lumamode->set_inconsistent (true);
        }

        lastLumamode = lumamode->get_active ();
    }

    if (listener && getEnabled()) {
        if (lumamode->get_active ()) {
            listener->panelChanged (EvColorToningLumamode, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvColorToningLumamode, M("GENERAL_DISABLED"));
        }
    }
}

void ColorToning::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    redlow->setDefault (defParams->colorToning.redlow);
    greenlow->setDefault (defParams->colorToning.greenlow);
    bluelow->setDefault (defParams->colorToning.bluelow);
    balance->setDefault (defParams->colorToning.balance);
    redmed->setDefault (defParams->colorToning.redmed);
    greenmed->setDefault (defParams->colorToning.greenmed);
    bluemed->setDefault (defParams->colorToning.bluemed);
    redhigh->setDefault (defParams->colorToning.redhigh);
    greenhigh->setDefault (defParams->colorToning.greenhigh);
    bluehigh->setDefault (defParams->colorToning.bluehigh);
    satProtectionThreshold->setDefault (defParams->colorToning.satProtectionThreshold);
    saturatedOpacity->setDefault (defParams->colorToning.saturatedOpacity);
    hlColSat->setDefault<int> (defParams->colorToning.hlColSat);
    shadowsColSat->setDefault<int> (defParams->colorToning.shadowsColSat);
    strength->setDefault (defParams->colorToning.strength);
    labgrid->setDefault(defParams->colorToning.labgridALow / ColorToningParams::LABGRID_CORR_MAX, defParams->colorToning.labgridBLow / ColorToningParams::LABGRID_CORR_MAX, defParams->colorToning.labgridAHigh / ColorToningParams::LABGRID_CORR_MAX, defParams->colorToning.labgridBHigh / ColorToningParams::LABGRID_CORR_MAX, 0, 0, 0, 0, 0, 0);


    if (pedited) {
        redlow->setDefaultEditedState (pedited->colorToning.redlow ? Edited : UnEdited);
        greenlow->setDefaultEditedState (pedited->colorToning.greenlow ? Edited : UnEdited);
        bluelow->setDefaultEditedState (pedited->colorToning.bluelow ? Edited : UnEdited);
        balance->setDefaultEditedState (pedited->colorToning.balance ? Edited : UnEdited);
        redmed->setDefaultEditedState (pedited->colorToning.redmed ? Edited : UnEdited);
        greenmed->setDefaultEditedState (pedited->colorToning.greenmed ? Edited : UnEdited);
        bluemed->setDefaultEditedState (pedited->colorToning.bluemed ? Edited : UnEdited);
        redhigh->setDefaultEditedState (pedited->colorToning.redhigh ? Edited : UnEdited);
        greenhigh->setDefaultEditedState (pedited->colorToning.greenhigh ? Edited : UnEdited);
        bluehigh->setDefaultEditedState (pedited->colorToning.bluehigh ? Edited : UnEdited);
        satProtectionThreshold->setDefaultEditedState (pedited->colorToning.satprotectionthreshold ? Edited : UnEdited);
        saturatedOpacity->setDefaultEditedState (pedited->colorToning.saturatedopacity ? Edited : UnEdited);
        hlColSat->setDefaultEditedState (pedited->colorToning.hlColSat ? Edited : UnEdited);
        shadowsColSat->setDefaultEditedState (pedited->colorToning.shadowsColSat ? Edited : UnEdited);
        strength->setDefaultEditedState (pedited->colorToning.strength ? Edited : UnEdited);
        labgrid->setEdited((pedited->colorToning.labgridALow || pedited->colorToning.labgridBLow || pedited->colorToning.labgridAHigh || pedited->colorToning.labgridBHigh) ? Edited : UnEdited);

        labRegionAB->setEdited(pedited->colorToning.labregions ? Edited : UnEdited);
    } else {
        redlow->setDefaultEditedState (Irrelevant);
        greenlow->setDefaultEditedState (Irrelevant);
        bluelow->setDefaultEditedState (Irrelevant);
        balance->setDefaultEditedState (Irrelevant);
        redmed->setDefaultEditedState (Irrelevant);
        greenmed->setDefaultEditedState (Irrelevant);
        bluemed->setDefaultEditedState (Irrelevant);
        redhigh->setDefaultEditedState (Irrelevant);
        greenhigh->setDefaultEditedState (Irrelevant);
        bluehigh->setDefaultEditedState (Irrelevant);
        satProtectionThreshold->setDefaultEditedState (Irrelevant);
        saturatedOpacity->setDefaultEditedState (Irrelevant);
        hlColSat->setDefaultEditedState (Irrelevant);
        shadowsColSat->setDefaultEditedState (Irrelevant);
        strength->setDefaultEditedState (Irrelevant);
        labgrid->setEdited(Edited);
        labRegionAB->setEdited(Edited);
    }
}

void ColorToning::setAdjusterBehavior (bool splitAdd, bool satThresholdAdd, bool satOpacityAdd, bool strprotectAdd, bool balanceAdd)
{
    redlow->setAddMode(splitAdd);
    greenlow->setAddMode(splitAdd);
    bluelow->setAddMode(splitAdd);
    balance->setAddMode(splitAdd);
    redmed->setAddMode(splitAdd);
    greenmed->setAddMode(splitAdd);
    bluemed->setAddMode(splitAdd);
    redhigh->setAddMode(splitAdd);
    greenhigh->setAddMode(splitAdd);
    bluehigh->setAddMode(splitAdd);
    satProtectionThreshold->setAddMode(satThresholdAdd);
    saturatedOpacity->setAddMode(satOpacityAdd);
    balance->setAddMode(balanceAdd);
    strength->setAddMode(strprotectAdd);

}

void ColorToning::autoColorTonChanged(int bwct, int satthres, int satprot)
{
    nextbw = bwct;
    nextsatth = satthres;
    nextsatpr = satprot;

    idle_register.add(
        [this]() -> bool
        {
            disableListener();
            saturatedOpacity->setValue(nextsatpr);
            satProtectionThreshold->setValue(nextsatth);
            enableListener();
            return false;
        }
    );
}

void ColorToning::adjusterChanged (ThresholdAdjuster* a, double newBottom, double newTop)
{
    if (listener && getEnabled()) {
        listener->panelChanged(
            a == hlColSat
                ? EvColorToningHighights
                : EvColorToningShadows,
            Glib::ustring::compose(Glib::ustring(M("TP_COLORTONING_HUE") + ": %1" + "\n" + M("TP_COLORTONING_STRENGTH") + ": %2"), int(newTop), int(newBottom))
        );
    }
}

void ColorToning::adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight)
{
}

void ColorToning::adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop)
{
}

void ColorToning::adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight)
{
}

void ColorToning::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
}

//Two Color changed
void ColorToning::twocolorChanged (bool changedbymethod)
{
    if (!batchMode) {
        if(method->get_active_row_number() == 0) {              // Lab
            if(twocolor->get_active_row_number() == 0) {
                colorCurveEditorG->show();   // visible
                opacityCurveEditorG->show(); // visible
                clCurveEditorG->hide();
                cl2CurveEditorG->hide();
            } else if(twocolor->get_active_row_number() == 1  || twocolor->get_active_row_number() == 3) {
                colorCurveEditorG->show();   // visible
                opacityCurveEditorG->hide();
                clCurveEditorG->show();      // visible
                cl2CurveEditorG->hide();
                irg->hide();

            } else if(twocolor->get_active_row_number() == 2) {
                colorCurveEditorG->show(); // visible
                opacityCurveEditorG->hide();
                clCurveEditorG->show();      // visible
                cl2CurveEditorG->show();     // visible
                irg->show();
            }
        } else if(method->get_active_row_number() == 1) {       // RGB Sliders
            colorCurveEditorG->hide();
            opacityCurveEditorG->hide();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
        } else if(method->get_active_row_number() == 2) {       // RGB Curves
            colorCurveEditorG->show();       // visible
            opacityCurveEditorG->show();     // visible
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
        } else if(method->get_active_row_number() == 3) {       // Split LR
            colorCurveEditorG->hide();
            opacityCurveEditorG->hide();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
        } else if(method->get_active_row_number() == 4) {       // Split color
            colorCurveEditorG->hide();
            opacityCurveEditorG->hide();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
        }
    }

    if (listener && getEnabled() && !changedbymethod) {
        listener->panelChanged (EvColorToningTwocolor, twocolor->get_active_text ());
    }
}

void ColorToning::twoColorChangedByGui()
{
    twocolorChanged(false);
}

void ColorToning::methodChanged ()
{

    if (!batchMode) {
        labgrid->hide();
        labRegionBox->hide();

        if (method->get_active_row_number() == 0) { // Lab
            colorSep->show();
            colorCurveEditorG->show();
            twocolor->show();
            opacityCurveEditorG->hide();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
            //neutralCurves->show();
            hlColSat->hide();
            shadowsColSat->hide();
            balance->hide();

//          satLimiterSep->show();
            if(autosat->get_active()) {
                saturatedOpacity->set_sensitive(false);
                satProtectionThreshold->set_sensitive(false);
                satProtectionThreshold->show();
                saturatedOpacity->show();
            } else {
                satProtectionThreshold->show();
                saturatedOpacity->show();
                saturatedOpacity->set_sensitive(true);
                satProtectionThreshold->set_sensitive(true);
            }

            autosat->show();
            p1Frame->show();

            strength->hide();
            chanMixerBox->hide();
            neutrHBox->hide();
            lumamode->hide();
            //splitSep->hide();
            //satlow->hide();
            //sathigh->hide();

            twocolorChanged(true);
        } else if (method->get_active_row_number() == 1) { // RGB Sliders
            colorSep->hide();
            colorCurveEditorG->hide();
            twocolor->hide();
            twocolor->set_active (false);
            opacityCurveEditorG->hide();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
            //neutralCurves->hide();
            hlColSat->show();
            shadowsColSat->show();
            balance->show();
//          satLimiterSep->show();
            autosat->show();
            p1Frame->show();

            if(autosat->get_active()) {
                saturatedOpacity->set_sensitive(false);
                satProtectionThreshold->set_sensitive(false);
                satProtectionThreshold->show();
                saturatedOpacity->show();
            } else {
                satProtectionThreshold->show();
                saturatedOpacity->show();
                saturatedOpacity->set_sensitive(true);
                satProtectionThreshold->set_sensitive(true);
            }

            strength->hide();
            chanMixerBox->hide();
            neutrHBox->hide();
            lumamode->hide();

        } else if (method->get_active_row_number() == 2) { // RGB Curves
            colorSep->hide();
            colorCurveEditorG->show();
            twocolor->hide();
            twocolor->set_active (false);
            opacityCurveEditorG->show();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
            //neutralCurves->show();
            hlColSat->hide();
            shadowsColSat->hide();
            balance->hide();
//          satLimiterSep->show();
            p1Frame->show();

            autosat->show();

            if(autosat->get_active()) {
                saturatedOpacity->set_sensitive(false);
                satProtectionThreshold->set_sensitive(false);
                satProtectionThreshold->show();
                saturatedOpacity->show();
            } else {
                satProtectionThreshold->show();
                saturatedOpacity->show();
                saturatedOpacity->set_sensitive(true);
                satProtectionThreshold->set_sensitive(true);
            }

            strength->hide();
            chanMixerBox->hide();
            neutrHBox->hide();
            lumamode->hide();
        } else if (method->get_active_row_number() == 3) { // Split LR
            colorSep->hide();
            colorCurveEditorG->hide();
            twocolor->hide();
            twocolor->set_active (false);
            opacityCurveEditorG->hide();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
            //neutralCurves->hide();
            hlColSat->hide();
            shadowsColSat->hide();
            balance->hide();
            p1Frame->hide();

//          satLimiterSep->hide();
            autosat->hide();
            satProtectionThreshold->hide();
            saturatedOpacity->hide();
            strength->show();
            chanMixerBox->show();
            neutrHBox->show();
            lumamode->show();
        } else if (method->get_active_row_number() == 4) { // Split Color
            colorSep->show();
            colorCurveEditorG->hide();
            twocolor->hide();
            opacityCurveEditorG->hide();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
            //neutralCurves->hide();
            hlColSat->show();
            shadowsColSat->show();
            balance->show();
//          satLimiterSep->hide();
            p1Frame->hide();

            autosat->hide();
            satProtectionThreshold->hide();
            saturatedOpacity->hide();
            strength->show();

            chanMixerBox->hide();
            neutrHBox->hide();
            lumamode->show();
        } else if (method->get_active_row_number() == 5 || method->get_active_row_number() == 6) { // Lab Grid or Lab Regions
            colorSep->hide();
            colorCurveEditorG->hide();
            twocolor->hide();
            opacityCurveEditorG->hide();
            clCurveEditorG->hide();
            cl2CurveEditorG->hide();
            hlColSat->hide();
            shadowsColSat->hide();
            balance->hide();
            p1Frame->hide();
            autosat->hide();
            satProtectionThreshold->hide();
            saturatedOpacity->hide();
            strength->hide();
            chanMixerBox->hide();
            neutrHBox->hide();
            lumamode->hide();

            if (method->get_active_row_number() == 5) {
                labgrid->show();
            } else {
                labRegionBox->show();
            }
        }
    }

    if (listener && getEnabled()) {
        listener->panelChanged (EvColorToningMethod, method->get_active_text ());
    }
}


void ColorToning::autoOpenCurve  ()
{
    colorShape->openIfNonlinear();
    opacityShape->openIfNonlinear();
    clshape->openIfNonlinear();
    cl2shape->openIfNonlinear();
}

void ColorToning::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R = 0.f, G = 0.f, B = 0.f;

    if (callerId == 1) {         // opacity curve left bar(s)
        Color::hsv2rgb01(float(valY*0.8), 1.0f, 0.5f, R, G, B);
    } else if (callerId == 2) {  // Slider 1 background
        if (valY <= 0.5)
            // the hue range
        {
            Color::hsv2rgb01(float(valX), 1.0f, 0.5f, R, G, B);
        } else {
            // the strength applied to the current hue
            double strength, hue;
            hlColSat->getValue(strength, hue);
            Color::hsv2rgb01(hue / 360.0, 1.f, 1.f, R, G, B);
            const double gray = 0.46;
            R = (gray * (1.0 - valX)) + static_cast<double>(R) * valX;
            G = (gray * (1.0 - valX)) + static_cast<double>(G) * valX;
            B = (gray * (1.0 - valX)) + static_cast<double>(B) * valX;
        }
    } else if (callerId == 3) {  // Slider 2 background
        if (valY <= 0.5)
            // the hue range
        {
            Color::hsv2rgb01(float(valX), 1.0f, 0.5f, R, G, B);
        } else {
            // the strength applied to the current hue
            double strength, hue;
            shadowsColSat->getValue(strength, hue);
            Color::hsv2rgb01(hue / 360.0, 1.f, 1.f, R, G, B);
            const double gray = 0.46;
            R = (gray * (1.0 - valX)) + static_cast<double>(R) * valX;
            G = (gray * (1.0 - valX)) + static_cast<double>(G) * valX;
            B = (gray * (1.0 - valX)) + static_cast<double>(B) * valX;
        }
    } else if (callerId == 4) {  // color curve vertical and horizontal crosshair
        Color::hsv2rgb01(float(valY), 1.0f, 0.5f, R, G, B);
    } else if (callerId == ID_LABREGION_HUE) {
        // TODO
        float x = valX - 1.0/6.0;
        if (x < 0.f) {
            x += 1.f;
        }
        x = log2lin(x, 3.f);
        // float x = valX;
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
    } else if (callerId == ID_LABREGION_HUE+1) {
        Color::hsv2rgb01(float(valY), float(valX), 0.5f, R, G, B);
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
}

void ColorToning::curveChanged (CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == colorShape) {
            listener->panelChanged (EvColorToningColor, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == opacityShape) {
            listener->panelChanged (EvColorToningOpacity, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == clshape) {
            listener->panelChanged (EvColorToningCLCurve, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == cl2shape) {
            listener->panelChanged (EvColorToningLLCurve, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == labRegionHueMask) {
            listener->panelChanged(EvLabRegionHueMask, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == labRegionChromaticityMask) {
            listener->panelChanged(EvLabRegionChromaticityMask, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == labRegionLightnessMask) {
            listener->panelChanged(EvLabRegionLightnessMask, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void ColorToning::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvColorToningEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvColorToningEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvColorToningEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void ColorToning::autosatChanged ()
{

    if (batchMode) {
        if (autosat->get_inconsistent()) {
            autosat->set_inconsistent (false);
            autosatConn.block (true);
            autosat->set_active (false);
            autosatConn.block (false);
        } else if (lastautosat) {
            autosat->set_inconsistent (true);
        }

        lastautosat = autosat->get_active ();
    }

    if (listener) {
        if (autosat->get_active()) {
            if (getEnabled()) {
                listener->panelChanged (EvColorToningautosat, M("GENERAL_ENABLED"));
            }

            saturatedOpacity->set_sensitive(false);
            satProtectionThreshold->set_sensitive(false);
        } else {
            if (getEnabled()) {
                listener->panelChanged (EvColorToningautosat, M("GENERAL_DISABLED"));
            }

            saturatedOpacity->set_sensitive(true);
            satProtectionThreshold->set_sensitive(true);
        }

    }
}

void ColorToning::trimValues (rtengine::procparams::ProcParams* pp)
{

    redlow->trimValue(pp->colorToning.redlow);
    greenlow->trimValue(pp->colorToning.greenlow);
    bluelow->trimValue(pp->colorToning.bluelow);
    balance->trimValue(pp->colorToning.balance);
    redmed->trimValue(pp->colorToning.redmed);
    greenmed->trimValue(pp->colorToning.greenmed);
    bluemed->trimValue(pp->colorToning.bluemed);
    redhigh->trimValue(pp->colorToning.redhigh);
    greenhigh->trimValue(pp->colorToning.greenhigh);
    bluehigh->trimValue(pp->colorToning.bluehigh);
}

void ColorToning::adjusterChanged(Adjuster* a, double newval)
{
    if (!listener || !getEnabled()) {
        return;
    }

    if (a == redlow) {
        listener->panelChanged (EvColorToningredlow, redlow->getTextValue());
    } else if (a == greenlow) {
        listener->panelChanged (EvColorToninggreenlow, greenlow->getTextValue());
    } else if (a == bluelow) {
        listener->panelChanged (EvColorToningbluelow, bluelow->getTextValue());
    } else if (a == redmed) {
        listener->panelChanged (EvColorToningredmed, redmed->getTextValue());
    } else if (a == greenmed) {
        listener->panelChanged (EvColorToninggreenmed, greenmed->getTextValue());
    } else if (a == bluemed) {
        listener->panelChanged (EvColorToningbluemed, bluemed->getTextValue());
    } else if (a == redhigh) {
        listener->panelChanged (EvColorToningredhigh, redhigh->getTextValue());
    } else if (a == greenhigh) {
        listener->panelChanged (EvColorToninggreenhigh, greenhigh->getTextValue());
    } else if (a == bluehigh) {
        listener->panelChanged (EvColorToningbluehigh, bluehigh->getTextValue());
    } else if (a == balance) {
        listener->panelChanged (EvColorToningbalance, balance->getTextValue());
    } else if (a == satProtectionThreshold) {
        listener->panelChanged (EvColorToningSatThreshold, a->getTextValue());
    } else if (a == saturatedOpacity) {
        listener->panelChanged (EvColorToningSatProtection, a->getTextValue());
    } else if (a == strength) {
        listener->panelChanged (EvColorToningStrength, a->getTextValue());
    } else if (a == labRegionSaturation) {
        listener->panelChanged(EvLabRegionSaturation, a->getTextValue());
    } else if (a == labRegionSlope) {
        listener->panelChanged(EvLabRegionSlope, a->getTextValue());
    } else if (a == labRegionOffset) {
        listener->panelChanged(EvLabRegionOffset, a->getTextValue());
    } else if (a == labRegionPower) {
        listener->panelChanged(EvLabRegionPower, a->getTextValue());
    } else if (a == labRegionMaskBlur) {
        listener->panelChanged(EvLabRegionMaskBlur, a->getTextValue());
    }
}

void ColorToning::setBatchMode (bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);
    method->append (M("GENERAL_UNCHANGED"));
    twocolor->append (M("GENERAL_UNCHANGED"));
    hlColSat->showEditedCB ();
    shadowsColSat->showEditedCB ();
    redlow->showEditedCB ();
    greenlow->showEditedCB ();
    bluelow->showEditedCB ();
    balance->showEditedCB ();
    redmed->showEditedCB ();
    greenmed->showEditedCB ();
    bluemed->showEditedCB ();
    redhigh->showEditedCB ();
    greenhigh->showEditedCB ();
    bluehigh->showEditedCB ();

    colorCurveEditorG->setBatchMode (batchMode);
    opacityCurveEditorG->setBatchMode (batchMode);
    clCurveEditorG->setBatchMode (batchMode);
    cl2CurveEditorG->setBatchMode (batchMode);

}


void ColorToning::onLabRegionSelectionChanged()
{
    auto s = labRegionList->get_selected();
    if (!s.empty()) {
        // update the selected values
        labRegionGet(labRegionSelected);
        labRegionSelected = s[0];
        labRegionShow(labRegionSelected);
        if (labRegionShowMask->get_active()) {
            labRegionShowMaskChanged();
        }
    }
}


void ColorToning::labRegionGet(int idx)
{
    if (idx < 0 || size_t(idx) >= labRegionData.size()) {
        return;
    }

    auto &r = labRegionData[idx];
    double la, lb;
    double zerox = 0.;
    double zeroy = 0.;
    labRegionAB->getParams(la, lb, r.a, r.b, zerox, zeroy, zerox, zeroy, zerox, zeroy);
    r.saturation = labRegionSaturation->getValue();
    r.slope = labRegionSlope->getValue();
    r.offset = labRegionOffset->getValue();
    r.power = labRegionPower->getValue();
    r.hueMask = labRegionHueMask->getCurve();
    r.chromaticityMask = labRegionChromaticityMask->getCurve();
    r.lightnessMask = labRegionLightnessMask->getCurve();
    r.maskBlur = labRegionMaskBlur->getValue();
    r.channel = labRegionChannel->get_active_row_number() - 1;
}


void ColorToning::labRegionAddPressed()
{
    labRegionSelected = labRegionData.size();
    labRegionData.push_back(rtengine::procparams::ColorToningParams::LabCorrectionRegion());
    labRegionPopulateList();
    labRegionShow(labRegionSelected);

    if (listener) {
        listener->panelChanged(EvLabRegionList, M("HISTORY_CHANGED"));
    }
}


void ColorToning::labRegionRemovePressed()
{
    if (labRegionList->size() > 1) {
        labRegionData.erase(labRegionData.begin() + labRegionSelected);
        labRegionSelected = LIM(labRegionSelected-1, 0, int(labRegionData.size()-1));
        labRegionPopulateList();
        labRegionShow(labRegionSelected);

        if (listener) {
            listener->panelChanged(EvLabRegionList, M("HISTORY_CHANGED"));
        }
    }
}


void ColorToning::labRegionUpPressed()
{
    if (labRegionSelected > 0) {
        auto r = labRegionData[labRegionSelected];
        labRegionData.erase(labRegionData.begin() + labRegionSelected);
        --labRegionSelected;
        labRegionData.insert(labRegionData.begin() + labRegionSelected, r);
        labRegionPopulateList();

        if (listener) {
            listener->panelChanged(EvLabRegionList, M("HISTORY_CHANGED"));
        }
    }
}


void ColorToning::labRegionDownPressed()
{
    if (labRegionSelected < int(labRegionData.size()-1)) {
        auto r = labRegionData[labRegionSelected];
        labRegionData.erase(labRegionData.begin() + labRegionSelected);
        ++labRegionSelected;
        labRegionData.insert(labRegionData.begin() + labRegionSelected, r);
        labRegionPopulateList();

        if (listener) {
            listener->panelChanged(EvLabRegionList, M("HISTORY_CHANGED"));
        }
    }
}


void ColorToning::labRegionCopyPressed()
{
    if (labRegionSelected < int(labRegionData.size())) {
        auto r = labRegionData[labRegionSelected];
        labRegionData.push_back(r);
        labRegionSelected = labRegionData.size()-1;
        labRegionPopulateList();

        if (listener) {
            listener->panelChanged(EvLabRegionList, M("HISTORY_CHANGED"));
        }
    }
}


void ColorToning::labRegionShowMaskChanged()
{
    if (listener) {
        listener->panelChanged(EvLabRegionShowMask, labRegionShowMask->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void ColorToning::labRegionPopulateList()
{
    ConnectionBlocker b(labRegionSelectionConn);
    labRegionList->clear_items();
    rtengine::procparams::ColorToningParams::LabCorrectionRegion dflt;

    for (size_t i = 0; i < labRegionData.size(); ++i) {
        auto &r = labRegionData[i];
        auto j = labRegionList->append(std::to_string(i+1));
        labRegionList->set_text(j, 1, Glib::ustring::compose("a=%1 b=%2 S=%3\ns=%4 o=%5 p=%6", round_ab(r.a), round_ab(r.b), r.saturation, r.slope, r.offset, r.power));
        const char *ch = "";
        switch (r.channel) {
        case rtengine::procparams::ColorToningParams::LabCorrectionRegion::CHAN_R:
            ch = "\n[Red]"; break;
        case rtengine::procparams::ColorToningParams::LabCorrectionRegion::CHAN_G:
            ch = "\n[Green]"; break;
        case rtengine::procparams::ColorToningParams::LabCorrectionRegion::CHAN_B:
            ch = "\n[Blue]"; break;
        default:
            ch = "";
        }
        labRegionList->set_text(
            j, 2, Glib::ustring::compose(
                "%1%2%3%4%5",
                hasMask(dflt.hueMask, r.hueMask) ? "H" : "",
                hasMask(dflt.chromaticityMask, r.chromaticityMask) ? "C" : "",
                hasMask(dflt.lightnessMask, r.lightnessMask) ? "L" : "",
                r.maskBlur ? Glib::ustring::compose(" b=%1", r.maskBlur) : "",
                ch));
    }
}


void ColorToning::labRegionShow(int idx, bool list_only)
{
    const bool disable = listener;
    if (disable) {
        disableListener();
    }
    rtengine::procparams::ColorToningParams::LabCorrectionRegion dflt;
    auto &r = labRegionData[idx];
    if (!list_only) {
        labRegionAB->setParams(0, 0, r.a, r.b,0, 0, 0, 0, 0, 0, false);
        labRegionSaturation->setValue(r.saturation);
        labRegionSlope->setValue(r.slope);
        labRegionOffset->setValue(r.offset);
        labRegionPower->setValue(r.power);
        labRegionHueMask->setCurve(r.hueMask);
        labRegionChromaticityMask->setCurve(r.chromaticityMask);
        labRegionLightnessMask->setCurve(r.lightnessMask);
        labRegionMaskBlur->setValue(r.maskBlur);
        labRegionChannel->set_active(r.channel+1);
    }
    labRegionList->set_text(idx, 1, Glib::ustring::compose("a=%1 b=%2 S=%3\ns=%4 o=%5 p=%6", round_ab(r.a), round_ab(r.b), r.saturation, r.slope, r.offset, r.power));
    const char *ch = "";
    switch (r.channel) {
    case rtengine::procparams::ColorToningParams::LabCorrectionRegion::CHAN_R:
        ch = "\n[Red]"; break;
    case rtengine::procparams::ColorToningParams::LabCorrectionRegion::CHAN_G:
        ch = "\n[Green]"; break;
    case rtengine::procparams::ColorToningParams::LabCorrectionRegion::CHAN_B:
        ch = "\n[Blue]"; break;
    default:
        ch = "";
    }
    labRegionList->set_text(
        idx, 2, Glib::ustring::compose(
            "%1%2%3%4%5",
            hasMask(dflt.hueMask, r.hueMask) ? "H" : "",
            hasMask(dflt.chromaticityMask, r.chromaticityMask) ? "C" : "",
            hasMask(dflt.lightnessMask, r.lightnessMask) ? "L" : "",
            r.maskBlur ? Glib::ustring::compose(" b=%1", r.maskBlur) : "", ch));
    Gtk::TreePath pth;
    pth.push_back(idx);
    labRegionList->get_selection()->select(pth);
    if (disable) {
        enableListener();
    }
}


void ColorToning::labRegionChannelChanged()
{
    if (listener) {
        listener->panelChanged(EvLabRegionChannel, labRegionChannel->get_active_text());
    }
}


void ColorToning::setEditProvider(EditDataProvider *provider)
{
    labRegionHueMask->setEditProvider(provider);
    labRegionChromaticityMask->setEditProvider(provider);
    labRegionLightnessMask->setEditProvider(provider);
}


float ColorToning::blendPipetteValues(CurveEditor *ce, float chan1, float chan2, float chan3)
{
    if (ce == labRegionChromaticityMask && chan1 > 0.f) {
        return lin2log(chan1, 10.f);
    } else if (ce == labRegionHueMask && chan1 > 0.f) {
        float x = chan1 + 1.f/6.f;
        if (x > 1.f) {
            x -= 1.f;
        }
        return lin2log(x, 3.f);
    }
    return CurveListener::blendPipetteValues(ce, chan1, chan2, chan3);
}
