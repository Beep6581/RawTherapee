/*
 *  This file is part of RawTherapee.
 */
#include "colortoning.h"
#include "mycurve.h"

using namespace rtengine;
using namespace rtengine::procparams;

ColorToning::ColorToning () : FoldableToolPanel(this, "colortoning", M("TP_COLORTONING_LABEL"), false, true)
{	
	nextbw=0;
	CurveListener::setMulti(true);

	//---------------method
	
	method = Gtk::manage (new MyComboBoxText ());
	method->append_text (M("TP_COLORTONING_LAB"));
	method->append_text (M("TP_COLORTONING_RGBSLIDERS"));
	method->append_text (M("TP_COLORTONING_RGBCURVES"));
	method->append_text (M("TP_COLORTONING_SPLITCOCO"));
	method->append_text (M("TP_COLORTONING_SPLITLR"));
	method->set_active (0);
	method->set_tooltip_text (M("TP_COLORTONING_METHOD_TOOLTIP"));

	ctbox = Gtk::manage (new Gtk::HBox ());
	Gtk::Label* lab = Gtk::manage (new Gtk::Label (M("TP_COLORTONING_METHOD")));
	ctbox->pack_start (*lab, Gtk::PACK_SHRINK, 4);
	ctbox->pack_start (*method);
	pack_start (*ctbox);

	methodconn = method->signal_changed().connect ( sigc::mem_fun(*this, &ColorToning::methodChanged) );
	
	//----------- Color curve ------------------------------

	colorSep = Gtk::manage (new  Gtk::HSeparator());
	pack_start (*colorSep);

	colLabel = Gtk::manage (new Gtk::Label (M("TP_COLORTONING_LABCOL")));
	colLabel->set_tooltip_text (M("TP_COLORTONING_LABCOL_TOOLTIP"));

	interLabel = Gtk::manage (new Gtk::Label (M("TP_COLORTONING_LABINT")));
	interLabel->set_tooltip_text (M("TP_COLORTONING_LABINT_TOOLTIP"));
	pack_start (*colLabel, Gtk::PACK_SHRINK, 4);
	pack_start (*interLabel, Gtk::PACK_SHRINK, 4);

	colorCurveEditorG = new CurveEditorGroup (options.lastColorToningCurvesDir, M("TP_COLORTONING_COLOR"));
	colorCurveEditorG->setCurveListener (this);

	colorShape = static_cast<FlatCurveEditor*>(colorCurveEditorG->addCurve(CT_Flat, "", NULL, false));
	colorShape->setCurveColorProvider(this, 1);
	std::vector<GradientMilestone> milestones;
	// whole hue range
	for (int i=0; i<7; i++) {
		float R, G, B;
		float x = float(i)*(1.0f/6.0);
		Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
		milestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
	}
	colorShape->setLeftBarBgGradient(milestones);

	// luminance gradient
	milestones.clear();
	milestones.push_back( GradientMilestone(0., 0., 0., 0.) );
	milestones.push_back( GradientMilestone(1., 1., 1., 1.) );
	colorShape->setBottomBarBgGradient(milestones);
	std::vector<double> defaultCurve;
	rtengine::ColorToningParams::getDefaultColorCurve(defaultCurve);
	colorShape->setResetCurve(FCT_MinMaxCPoints, defaultCurve);

	// This will add the reset button at the end of the curveType buttons
	colorCurveEditorG->curveListComplete();
	colorCurveEditorG->show();

	pack_start( *colorCurveEditorG, Gtk::PACK_SHRINK, 2);

	//----------------------red green  blue yellow colours

	twocolor = Gtk::manage (new MyComboBoxText ());
	twocolor->append_text (M("TP_COLORTONING_TWOSTD"));
	twocolor->append_text (M("TP_COLORTONING_TWOALL"));
	twocolor->append_text (M("TP_COLORTONING_TWOBY"));
	twocolor->append_text (M("TP_COLORTONING_TWO2"));
	twocolor->set_tooltip_text (M("TP_COLORTONING_TWOCOLOR_TOOLTIP"));
	twocolor->set_active (0);

	twocconn = twocolor->signal_changed().connect( sigc::mem_fun(*this, &ColorToning::twoColorChangedByGui) );

	pack_start (*twocolor, Gtk::PACK_SHRINK, 4);

	//----------- Opacity curve ------------------------------

	opacityCurveEditorG = new CurveEditorGroup (options.lastColorToningCurvesDir, M("TP_COLORTONING_OPACITY"));
	opacityCurveEditorG->setCurveListener (this);

	rtengine::ColorToningParams::getDefaultOpacityCurve(defaultCurve);
	opacityShape = static_cast<FlatCurveEditor*>(opacityCurveEditorG->addCurve(CT_Flat, "", NULL, false));
	opacityShape->setIdentityValue(0.);
	opacityShape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
	opacityShape->setBottomBarBgGradient(milestones);

	// This will add the reset button at the end of the curveType buttons
	opacityCurveEditorG->curveListComplete();
	opacityCurveEditorG->show();

	pack_start( *opacityCurveEditorG, Gtk::PACK_SHRINK, 2);

	//---------Chroma curve 1 --------------------
	iby   = Gtk::manage (new RTImage ("Chanmixer-BY.png"));
	irg   = Gtk::manage (new RTImage ("Chanmixer-RG.png"));
	
	clCurveEditorG = new CurveEditorGroup (options.lastColorToningCurvesDir, M("TP_COLORTONING_CHROMAC"));
	clCurveEditorG->setCurveListener (this);

	rtengine::ColorToningParams::getDefaultCLCurve(defaultCurve);
	clshape = static_cast<DiagonalCurveEditor*>(clCurveEditorG->addCurve(CT_Diagonal, M("TP_COLORTONING_AB"),irg));
	clshape->setResetCurve(DiagonalCurveType(defaultCurve.at(0)), defaultCurve);
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

	rtengine::ColorToningParams::getDefaultCL2Curve(defaultCurve);
	cl2shape = static_cast<DiagonalCurveEditor*>(cl2CurveEditorG->addCurve(CT_Diagonal, M("TP_COLORTONING_BY"),iby));
	cl2shape->setResetCurve(DiagonalCurveType(defaultCurve.at(0)), defaultCurve);
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
	RTImage *resetImgc = Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
	neutralCurves->set_image(*resetImgc);
	neutralCurves->set_tooltip_text (M("TP_COLORTONING_NEUTRALCUR_TIP"));
	neutralcurvesconn = neutralCurves->signal_pressed().connect( sigc::mem_fun(*this, &ColorToning::neutralCurves_pressed) );
	neutralCurves->show();

	pack_start (*neutralCurves);
	*/

	//----------- Sliders + balance ------------------------------

	hlColSat = Gtk::manage (new ThresholdAdjuster (M("TP_COLORTONING_HIGHLIGHT"), 0., 100., 60., M("TP_COLORTONING_STRENGTH"), 1., 0., 360., 80., M("TP_COLORTONING_HUE"), 1., NULL, false));
	hlColSat->setAdjusterListener (this);
	hlColSat->setBgColorProvider(this, 2);
	hlColSat->setUpdatePolicy(RTUP_DYNAMIC);

	pack_start( *hlColSat, Gtk::PACK_SHRINK, 0);

	shadowsColSat = Gtk::manage (new ThresholdAdjuster (M("TP_COLORTONING_SHADOWS"), 0., 100., 80., M("TP_COLORTONING_STRENGTH"), 1., 0., 360., 208., M("TP_COLORTONING_HUE"), 1., NULL, false));
	shadowsColSat->setAdjusterListener (this);
	shadowsColSat->setBgColorProvider(this, 3);
	shadowsColSat->setUpdatePolicy(RTUP_DYNAMIC);

	pack_start( *shadowsColSat, Gtk::PACK_SHRINK, 0);


	balance = Gtk::manage( new Adjuster(M("TP_COLORTONING_BALANCE"), -100., 100., 1., 0.) );
	balance->setAdjusterListener(this);

	pack_start (*balance, Gtk::PACK_SHRINK, 2);

	//----------- Saturation and strength ------------------------------

//	satLimiterSep = Gtk::manage (new Gtk::HSeparator());
	
	
//	pack_start (*satLimiterSep, Gtk::PACK_SHRINK);
	
//	Gtk::Frame *p1Frame;
	// Vertical box container for the content of the Process 1 frame
	Gtk::VBox *p1VBox;
	p1Frame = Gtk::manage (new Gtk::Frame(M("TP_COLORTONING_SA")) );
	p1Frame->set_border_width(0);
	p1Frame->set_label_align(0.025, 0.5);
	
	p1VBox = Gtk::manage ( new Gtk::VBox());
	p1VBox->set_border_width(4);
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

	
	
	//	--------------------Sliders BW Colortoning -------------------

	chanMixerBox = Gtk::manage (new Gtk::VBox());
	Gtk::VBox *chanMixerHLBox = Gtk::manage (new Gtk::VBox());
	Gtk::VBox *chanMixerMidBox = Gtk::manage (new Gtk::VBox());
	Gtk::VBox *chanMixerShadowsBox = Gtk::manage (new Gtk::VBox());

	Gtk::Image* iblueR   = Gtk::manage (new RTImage ("ajd-wb-temp1.png"));
	Gtk::Image* iyelL    = Gtk::manage (new RTImage ("ajd-wb-temp2.png"));
	Gtk::Image* imagL    = Gtk::manage (new RTImage ("ajd-wb-green1.png"));
	Gtk::Image* igreenR  = Gtk::manage (new RTImage ("ajd-wb-green2.png"));
	Gtk::Image* icyanL   = Gtk::manage (new RTImage ("ajd-wb-bluered1.png"));
	Gtk::Image* iredR    = Gtk::manage (new RTImage ("ajd-wb-bluered2.png"));

	Gtk::Image* iblueRm  = Gtk::manage (new RTImage ("ajd-wb-temp1.png"));
	Gtk::Image* iyelLm   = Gtk::manage (new RTImage ("ajd-wb-temp2.png"));
	Gtk::Image* imagLm   = Gtk::manage (new RTImage ("ajd-wb-green1.png"));
	Gtk::Image* igreenRm = Gtk::manage (new RTImage ("ajd-wb-green2.png"));
	Gtk::Image* icyanLm  = Gtk::manage (new RTImage ("ajd-wb-bluered1.png"));
	Gtk::Image* iredRm   = Gtk::manage (new RTImage ("ajd-wb-bluered2.png"));

	Gtk::Image* iblueRh  = Gtk::manage (new RTImage ("ajd-wb-temp1.png"));
	Gtk::Image* iyelLh   = Gtk::manage (new RTImage ("ajd-wb-temp2.png"));
	Gtk::Image* imagLh   = Gtk::manage (new RTImage ("ajd-wb-green1.png"));
	Gtk::Image* igreenRh = Gtk::manage (new RTImage ("ajd-wb-green2.png"));
	Gtk::Image* icyanLh  = Gtk::manage (new RTImage ("ajd-wb-bluered1.png"));
	Gtk::Image* iredRh   = Gtk::manage (new RTImage ("ajd-wb-bluered2.png"));

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

	chanMixerHLFrame->add(*chanMixerHLBox);
	chanMixerMidFrame->add(*chanMixerMidBox);
	chanMixerShadowsFrame->add(*chanMixerShadowsBox);

	chanMixerBox->pack_start(*chanMixerHLFrame, Gtk::PACK_SHRINK);
	chanMixerBox->pack_start(*chanMixerMidFrame, Gtk::PACK_SHRINK);
	chanMixerBox->pack_start(*chanMixerShadowsFrame, Gtk::PACK_SHRINK);

	pack_start(*chanMixerBox, Gtk::PACK_SHRINK);
	pack_start( *strength, Gtk::PACK_SHRINK, 2); //I have moved after Chanmixer 

	//--------------------- Reset sliders  ---------------------------
	neutrHBox = Gtk::manage (new Gtk::HBox ());
	neutrHBox->set_border_width (2);

	neutral = Gtk::manage (new Gtk::Button (M("TP_COLORTONING_NEUTRAL")));
	RTImage *resetImg = Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
	neutral->set_image(*resetImg);
	neutral->set_tooltip_text (M("TP_COLORTONING_NEUTRAL_TIP"));
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

	show_all();

	disableListener();
	methodChanged();
	enableListener();
}

ColorToning::~ColorToning() {
	delete colorCurveEditorG;
	delete opacityCurveEditorG;
	delete clCurveEditorG;
	delete cl2CurveEditorG;
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
		listener->panelChanged (EvColorToningNeutralcur, M("ADJUSTER_RESET_TO_DEFAULT"));
}
*/

// Will only reset the chanel mixer
void ColorToning::neutral_pressed () {
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
	if (listener && getEnabled())
		listener->panelChanged (EvColorToningNeutral, M("ADJUSTER_RESET_TO_DEFAULT"));
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
	
	if (pedited && !pedited->colorToning.method)
		method->set_active (5);
	else if (pp->colorToning.method=="Lab")
		method->set_active (0);
	else if (pp->colorToning.method=="RGBSliders")
		method->set_active (1);
	else if (pp->colorToning.method=="RGBCurves")
		method->set_active (2);
	else if (pp->colorToning.method=="Splitco")
		method->set_active (3);
	else if (pp->colorToning.method=="Splitlr")
		method->set_active (4);
	methodChanged();
	methodconn.block(false);


	if (pedited && !pedited->colorToning.twocolor)
		twocolor->set_active (4);
	else if (pp->colorToning.twocolor=="Std")
		twocolor->set_active (0);
	else if (pp->colorToning.twocolor=="All")
		twocolor->set_active (1);
	else if (pp->colorToning.twocolor=="Separ")
		twocolor->set_active (2);			
	else if (pp->colorToning.twocolor=="Two")
		twocolor->set_active (3);

	twocolorChanged(true);

	twocconn.block(false);

	enableListener ();
}

void ColorToning::write (ProcParams* pp, ParamsEdited* pedited) {
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
		pedited->colorToning.method     = method->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->colorToning.twocolor   = twocolor->get_active_text()!=M("GENERAL_UNCHANGED");

		pedited->colorToning.enabled       = !get_inconsistent();
		pedited->colorToning.autosat       = !autosat->get_inconsistent();
		pedited->colorToning.colorCurve    = !colorShape->isUnChanged ();
		pedited->colorToning.opacityCurve  = !opacityShape->isUnChanged ();
		pedited->colorToning.clcurve       = !clshape->isUnChanged ();
		pedited->colorToning.cl2curve      = !cl2shape->isUnChanged ();
		pedited->colorToning.lumamode      = !lumamode->get_inconsistent();

		pedited->colorToning.hlColSat      = hlColSat->getEditedState ();
		pedited->colorToning.shadowsColSat = shadowsColSat->getEditedState ();
	}

	if (method->get_active_row_number()==0)
		pp->colorToning.method = "Lab";
	else if (method->get_active_row_number()==1)
		pp->colorToning.method = "RGBSliders";
	else if (method->get_active_row_number()==2)
		pp->colorToning.method = "RGBCurves";
	else if (method->get_active_row_number()==3)
		pp->colorToning.method = "Splitco";
	else if (method->get_active_row_number()==4)
		pp->colorToning.method = "Splitlr";

	if (twocolor->get_active_row_number()==0)
		pp->colorToning.twocolor = "Std";
	else if (twocolor->get_active_row_number()==1)
		pp->colorToning.twocolor = "All";
	else if (twocolor->get_active_row_number()==2)
		pp->colorToning.twocolor = "Separ";
	else if (twocolor->get_active_row_number()==3)
		pp->colorToning.twocolor = "Two";
}

void ColorToning::lumamodeChanged () {

	if (batchMode) {
		if (lumamode->get_inconsistent()) {
			lumamode->set_inconsistent (false);
			lumamodeConn.block (true);
			lumamode->set_active (false);
			lumamodeConn.block (false);
		}
		else if (lastLumamode)
			lumamode->set_inconsistent (true);

		lastLumamode = lumamode->get_active ();
	}

	if (listener && getEnabled()) {
		if (lumamode->get_active ())
			listener->panelChanged (EvColorToningLumamode, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvColorToningLumamode, M("GENERAL_DISABLED"));
	}
}

void ColorToning::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

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
	}
	else {
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
	}
}

void ColorToning::setAdjusterBehavior (bool splitAdd, bool satThresholdAdd, bool satOpacityAdd, bool strprotectAdd, bool balanceAdd) {
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

void ColorToning::adjusterChanged (ThresholdAdjuster* a, double newBottom, double newTop) {
	if (listener && getEnabled())
		listener->panelChanged (a==hlColSat ? EvColorToningHighights : EvColorToningShadows,
								Glib::ustring::compose(Glib::ustring(M("TP_COLORTONING_HUE")+": %1"+"\n"+M("TP_COLORTONING_STRENGTH")+": %2"), int(newTop), int(newBottom)));
}

int CTChanged_UI (void* data) {
	GThreadLock lock;
	(static_cast<ColorToning*>(data))->CTComp_ ();
	return 0;
}


void ColorToning::autoColorTonChanged(int bwct, int satthres, int satprot){
	nextbw = bwct;
	nextsatth=satthres;
	nextsatpr=satprot;
	g_idle_add (CTChanged_UI, this);
}

bool ColorToning::CTComp_ () {

	disableListener ();
	saturatedOpacity->setValue (nextsatpr);
	satProtectionThreshold->setValue (nextsatth);
/*	if(nextbw==1) {
		saturatedOpacity->show();
		satProtectionThreshold->show();
		autosat->show();
	}
	else {
		saturatedOpacity->hide();
		satProtectionThreshold->hide();
		autosat->hide();
	}
*/	
	enableListener ();

	return false;
}

void ColorToning::adjusterChanged (Adjuster* a, double newval) {

	if (!listener || !getEnabled())
		return;

	if (a==redlow)
		listener->panelChanged (EvColorToningredlow, redlow->getTextValue());
	else if (a==greenlow)
		listener->panelChanged (EvColorToninggreenlow, greenlow->getTextValue());
	else if (a==bluelow)
		listener->panelChanged (EvColorToningbluelow, bluelow->getTextValue());
	else if (a==redmed)
		listener->panelChanged (EvColorToningredmed, redmed->getTextValue());
	else if (a==greenmed)
		listener->panelChanged (EvColorToninggreenmed, greenmed->getTextValue());
	else if (a==bluemed)
		listener->panelChanged (EvColorToningbluemed, bluemed->getTextValue());
	else if (a==redhigh)
		listener->panelChanged (EvColorToningredhigh, redhigh->getTextValue());
	else if (a==greenhigh)
		listener->panelChanged (EvColorToninggreenhigh, greenhigh->getTextValue());
	else if (a==bluehigh)
		listener->panelChanged (EvColorToningbluehigh, bluehigh->getTextValue());
	else if (a==balance)
		listener->panelChanged (EvColorToningbalance, balance->getTextValue());
	else if (a==satProtectionThreshold)
		listener->panelChanged (EvColorToningSatThreshold, a->getTextValue());
	else if (a==saturatedOpacity)
		listener->panelChanged (EvColorToningSatProtection, a->getTextValue());
	else if (a==strength)
		listener->panelChanged (EvColorToningStrength, a->getTextValue());
}

//Two Color changed
void ColorToning::twocolorChanged (bool changedbymethod) {
	if (!batchMode) {
		if(method->get_active_row_number()==0) {				// Lab
			if(twocolor->get_active_row_number()==0) {
				colorCurveEditorG->show();   // visible
				opacityCurveEditorG->show(); // visible
				clCurveEditorG->hide();
				cl2CurveEditorG->hide();
			}
			else if(twocolor->get_active_row_number()==1  || twocolor->get_active_row_number()==3) {
				colorCurveEditorG->show();   // visible
				opacityCurveEditorG->hide();
				clCurveEditorG->show();      // visible
				cl2CurveEditorG->hide();
				irg->hide();
				
			}
			else if(twocolor->get_active_row_number()==2) {
				colorCurveEditorG->show(); // visible
				opacityCurveEditorG->hide();
				clCurveEditorG->show();      // visible
				cl2CurveEditorG->show();     // visible
				irg->show();
			}
		}
		else if(method->get_active_row_number()==1) {			// RGB Sliders
			colorCurveEditorG->hide();
			opacityCurveEditorG->hide();
			clCurveEditorG->hide();
			cl2CurveEditorG->hide();
		}
		else if(method->get_active_row_number()==2) {			// RGB Curves
			colorCurveEditorG->show();       // visible
			opacityCurveEditorG->show();     // visible
			clCurveEditorG->hide();
			cl2CurveEditorG->hide();
		}
		else if(method->get_active_row_number()==3) {			// Split LR
			colorCurveEditorG->hide();
			opacityCurveEditorG->hide();
			clCurveEditorG->hide();
			cl2CurveEditorG->hide();
		}
		else if(method->get_active_row_number()==4) {			// Split color
			colorCurveEditorG->hide();
			opacityCurveEditorG->hide();
			clCurveEditorG->hide();
			cl2CurveEditorG->hide();
		}
	}

	if (listener && getEnabled() && !changedbymethod)
		listener->panelChanged (EvColorToningTwocolor, twocolor->get_active_text ());
}

void ColorToning::twoColorChangedByGui() {
	twocolorChanged(false);
}

void ColorToning::methodChanged () {

	if (!batchMode) {
		if (method->get_active_row_number()==0) {  // Lab
			colorSep->show();
			colLabel->hide();
			interLabel->hide();
			colorCurveEditorG->show();
			twocolor->show();
			opacityCurveEditorG->hide();
			clCurveEditorG->hide();
			cl2CurveEditorG->hide();
			//neutralCurves->show();
			hlColSat->hide();
			shadowsColSat->hide();
			balance->hide();
//			satLimiterSep->show();
			if(autosat->get_active()) {
			saturatedOpacity->set_sensitive(false);
			satProtectionThreshold->set_sensitive(false);
			satProtectionThreshold->show();
			saturatedOpacity->show();			
			}
			else {
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
		}
		else if (method->get_active_row_number()==1) {  // RGB Sliders
			colorSep->hide();
			colLabel->hide();
			interLabel->hide();
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
//			satLimiterSep->show();
			autosat->show();
			p1Frame->show();
			if(autosat->get_active()) {
			saturatedOpacity->set_sensitive(false);
			satProtectionThreshold->set_sensitive(false);
			satProtectionThreshold->show();
			saturatedOpacity->show();			
			}
			else {
			satProtectionThreshold->show();
			saturatedOpacity->show();
			saturatedOpacity->set_sensitive(true);
			satProtectionThreshold->set_sensitive(true);
			}
			
			strength->hide();
			chanMixerBox->hide();
			neutrHBox->hide();
			lumamode->hide();

		}
		else if (method->get_active_row_number()==2) {  // RGB Curves
			colorSep->hide();
			colLabel->hide();
			interLabel->hide();
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
//			satLimiterSep->show();
			p1Frame->show();

			autosat->show();
			if(autosat->get_active()) {
			saturatedOpacity->set_sensitive(false);
			satProtectionThreshold->set_sensitive(false);
			satProtectionThreshold->show();
			saturatedOpacity->show();			
			}
			else {
			satProtectionThreshold->show();
			saturatedOpacity->show();
			saturatedOpacity->set_sensitive(true);
			satProtectionThreshold->set_sensitive(true);
			}
			strength->hide();
			chanMixerBox->hide();
			neutrHBox->hide();
			lumamode->hide();
		}
		else if (method->get_active_row_number()==3) {  // Split LR
			colorSep->hide();
			colLabel->hide();
			interLabel->hide();
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
			
//			satLimiterSep->hide();
			autosat->hide();
			satProtectionThreshold->hide();
			saturatedOpacity->hide();
			strength->show();
			chanMixerBox->show();
			neutrHBox->show();
			lumamode->show();
		}
		else if (method->get_active_row_number()==4) {  // Split Color
			colorSep->show();
			colLabel->hide();
			interLabel->hide();
			colorCurveEditorG->hide();
			twocolor->hide();
			opacityCurveEditorG->hide();
			clCurveEditorG->hide();
			cl2CurveEditorG->hide();
			//neutralCurves->hide();
			hlColSat->show();
			shadowsColSat->show();
			balance->show();
//			satLimiterSep->hide();
			p1Frame->hide();

			autosat->hide();
			satProtectionThreshold->hide();
			saturatedOpacity->hide();
			strength->show();

			chanMixerBox->hide();
			neutrHBox->hide();
			lumamode->show();
		}
	}

	if (listener && getEnabled())
		listener->panelChanged (EvColorToningMethod, method->get_active_text ());
}


void ColorToning::autoOpenCurve  () {
	colorShape->openIfNonlinear();
	opacityShape->openIfNonlinear();
	clshape->openIfNonlinear();
	cl2shape->openIfNonlinear();
}

void ColorToning::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller) {

	float R, G, B;
	if (callerId == 1) {         // ch - main curve
		Color::hsv2rgb01(float(valY), 1.0f, 0.5f, R, G, B);
	}
	else if (callerId == 2) {    // Slider 1 background
		if (valY > 0.5)
			// the hue range
			Color::hsv2rgb01(float(valX), 1.0f, 0.5f, R, G, B);
		else {
			// the strength applied to the current hue
			double strength, hue;
			float r_, g_, b_;
			hlColSat->getValue(strength, hue);
			Color::hsv2rgb01(valY*2.f, 1.f, 1.f, r_, g_, b_);
			Color::hsv2rgb01(hue/360.f, 1.f, 1.f, R, G, B);
			R = r_+(R-r_)*valX;
			G = g_+(G-g_)*valX;
			B = b_+(B-b_)*valX;
		}
	}
	else if (callerId == 3) {    // Slider 2 background
		if (valY > 0.5)
			// the hue range
			Color::hsv2rgb01(float(valX), 1.0f, 0.5f, R, G, B);
		else {
			// the strength applied to the current hue
			double strength, hue;
			float r_, g_, b_;
			shadowsColSat->getValue(strength, hue);
			Color::hsv2rgb01(valY*2.f, 1.f, 1.f, r_, g_, b_);
			Color::hsv2rgb01(hue/360.f, 1.f, 1.f, R, G, B);
			R = r_+(R-r_)*valX;
			G = g_+(G-g_)*valX;
			B = b_+(B-b_)*valX;
		}
	}
	caller->ccRed = double(R);
	caller->ccGreen = double(G);
	caller->ccBlue = double(B);
}

void ColorToning::curveChanged (CurveEditor* ce) {

	if (listener && getEnabled()) {
		if (ce == colorShape)
			listener->panelChanged (EvColorToningColor, M("HISTORY_CUSTOMCURVE"));
		else if (ce == opacityShape)
			listener->panelChanged (EvColorToningOpacity, M("HISTORY_CUSTOMCURVE"));
		else if (ce == clshape)
			listener->panelChanged (EvColorToningCLCurve, M("HISTORY_CUSTOMCURVE"));
		else if (ce == cl2shape)
			listener->panelChanged (EvColorToningLLCurve, M("HISTORY_CUSTOMCURVE"));
	}
}

void ColorToning::enabledChanged () {

	if (listener) {
		if (get_inconsistent())
			listener->panelChanged (EvColorToningEnabled, M("GENERAL_UNCHANGED"));
		else if (getEnabled())
			listener->panelChanged (EvColorToningEnabled, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvColorToningEnabled, M("GENERAL_DISABLED"));
	}
}

void ColorToning::autosatChanged () {

	if (batchMode) {
		if (autosat->get_inconsistent()) {
			autosat->set_inconsistent (false);
			autosatConn.block (true);
			autosat->set_active (false);
			autosatConn.block (false);
		}
		else if (lastautosat)
			autosat->set_inconsistent (true);

		lastautosat = autosat->get_active ();
	}
	if (listener) {
		if (autosat->get_active()) {
			if (getEnabled())
				listener->panelChanged (EvColorToningautosat, M("GENERAL_ENABLED"));
			saturatedOpacity->set_sensitive(false);
			satProtectionThreshold->set_sensitive(false);
		}
		else {
			if (getEnabled())
				listener->panelChanged (EvColorToningautosat, M("GENERAL_DISABLED"));
			saturatedOpacity->set_sensitive(true);
			satProtectionThreshold->set_sensitive(true);
		}

	}
}

void ColorToning::trimValues (rtengine::procparams::ProcParams* pp) {

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

void ColorToning::setBatchMode (bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	method->append_text (M("GENERAL_UNCHANGED"));
	twocolor->append_text (M("GENERAL_UNCHANGED"));
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
