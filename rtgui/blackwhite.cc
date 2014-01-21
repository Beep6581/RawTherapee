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
 *  GNU General Public License for more details
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "blackwhite.h"
#include "rtimage.h"
#include "../rtengine/color.h"
#include <iomanip>
#include <cmath>
#include "guiutils.h"

using namespace rtengine;
using namespace rtengine::procparams;


BlackWhite::BlackWhite (): Gtk::VBox(), FoldableToolPanel(this) {
	CurveListener::setMulti(true);
	set_border_width(4);
	set_spacing(4);

	nextredbw = 0.3333;
	nextgreenbw = 0.3333;
	nextbluebw = 0.3333;

	//----------- Enables checkbox ------------------------------

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);
	
	pack_start(*enabled, Gtk::PACK_SHRINK, 0);
	enabled->show ();
	enaconn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &BlackWhite::enabled_toggled) );

	pack_start (*Gtk::manage (new  Gtk::HSeparator()));

	//----------- Method combobox ------------------------------

	Gtk::HBox* metHBox = Gtk::manage (new Gtk::HBox ());
	metHBox->set_border_width (0);
	metHBox->set_spacing (2);
	Gtk::Label* metLabel = Gtk::manage (new Gtk::Label (M("TP_BWMIX_MET")+":"));
	metHBox->pack_start (*metLabel, Gtk::PACK_SHRINK);
	method = Gtk::manage (new MyComboBoxText ());
	method->append_text (M("TP_BWMIX_MET_DESAT"));
	method->append_text (M("TP_BWMIX_MET_LUMEQUAL"));
	method->append_text (M("TP_BWMIX_MET_CHANMIX"));

	method->set_active (0);
	metHBox->pack_start (*method);
	pack_start (*metHBox);
	methodconn = method->signal_changed().connect ( sigc::mem_fun(*this, &BlackWhite::methodChanged) );


	//----------- Luminance equalizer ------------------------------

	luminanceSep = Gtk::manage (new  Gtk::HSeparator());
	pack_start (*luminanceSep);

	std::vector<GradientMilestone> bottomMilestones;
	float R, G, B;
	// -0.1 rad < Hue < 1.6 rad
	for (int i=0; i<7; i++) {
		float x = float(i)*(1.0f/6.0);
		Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
		bottomMilestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
	}
	luminanceCEG = new CurveEditorGroup (options.lastBWCurvesDir, M("TP_BWMIX_CHANNEL"));
	luminanceCEG->setCurveListener (this);
	luminanceCurve = static_cast<FlatCurveEditor*>(luminanceCEG->addCurve(CT_Flat, M("TP_BWMIX_VAL")));
	luminanceCurve->setBottomBarBgGradient(bottomMilestones);
	luminanceCurve->setCurveColorProvider(this, 3);
	luminanceCurve->setTooltip(M("TP_BWMIX_CURVEEDITOR_LH_TOOLTIP"));
	
	luminanceCEG->curveListComplete();
	pack_start (*luminanceCEG, Gtk::PACK_SHRINK, 4);

	//----------- Auto and Reset buttons ------------------------------

	mixerFrame = Gtk::manage (new Gtk::Frame (M("TP_BWMIX_MET_CHANMIX")));
	pack_start (*mixerFrame, Gtk::PACK_SHRINK, 0);

	mixerVBox = Gtk::manage (new Gtk::VBox ());
	mixerVBox->set_border_width(4);
	mixerVBox->set_spacing(4);

	autoHBox = Gtk::manage (new Gtk::HBox ());
	autoHBox->set_border_width (2);

	autoch = Gtk::manage (new Gtk::ToggleButton (M("TP_BWMIX_AUTOCH")));
	autoch->set_tooltip_markup (M("TP_BWMIX_AUTOCH_TIP"));
	autoconn = autoch->signal_toggled().connect( sigc::mem_fun(*this, &BlackWhite::autoch_toggled) );

	neutral = Gtk::manage (new Gtk::Button (M("TP_BWMIX_NEUTRAL")));
	RTImage *resetImg = Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
	neutral->set_image(*resetImg);
	neutral->set_tooltip_text (M("TP_BWMIX_NEUTRAL_TIP"));
	neutralconn = neutral->signal_pressed().connect( sigc::mem_fun(*this, &BlackWhite::neutral_pressed) );
	neutral->show();

	autoHBox->pack_start (*autoch);
	autoHBox->pack_end (*neutral);
	autoHBox->pack_end (*Gtk::manage (new Gtk::Label (" "))); //spacer
	mixerVBox->pack_start (*autoHBox);

	//----------- Presets combobox ------------------------------

	mixerVBox->pack_start (*Gtk::manage (new  Gtk::HSeparator()));

	settingHBox = Gtk::manage (new Gtk::HBox ());
	settingHBox->set_border_width (0);
	settingHBox->set_spacing (2);
	settingHBox->set_tooltip_markup (M("TP_BWMIX_SETTING_TOOLTIP"));
	Gtk::Label *settingLabel = Gtk::manage (new Gtk::Label (M("TP_BWMIX_SETTING")+":"));

	settingHBox->pack_start (*settingLabel, Gtk::PACK_SHRINK);
	setting = Gtk::manage (new MyComboBoxText ());
	setting->append_text (M("TP_BWMIX_SET_NORMCONTAST"));
	setting->append_text (M("TP_BWMIX_SET_HIGHCONTAST"));
	setting->append_text (M("TP_BWMIX_SET_LUMINANCE"));
	setting->append_text (M("TP_BWMIX_SET_LANDSCAPE"));
	setting->append_text (M("TP_BWMIX_SET_PORTRAIT"));
	setting->append_text (M("TP_BWMIX_SET_LOWSENSIT"));
	setting->append_text (M("TP_BWMIX_SET_HIGHSENSIT"));
	setting->append_text (M("TP_BWMIX_SET_PANCHRO"));
	setting->append_text (M("TP_BWMIX_SET_HYPERPANCHRO"));
	setting->append_text (M("TP_BWMIX_SET_ORTHOCHRO"));
	setting->append_text (M("TP_BWMIX_SET_RGBABS"));
	setting->append_text (M("TP_BWMIX_SET_RGBREL"));
	setting->append_text (M("TP_BWMIX_SET_ROYGCBPMABS"));
	setting->append_text (M("TP_BWMIX_SET_ROYGCBPMREL"));
	setting->append_text (M("TP_BWMIX_SET_INFRARED"));

	setting->set_active (0);
	settingHBox->pack_start (*setting);
	mixerVBox->pack_start (*settingHBox);
	settingconn = setting->signal_changed().connect ( sigc::mem_fun(*this, &BlackWhite::settingChanged) );

	RGBLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
	RGBLabels->set_tooltip_text(M("TP_BWMIX_RGBLABEL_HINT"));
	mixerVBox->pack_start (*RGBLabels);

	//----------- Complementary Color checkbox ------------------------------

	enabledccSep= Gtk::manage (new  Gtk::HSeparator());
	mixerVBox->pack_start (*enabledccSep);

	enabledcc = Gtk::manage (new Gtk::CheckButton (M("TP_BWMIX_CC_ENABLED")));

	enabledcc->set_active (true);
	enabledcc->set_tooltip_markup (M("TP_BWMIX_CC_TOOLTIP"));

	mixerVBox->pack_start(*enabledcc, Gtk::PACK_SHRINK, 0);
	enabledcc->show ();
	enaccconn = enabledcc->signal_toggled().connect( sigc::mem_fun(*this, &BlackWhite::enabledcc_toggled) );

	//----------- Color Filters ------------------------------

	filterSep = Gtk::manage (new  Gtk::HSeparator());
	mixerVBox->pack_start (*filterSep);

	filterHBox = Gtk::manage (new Gtk::HBox ());
	filterHBox->set_border_width (0);
	filterHBox->set_spacing (2);
	filterHBox->set_tooltip_markup (M("TP_BWMIX_FILTER_TOOLTIP"));
	Gtk::Label *filterLabel = Gtk::manage (new Gtk::Label (M("TP_BWMIX_FILTER")+":"));
	filterHBox->pack_start (*filterLabel, Gtk::PACK_SHRINK);
	filter = Gtk::manage (new MyComboBoxText ());
	filter->append_text (M("TP_BWMIX_FILTER_NONE"));
	filter->append_text (M("TP_BWMIX_FILTER_RED"));
	filter->append_text (M("TP_BWMIX_FILTER_REDYELLOW"));
	filter->append_text (M("TP_BWMIX_FILTER_YELLOW"));
	filter->append_text (M("TP_BWMIX_FILTER_GREENYELLOW"));
	filter->append_text (M("TP_BWMIX_FILTER_GREEN"));
	filter->append_text (M("TP_BWMIX_FILTER_BLUEGREEN"));
	filter->append_text (M("TP_BWMIX_FILTER_BLUE"));
	filter->append_text (M("TP_BWMIX_FILTER_PURPLE"));

	filter->set_active (0);
	filterHBox->pack_start (*filter);
	mixerVBox->pack_start (*filterHBox);
	filterconn = filter->signal_changed().connect ( sigc::mem_fun(*this, &BlackWhite::filterChanged) );

	//----------- RGB / ROYGCBPM Mixer ------------------------------

	imgIcon[0] = Gtk::manage (new RTImage ("Chanmixer-R.png"));
	imgIcon[1] = Gtk::manage (new RTImage ("Chanmixer-O.png"));
	imgIcon[2] = Gtk::manage (new RTImage ("Chanmixer-Y.png"));
	imgIcon[3] = Gtk::manage (new RTImage ("Chanmixer-G.png"));
	imgIcon[4] = Gtk::manage (new RTImage ("Chanmixer-C.png"));
	imgIcon[5] = Gtk::manage (new RTImage ("Chanmixer-B.png"));
	imgIcon[6] = Gtk::manage (new RTImage ("Chanmixer-P.png"));
	imgIcon[7] = Gtk::manage (new RTImage ("Chanmixer-M.png"));

	imgIcon[8]  = Gtk::manage (new RTImage ("Chanmixer-Rgamma.png"));
	imgIcon[9]  = Gtk::manage (new RTImage ("Chanmixer-Ggamma.png"));
	imgIcon[10] = Gtk::manage (new RTImage ("Chanmixer-Bgamma.png"));

	mixerVBox->pack_start (*Gtk::manage (new  Gtk::HSeparator()));

	mixerRed= Gtk::manage(new Adjuster (/*M("TP_BWMIX_RED")*/"", -100, 200, 1, 33, imgIcon[0]));
	if (mixerRed->delay < 50) mixerRed->delay = 50;
	mixerRed->setAdjusterListener (this);
	mixerRed->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	mixerRed->show();
	mixerVBox->pack_start( *mixerRed, Gtk::PACK_SHRINK, 0);

	mixerGreen= Gtk::manage(new Adjuster (/*M("TP_BWMIX_GREEN")*/"", -100, 200, 1, 33, imgIcon[3]));
	if (mixerGreen->delay < 50) mixerGreen->delay = 50;
	mixerGreen->setAdjusterListener (this);
	mixerGreen->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	mixerGreen->show();
	mixerVBox->pack_start( *mixerGreen, Gtk::PACK_SHRINK, 0);

	mixerBlue= Gtk::manage(new Adjuster (/*M("TP_BWMIX_BLUE")*/"", -100, 200, 1, 33, imgIcon[5]));
	if (mixerBlue->delay < 50) mixerBlue->delay = 50;
	mixerBlue->setAdjusterListener (this);
	mixerBlue->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	mixerBlue->show();
	mixerVBox->pack_start( *mixerBlue, Gtk::PACK_SHRINK, 0);

	filterSep2 = Gtk::manage (new  Gtk::HSeparator());
	mixerVBox->pack_start (*filterSep2);

	mixerOrange= Gtk::manage(new Adjuster (/*M("TP_BWMIX_ORANGE")*/"", -100, 200, 1, 33, imgIcon[1]));
	if (mixerOrange->delay < 50) mixerOrange->delay = 50;
	mixerOrange->setAdjusterListener (this);
	mixerOrange->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	mixerOrange->show();
	mixerVBox->pack_start( *mixerOrange, Gtk::PACK_SHRINK, 0);

	mixerYellow= Gtk::manage(new Adjuster (/*M("TP_BWMIX_YELLOW")*/"", -100, 200, 1, 33, imgIcon[2]));
	if (mixerYellow->delay < 50) mixerYellow->delay = 50;
	mixerYellow->setAdjusterListener (this);
	mixerYellow->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	mixerYellow->show();
	mixerVBox->pack_start( *mixerYellow, Gtk::PACK_SHRINK, 0);

	mixerCyan= Gtk::manage(new Adjuster (/*M("TP_BWMIX_CYAN")*/"", -100, 200, 1, 33, imgIcon[4]));
	if (mixerCyan->delay < 50) mixerCyan->delay = 50;
	mixerCyan->setAdjusterListener (this);
	mixerCyan->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	mixerCyan->show();
	mixerVBox->pack_start( *mixerCyan, Gtk::PACK_SHRINK, 0);

	mixerPurple= Gtk::manage(new Adjuster (/*M("TP_BWMIX_PURPLE")*/"", -100, 200, 1, 33, imgIcon[6]));
	if (mixerPurple->delay < 50) mixerPurple->delay = 50;
	mixerPurple->setAdjusterListener (this);
	mixerPurple->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	mixerPurple->show();
	mixerVBox->pack_start( *mixerPurple, Gtk::PACK_SHRINK, 0);

	mixerMagenta= Gtk::manage(new Adjuster (/*M("TP_BWMIX_MAGENTA")*/"", -100, 200, 1, 33, imgIcon[7]));
	if (mixerMagenta->delay < 50) mixerMagenta->delay = 50;
	mixerMagenta->setAdjusterListener (this);
	mixerMagenta->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	mixerMagenta->show();
	mixerVBox->pack_start( *mixerMagenta, Gtk::PACK_SHRINK, 0);

	mixerFrame->add(*mixerVBox);

	//----------- Gamma sliders ------------------------------

	gammaFrame = Gtk::manage (new Gtk::Frame (M("TP_BWMIX_GAMMA")));
	pack_start (*gammaFrame, Gtk::PACK_SHRINK, 0);

	Gtk::VBox *gammaVBox = Gtk::manage (new Gtk::VBox());
	gammaVBox->set_spacing(4);
	gammaVBox->set_border_width(4);


	gammaRed= Gtk::manage(new Adjuster (/*M("TP_BWMIX_GAM_RED")*/"", -100, 100, 1, 0, imgIcon[8]));
	if (gammaRed->delay < 50) gammaRed->delay = 50;

	gammaRed->setAdjusterListener (this);
	gammaRed->set_tooltip_markup (M("TP_BWMIX_GAM_TOOLTIP"));
	gammaRed->show();
	gammaVBox->pack_start( *gammaRed, Gtk::PACK_SHRINK, 0);

	gammaGreen= Gtk::manage(new Adjuster (/*M("TP_BWMIX_GAM_GREEN")*/"", -100, 100, 1, 0, imgIcon[9]));
	if (gammaGreen->delay < 50) gammaGreen->delay = 50;
	gammaGreen->setAdjusterListener (this);
	gammaGreen->set_tooltip_markup (M("TP_BWMIX_GAM_TOOLTIP"));	
	gammaGreen->show();
	gammaVBox->pack_start( *gammaGreen, Gtk::PACK_SHRINK, 0);

	gammaBlue= Gtk::manage(new Adjuster (/*M("TP_BWMIX_GAM_BLUE")*/"", -100, 100, 1, 0, imgIcon[10]));
	if (gammaBlue->delay < 50) gammaBlue->delay = 50;
	gammaBlue->setAdjusterListener (this);
	gammaBlue->set_tooltip_markup (M("TP_BWMIX_GAM_TOOLTIP"));
	gammaBlue->show();
	gammaVBox->pack_start( *gammaBlue, Gtk::PACK_SHRINK, 0);

	gammaFrame->add(*gammaVBox);

	//----------- Curve 1 ------------------------------

	std::vector<GradientMilestone> bottomMilestonesbw;
	bottomMilestonesbw.push_back( GradientMilestone(0., 0., 0., 0.) );
	bottomMilestonesbw.push_back( GradientMilestone(1., 1., 1., 1.) );

	beforeCurveMode = Gtk::manage (new MyComboBoxText ());
	beforeCurveMode->append_text (M("TP_BWMIX_TCMODE_STANDARD"));
	beforeCurveMode->append_text (M("TP_BWMIX_TCMODE_WEIGHTEDSTD"));
	beforeCurveMode->append_text (M("TP_BWMIX_TCMODE_FILMLIKE"));
	beforeCurveMode->append_text (M("TP_BWMIX_TCMODE_SATANDVALBLENDING"));
	beforeCurveMode->set_active (0);

	beforeCurveCEG = new CurveEditorGroup (options.lastBWCurvesDir, M("TP_BWMIX_CURVEEDITOR1"));
	beforeCurveCEG->setCurveListener (this);

	beforeCurve = static_cast<DiagonalCurveEditor*>(beforeCurveCEG->addCurve(CT_Diagonal, "", beforeCurveMode));
	beforeCurve->setBottomBarBgGradient(bottomMilestonesbw);
	beforeCurve->setLeftBarBgGradient(bottomMilestonesbw);
	beforeCurve->setTooltip(M("TP_BWMIX_CURVEEDITOR_BEFORE_TOOLTIP"));

	// This will add the reset button at the end of the curveType buttons
	beforeCurveCEG->curveListComplete();

	pack_start( *beforeCurveCEG, Gtk::PACK_SHRINK, 2);

	tcmodeconn = beforeCurveMode->signal_changed().connect( sigc::mem_fun(*this, &BlackWhite::curveMode1Changed), true );

	//----------- Curve 2 ------------------------------
/*
	afterCurveMode = Gtk::manage (new MyComboBoxText ());
	afterCurveMode->append_text (M("TP_BWMIX_TCMODE_STANDARD"));
	//  afterCurveMode->append_text (M("TP_BWMIX_TCMODE_WEIGHTEDSTD"));
	afterCurveMode->set_active (0);
*/
	afterCurveCEG = new CurveEditorGroup (options.lastBWCurvesDir, M("TP_BWMIX_CURVEEDITOR2"));
	afterCurveCEG->setCurveListener (this);

//	afterCurve = static_cast<DiagonalCurveEditor*>(afterCurveCEG->addCurve(CT_Diagonal, "", afterCurveMode));
	afterCurve = static_cast<DiagonalCurveEditor*>(afterCurveCEG->addCurve(CT_Diagonal, ""));
	afterCurve->setBottomBarBgGradient(bottomMilestonesbw);
	afterCurve->setLeftBarBgGradient(bottomMilestonesbw);
	afterCurve->setTooltip(M("TP_BWMIX_CURVEEDITOR_AFTER_TOOLTIP"));

	afterCurveCEG->curveListComplete();

	pack_start( *afterCurveCEG, Gtk::PACK_SHRINK, 2);

//	tcmodeconn2 = afterCurveMode->signal_changed().connect( sigc::mem_fun(*this, &BlackWhite::curveMode1Changed2), true );

	show_all();
}
BlackWhite::~BlackWhite () {
	delete luminanceCEG;
	delete beforeCurveCEG;
	delete afterCurveCEG;
}

int BWChangedUI (void* data) {
	GThreadLock lock;
	(static_cast<BlackWhite*>(data))->BWComputed_ ();
	return 0;
}

void BlackWhite::BWChanged  (double redbw, double greenbw, double bluebw){
	nextredbw = redbw;
	nextgreenbw = greenbw;
	nextbluebw = bluebw;
	g_idle_add (BWChangedUI, this);
}

bool BlackWhite::BWComputed_ () {

	disableListener ();
	mixerRed->setValue (nextredbw);
	mixerGreen->setValue (nextgreenbw);
	mixerBlue->setValue (nextbluebw);
	enableListener ();

	updateRGBLabel();

	return false;
}

void BlackWhite::read (const ProcParams* pp, const ParamsEdited* pedited) {

	disableListener ();
	methodconn.block(true);
	//autoconn.block (true);
	filterconn.block(true);
	settingconn.block(true);
	enaccconn.block (true);
	enaconn.block (true);

	if (pedited && !pedited->blackwhite.setting)
		setting->set_active (15); // "Unchanged"
	else if (pp->blackwhite.setting=="NormalContrast")
		setting->set_active (0);
	else if (pp->blackwhite.setting=="HighContrast")
		setting->set_active (1);
	else if (pp->blackwhite.setting=="Luminance")
		setting->set_active (2);
	else if (pp->blackwhite.setting=="Landscape")
		setting->set_active (3);
	else if (pp->blackwhite.setting=="Portrait")
		setting->set_active (4);
	else if (pp->blackwhite.setting=="LowSensitivity")
		setting->set_active (5);
	else if (pp->blackwhite.setting=="HighSensitivity")
		setting->set_active (6);
	else if (pp->blackwhite.setting=="Panchromatic")
		setting->set_active (7);
	else if (pp->blackwhite.setting=="HyperPanchromatic")
		setting->set_active (8);
	else if (pp->blackwhite.setting=="Orthochromatic")
		setting->set_active (9);
	else if (pp->blackwhite.setting=="RGB-Abs")
		setting->set_active (10);
	else if (pp->blackwhite.setting=="RGB-Rel")
		setting->set_active (11);
	else if (pp->blackwhite.setting=="ROYGCBPM-Abs")
		setting->set_active (12);
	else if (pp->blackwhite.setting=="ROYGCBPM-Rel")
		setting->set_active (13);
	else if (pp->blackwhite.setting=="InfraRed")
		setting->set_active (14);
	settingChanged();


	if (pedited && !pedited->blackwhite.method)
		method->set_active (3); // "Unchanged"
	else if (pp->blackwhite.method=="Desaturation")
		method->set_active (0);
	else if (pp->blackwhite.method=="LumEqualizer")
		method->set_active (1);
	else if (pp->blackwhite.method=="ChannelMixer")
		method->set_active (2);
	methodChanged();

	
	if (pedited && !pedited->blackwhite.filter)
		filter->set_active (9); // "Unchanged"
	else if (pp->blackwhite.filter=="None")
		filter->set_active (0);
	else if (pp->blackwhite.filter=="Red")
		filter->set_active (1);
	else if (pp->blackwhite.filter=="Orange")
		filter->set_active (2);
	else if (pp->blackwhite.filter=="Yellow")
		filter->set_active (3);
	else if (pp->blackwhite.filter=="YellowGreen")
		filter->set_active (4);
	else if (pp->blackwhite.filter=="Green")
		filter->set_active (5);
	else if (pp->blackwhite.filter=="Cyan")
		filter->set_active (6);
	else if (pp->blackwhite.filter=="Blue")
		filter->set_active (7);
	else if (pp->blackwhite.filter=="Purple")
		filter->set_active (8);
	filterChanged();

	enabledcc->set_active (pp->blackwhite.enabledcc);
	lastEnabledcc = pp->blackwhite.enabledcc;
	enabled->set_active (pp->blackwhite.enabled);
	lastEnabled = pp->blackwhite.enabled;

	mixerRed->setValue (pp->blackwhite.mixerRed);
	mixerGreen->setValue (pp->blackwhite.mixerGreen);
	mixerBlue->setValue (pp->blackwhite.mixerBlue);
	gammaRed->setValue (pp->blackwhite.gammaRed);
	gammaGreen->setValue (pp->blackwhite.gammaGreen);
	gammaBlue->setValue (pp->blackwhite.gammaBlue);
	mixerOrange->setValue (pp->blackwhite.mixerOrange);
	mixerYellow->setValue (pp->blackwhite.mixerYellow);
	mixerCyan->setValue (pp->blackwhite.mixerCyan);
	mixerMagenta->setValue (pp->blackwhite.mixerMagenta);
	mixerPurple->setValue (pp->blackwhite.mixerPurple);
	luminanceCurve->setCurve (pp->blackwhite.luminanceCurve);
	beforeCurve->setCurve (pp->blackwhite.beforeCurve);
	beforeCurveMode->set_active(pp->blackwhite.beforeCurveMode);
	afterCurve->setCurve (pp->blackwhite.afterCurve);
//	afterCurveMode->set_active(pp->blackwhite.afterCurveMode);

	autoch->set_active (pp->blackwhite.autoc);
	lastAuto = pp->blackwhite.autoc;

	if (pedited) {
		luminanceCurve->setUnChanged (!pedited->blackwhite.luminanceCurve);
		beforeCurve->setUnChanged (!pedited->blackwhite.beforeCurve);
		afterCurve->setUnChanged (!pedited->blackwhite.afterCurve);
		autoch->set_inconsistent (!pedited->blackwhite.autoc);
		enabled->set_inconsistent  (!pedited->blackwhite.enabled);
		enabledcc->set_inconsistent  (!pedited->blackwhite.enabledcc);
		mixerRed->setEditedState (pedited->blackwhite.mixerRed ? Edited : UnEdited);
		mixerGreen->setEditedState (pedited->blackwhite.mixerGreen ? Edited : UnEdited);
		mixerBlue->setEditedState (pedited->blackwhite.mixerBlue ? Edited : UnEdited);
		gammaRed->setEditedState (pedited->blackwhite.gammaRed ? Edited : UnEdited);
		gammaGreen->setEditedState (pedited->blackwhite.gammaGreen ? Edited : UnEdited);
		gammaBlue->setEditedState (pedited->blackwhite.gammaBlue ? Edited : UnEdited);
		mixerOrange->setEditedState (pedited->blackwhite.mixerOrange ? Edited : UnEdited);
		mixerYellow->setEditedState (pedited->blackwhite.mixerYellow ? Edited : UnEdited);
		mixerCyan->setEditedState (pedited->blackwhite.mixerCyan ? Edited : UnEdited);
		mixerMagenta->setEditedState (pedited->blackwhite.mixerMagenta ? Edited : UnEdited);
		mixerPurple->setEditedState (pedited->blackwhite.mixerPurple ? Edited : UnEdited);
		if (!pedited->blackwhite.beforeCurveMode) {
			beforeCurveMode->set_active(4); // "Unchanged"
		}
//		if (!pedited->blackwhite.afterCurveMode) {
//			afterCurveMode->set_active(1); // "Unchanged"
//		}
	}
	methodconn.block(false);
	filterconn.block(false);
	settingconn.block(false);
	//autoconn.block (false);
	enaconn.block (false);
	enaccconn.block (false);

	updateRGBLabel();

	enableListener ();
}

void BlackWhite::write (ProcParams* pp, ParamsEdited* pedited) {
	pp->blackwhite.enabled = enabled->get_active ();
	pp->blackwhite.luminanceCurve = luminanceCurve->getCurve ();
	pp->blackwhite.autoc = autoch->get_active();
	pp->blackwhite.enabledcc = enabledcc->get_active ();
	pp->blackwhite.mixerRed = mixerRed->getValue ();
	pp->blackwhite.mixerGreen = mixerGreen->getValue ();
	pp->blackwhite.mixerBlue = mixerBlue->getValue ();
	pp->blackwhite.gammaRed = gammaRed->getValue ();
	pp->blackwhite.gammaGreen = gammaGreen->getValue ();
	pp->blackwhite.gammaBlue = gammaBlue->getValue ();
	pp->blackwhite.mixerOrange = mixerOrange->getValue ();
	pp->blackwhite.mixerYellow = mixerYellow->getValue ();
	pp->blackwhite.mixerCyan = mixerCyan->getValue ();
	pp->blackwhite.mixerMagenta = mixerMagenta->getValue ();
	pp->blackwhite.mixerPurple = mixerPurple->getValue ();
	pp->blackwhite.beforeCurve = beforeCurve->getCurve ();
	pp->blackwhite.afterCurve = afterCurve->getCurve ();

	int tcMode = beforeCurveMode->get_active_row_number();
	if      (tcMode == 0) pp->blackwhite.beforeCurveMode = BlackWhiteParams::TC_MODE_STD_BW;
	else if (tcMode == 1) pp->blackwhite.beforeCurveMode = BlackWhiteParams::TC_MODE_WEIGHTEDSTD_BW;
	else if (tcMode == 2) pp->blackwhite.beforeCurveMode = BlackWhiteParams::TC_MODE_FILMLIKE_BW;
	else if (tcMode == 3) pp->blackwhite.beforeCurveMode = BlackWhiteParams::TC_MODE_SATANDVALBLENDING_BW;

//	tcMode = afterCurveMode->get_active_row_number();
//	if      (tcMode == 0) pp->blackwhite.afterCurveMode = BlackWhiteParams::TC_MODE_STD_BW;
	//  else if (tcMode == 1) pp->blackwhite.afterCurveMode = BlackWhiteParams::TC_MODE_WEIGHTEDSTD;
	
	if (pedited) {
		pedited->blackwhite.enabled = !enabled->get_inconsistent();
		pedited->blackwhite.luminanceCurve = !luminanceCurve->isUnChanged ();
		pedited->blackwhite.autoc = !autoch->get_inconsistent();
		pedited->blackwhite.enabledcc = !enabledcc->get_inconsistent();
		pedited->blackwhite.mixerRed = mixerRed->getEditedState ();
		pedited->blackwhite.mixerGreen = mixerGreen->getEditedState ();
		pedited->blackwhite.mixerBlue = mixerBlue->getEditedState ();
		pedited->blackwhite.gammaRed = gammaRed->getEditedState ();
		pedited->blackwhite.gammaGreen = gammaGreen->getEditedState ();
		pedited->blackwhite.gammaBlue = gammaBlue->getEditedState ();
		pedited->blackwhite.filter = filter->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->blackwhite.setting = setting->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->blackwhite.method = method->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->blackwhite.mixerOrange = mixerOrange->getEditedState ();
		pedited->blackwhite.mixerYellow = mixerYellow->getEditedState ();
		pedited->blackwhite.mixerCyan = mixerCyan->getEditedState ();
		pedited->blackwhite.mixerMagenta = mixerMagenta->getEditedState ();
		pedited->blackwhite.mixerPurple = mixerPurple->getEditedState ();
		pedited->blackwhite.beforeCurve = !beforeCurve->isUnChanged ();
		pedited->blackwhite.beforeCurveMode = beforeCurveMode->get_active_row_number() != 4;
		pedited->blackwhite.afterCurve = !afterCurve->isUnChanged ();
//		pedited->blackwhite.afterCurveMode = afterCurveMode->get_active_row_number() != 1;
	}
	if (method->get_active_row_number()==0)
		pp->blackwhite.method = "Desaturation";
	else if (method->get_active_row_number()==1)
		pp->blackwhite.method = "LumEqualizer";
	else if (method->get_active_row_number()==2)
		pp->blackwhite.method = "ChannelMixer";

	pp->blackwhite.setting = getSettingString();
	pp->blackwhite.filter = getFilterString();
}

void BlackWhite::curveChanged (CurveEditor* ce) {
	if (listener) {
		if (ce == beforeCurve)
			listener->panelChanged (EvBWBeforeCurve, M("HISTORY_CUSTOMCURVE"));
		if (ce == afterCurve)
			listener->panelChanged (EvBWAfterCurve, M("HISTORY_CUSTOMCURVE"));
		if (ce == luminanceCurve)
			listener->panelChanged (EvBWLuminanceEqual, M("HISTORY_CUSTOMCURVE"));
	}
}

void BlackWhite::curveMode1Changed () {
	if (listener)  Glib::signal_idle().connect (sigc::mem_fun(*this, &BlackWhite::curveMode1Changed_));
}
bool BlackWhite::curveMode1Changed_ () {
	if (listener) listener->panelChanged (EvBWBeforeCurveMode, escapeHtmlChars(beforeCurveMode->get_active_text()));
	return false;
}
/*
void BlackWhite::curveMode1Changed2 () {
	if (listener)  Glib::signal_idle().connect (sigc::mem_fun(*this, &BlackWhite::curveMode1Changed2_));
}
bool BlackWhite::curveMode1Changed2_ () {
	if (listener) listener->panelChanged (EvBWAfterCurveMode, escapeHtmlChars(afterCurveMode->get_active_text()));
	return false;
}
*/
void BlackWhite::colorForValue (double valX, double valY, int callerId, ColorCaller* caller) {

	float r, g, b;

	if (callerId == 1) {        // Hue = f(Hue)

		float h = float((valY - 0.5) * 2. + valX);
		if (h > 1.0f)
			h -= 1.0f;
		else if (h < 0.0f)
			h += 1.0f;
		Color::hsv2rgb01(h, 0.5f, 0.5f, r, g, b);
		caller->ccRed = double(r);
		caller->ccGreen = double(g);
		caller->ccBlue = double(b);
	}
	else if (callerId == 2) {   // Saturation = f(Hue)
		Color::hsv2rgb01(float(valX), float(valY), 0.5f, r, g, b);
		caller->ccRed = double(r);
		caller->ccGreen = double(g);
		caller->ccBlue = double(b);
	}
	else if (callerId == 3) {   // Value = f(Hue)
		Color::hsv2rgb01(float(valX), 0.5f, float(valY), r, g, b);
		caller->ccRed = double(r);
		caller->ccGreen = double(g);
		caller->ccBlue = double(b);
	}
	else {
		printf("Error: no curve displayed!\n");
	}
}


void BlackWhite::settingChanged () {

	if ( setting->get_active_row_number()==10 || setting->get_active_row_number()==11 ) {
		// RGB Channel Mixer
		showMixer(3);
		hideEnabledCC();
		showFilter();
	}
	else if ( setting->get_active_row_number()==12 || setting->get_active_row_number()==13 ) {
		// ROYGCBPM Channel Mixer
		showMixer(7);
		showEnabledCC();
		showFilter();
	}
	else if ( setting->get_active_row_number()==14 ) {
		// Infrared
		filter->set_active (0);
		showMixer(3, false);
		hideEnabledCC();
		hideFilter();
	}
	else {
		// RGB Presets
		showMixer(3, false);
		hideEnabledCC();
		showFilter();
	}

	// Checking "listener" to avoid "autoch" getting toggled off because it has to change the sliders when toggling on
	if (listener){
		if (multiImage && autoch->get_inconsistent())
			autoch->set_inconsistent (false);
		autoconn.block(true);
		autoch->set_active (false);
		autoconn.block(false);
		lastAuto = false;
	}

	updateRGBLabel();

	if (listener && (multiImage||enabled->get_active())) {
		listener->panelChanged (EvBWsetting, setting->get_active_text ());
	}
}


void BlackWhite::filterChanged () {
	// Checking "listener" to avoid "autoch" getting toggled off because it has to change the sliders when toggling on
	if (listener){
		if (multiImage && autoch->get_inconsistent())
			autoch->set_inconsistent (false);
		autoconn.block(true);
		autoch->set_active (false);
		autoconn.block(false);
		lastAuto = false;
	}

	updateRGBLabel();

	if (listener && (multiImage||enabled->get_active())) {	
		listener->panelChanged (EvBWfilter, filter->get_active_text ());
	}
}

void BlackWhite::methodChanged () {
	if(method->get_active_row_number()==2) {
		// Channel Mixer
		hideLuminance();

		if(setting->get_active_row_number()==10 || setting->get_active_row_number()==11){
			hideEnabledCC();
			showMixer(3);
		}
		else if(setting->get_active_row_number()==12 || setting->get_active_row_number()==13) {
			showEnabledCC();
			showMixer(7);
		}
		else {
			hideEnabledCC();
			showMixer(3, false);
		}
		beforeCurveCEG->show();
		afterCurveCEG->show();

		bool wasEnabled = disableListener();
		settingChanged();
		if (wasEnabled) enableListener();
	}
	 else if(method->get_active_row_number()==1) {
		// Luminance Equalizer
		showLuminance();
		hideMixer();
		beforeCurveCEG->show();
		afterCurveCEG->show();
	}
	else if(method->get_active_row_number()==0) {
		// Desaturation
		hideLuminance();
		hideMixer();
		beforeCurveCEG->show();
		afterCurveCEG->show();
	}
	if (listener && (multiImage||enabled->get_active())) {
		listener->panelChanged (EvBWmethod, method->get_active_text ());
	}
}

void BlackWhite::enabled_toggled () {

	if (multiImage) {
		if (enabled->get_inconsistent()) {
			enabled->set_inconsistent (false);
			enaconn.block (true);
			enabled->set_active (false);
			enaconn.block (false);
		}
		else if (lastEnabled)
			enabled->set_inconsistent (true);

		lastEnabled = enabled->get_active ();
	}

	if (listener) {
		if (enabled->get_inconsistent())
			listener->panelChanged (EvBWChmixEnabled, M("GENERAL_UNCHANGED"));
		else if (enabled->get_active ()) {
			listener->panelChanged (EvBWChmixEnabled, M("GENERAL_ENABLED"));
		}
		else {
			listener->panelChanged (EvBWChmixEnabled, M("GENERAL_DISABLED"));
		}
	}
}

void BlackWhite::neutral_pressed () {
	// This method deselects auto chmixer and sets "neutral" values to params
	disableListener();

	if (multiImage && autoch->get_inconsistent())
		autoch->set_inconsistent (false);
	autoconn.block(true);
	autoch->set_active (false);
	autoconn.block(false);
	lastAuto = false;

	int activeSetting = setting->get_active_row_number();
	if (activeSetting < 10 || activeSetting > 13)
		setting->set_active (11);
	filter->set_active (0);
	mixerRed->resetValue(false);
	mixerGreen->resetValue(false);
	mixerBlue->resetValue(false);
	mixerOrange->resetValue(false);
	mixerYellow->resetValue(false);
	mixerMagenta->resetValue(false);
	mixerPurple->resetValue(false);
	mixerCyan->resetValue(false);

	enableListener();

	updateRGBLabel();

	listener->panelChanged (EvNeutralBW, M("ADJUSTER_RESET_TO_DEFAULT"));
}

void BlackWhite::enabledcc_toggled () {

	// toggling off the Complementary Colors does switch off the Auto button
	if (multiImage) {
		// multiple image editing (batch)
		if (enabledcc->get_inconsistent()) {
			enabledcc->set_inconsistent (false);  // set consistent
			enaccconn.block (true);
			enabledcc->set_active (false);  // ... and deactivated
			enaccconn.block (false);
		}
		else if (lastEnabledcc)
			enabledcc->set_inconsistent (true);

		lastEnabledcc = enabledcc->get_active ();
	}

	if (multiImage && autoch->get_inconsistent())
		autoch->set_inconsistent (false);
	autoconn.block(true);
	autoch->set_active (false);
	autoconn.block(false);
	lastAuto = false;

	updateRGBLabel();

	if (listener) {
		if (enabledcc->get_inconsistent())
			listener->panelChanged (EvBWChmixEnabledLm, M("GENERAL_UNCHANGED"));
		else if (enabledcc->get_active ()) {
			listener->panelChanged (EvBWChmixEnabledLm, M("GENERAL_ENABLED"));
		}
		else {
			listener->panelChanged (EvBWChmixEnabledLm, M("GENERAL_DISABLED"));
		}
	}
}


void BlackWhite::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

	mixerRed->setDefault (defParams->blackwhite.mixerRed);
	mixerGreen->setDefault (defParams->blackwhite.mixerGreen);
	mixerBlue->setDefault (defParams->blackwhite.mixerBlue);
	gammaRed->setDefault (defParams->blackwhite.gammaRed);
	gammaGreen->setDefault (defParams->blackwhite.gammaGreen);
	gammaBlue->setDefault (defParams->blackwhite.gammaBlue);
	mixerOrange->setDefault (defParams->blackwhite.mixerOrange);
	mixerYellow->setDefault (defParams->blackwhite.mixerYellow);
	mixerCyan->setDefault (defParams->blackwhite.mixerCyan);
	mixerMagenta->setDefault (defParams->blackwhite.mixerMagenta);
	mixerPurple->setDefault (defParams->blackwhite.mixerPurple);

	if (pedited) {
		mixerRed->setDefaultEditedState (pedited->blackwhite.mixerRed ? Edited : UnEdited);
		mixerGreen->setDefaultEditedState (pedited->blackwhite.mixerGreen ? Edited : UnEdited);
		mixerBlue->setDefaultEditedState (pedited->blackwhite.mixerBlue ? Edited : UnEdited);
		gammaRed->setDefaultEditedState (pedited->blackwhite.gammaRed ? Edited : UnEdited);
		gammaGreen->setDefaultEditedState (pedited->blackwhite.gammaGreen ? Edited : UnEdited);
		gammaBlue->setDefaultEditedState (pedited->blackwhite.gammaBlue ? Edited : UnEdited);
		mixerOrange->setDefaultEditedState (pedited->blackwhite.mixerOrange ? Edited : UnEdited);
		mixerYellow->setDefaultEditedState (pedited->blackwhite.mixerYellow ? Edited : UnEdited);
		mixerCyan->setDefaultEditedState (pedited->blackwhite.mixerCyan ? Edited : UnEdited);
		mixerMagenta->setDefaultEditedState (pedited->blackwhite.mixerMagenta ? Edited : UnEdited);
		mixerPurple->setDefaultEditedState (pedited->blackwhite.mixerPurple ? Edited : UnEdited);
	}
	else {
		mixerRed->setDefaultEditedState (Irrelevant);
		mixerGreen->setDefaultEditedState (Irrelevant);
		mixerBlue->setDefaultEditedState (Irrelevant);
		gammaRed->setDefaultEditedState (Irrelevant);
		gammaGreen->setDefaultEditedState (Irrelevant);
		gammaBlue->setDefaultEditedState (Irrelevant);	
		mixerOrange->setDefaultEditedState (Irrelevant);
		mixerYellow->setDefaultEditedState (Irrelevant);
		mixerCyan->setDefaultEditedState (Irrelevant);
		mixerMagenta->setDefaultEditedState (Irrelevant);
		mixerPurple->setDefaultEditedState (Irrelevant);
	}
}

void BlackWhite::autoch_toggled () {

	if (batchMode) {
		if (multiImage) {
			if (autoch->get_inconsistent()) {
				autoch->set_inconsistent (false);
				autoconn.block (true);
				autoch->set_active (false);
				autoconn.block (false);
			}
			else if (lastAuto)
				autoch->set_inconsistent (true);
		}

		lastAuto = autoch->get_active ();

		mixerRed->setEditedState (UnEdited);
		mixerGreen->setEditedState (UnEdited);
		mixerBlue->setEditedState (UnEdited);
		mixerOrange->setEditedState (UnEdited);
		mixerYellow->setEditedState (UnEdited);
		mixerPurple->setEditedState (UnEdited);
		mixerMagenta->setEditedState (UnEdited);
		mixerCyan->setEditedState (UnEdited);
		
		bool wasEnabled = disableListener();
		if (mixerRed->getAddMode())
			mixerRed->resetValue(false);
		if (mixerGreen->getAddMode())
			mixerGreen->resetValue(true);
		if (mixerBlue->getAddMode())
			mixerBlue->resetValue(true);
		if (mixerOrange->getAddMode())
			mixerOrange->resetValue(true);
		if (mixerYellow->getAddMode())
			mixerYellow->resetValue(true);
		if (mixerMagenta->getAddMode())
			mixerMagenta->resetValue(true);
		if (mixerPurple->getAddMode())
			mixerPurple->resetValue(true);
		if (mixerCyan->getAddMode())
			mixerCyan->resetValue(true);
		setting->set_active (11);
		filter->set_active (0);
		if (wasEnabled) enableListener();

		if (listener) {
			if (autoch->get_inconsistent())
				listener->panelChanged (EvAutoch, M("GENERAL_UNCHANGED"));
			else if (autoch->get_active ())
				listener->panelChanged (EvAutoch, M("GENERAL_ENABLED"));
			else
				listener->panelChanged (EvAutoch, M("GENERAL_DISABLED"));
		}
	}
	else {
		if (autoch->get_active()) {
			bool wasEnabled = disableListener();
			mixerRed->setValue(33);
			mixerGreen->setValue(33);
			mixerBlue->setValue(33);
			mixerOrange->setValue(33);
			mixerYellow->setValue(33);
			mixerMagenta->setValue(33);
			mixerPurple->setValue(33);
			mixerCyan->setValue(33);
			setting->set_active (11);
			filter->set_active (0);
			if (wasEnabled) enableListener();

			updateRGBLabel();

			if (listener)
				listener->panelChanged (EvAutoch, M("GENERAL_ENABLED"));
		}
		else {
			if (listener)
				listener->panelChanged (EvAutoch, M("GENERAL_DISABLED"));
		}
	}

}

void BlackWhite::adjusterChanged (Adjuster* a, double newval) {

	// Checking "listener" to avoid "autoch" getting toggled off because it has to change the sliders when toggling on
	if (listener && (a==mixerRed || a==mixerGreen || a==mixerBlue || a==mixerOrange || a==mixerYellow || a==mixerMagenta || a==mixerPurple || a==mixerCyan) ) {
		if (multiImage && autoch->get_inconsistent())
			autoch->set_inconsistent (false);
		autoconn.block(true);
		autoch->set_active (false);
		autoconn.block(false);
		lastAuto = false;
	}

	if (a == mixerRed || a==mixerGreen || a== mixerBlue
		|| a == mixerOrange || a == mixerYellow ||a == mixerCyan
		|| a == mixerMagenta || a == mixerPurple)

		updateRGBLabel();

	if (listener  && (multiImage||enabled->get_active())) {
		Glib::ustring value = a->getTextValue();
		if (a == mixerRed)
			listener->panelChanged (EvBWred, value );
		else if (a == mixerGreen)
			listener->panelChanged (EvBWgreen, value );
		else if (a == mixerBlue)
			listener->panelChanged (EvBWblue, value );
		else if (a == gammaGreen)
			listener->panelChanged (EvBWgreengam, value );
		else if (a == gammaBlue)
			listener->panelChanged (EvBWbluegam, value );
		else if (a == gammaRed)
			listener->panelChanged (EvBWredgam, value );
		else if (a == mixerOrange)
			listener->panelChanged (EvBWoran, value );
		else if (a == mixerYellow)
			listener->panelChanged (EvBWyell, value );
		else if (a == mixerCyan)
			listener->panelChanged (EvBWcyan, value );
		else if (a == mixerMagenta)
			listener->panelChanged (EvBWmag, value );
		else if (a == mixerPurple)
			listener->panelChanged (EvBWpur, value );
	}
}

void BlackWhite::updateRGBLabel () {
	if (!batchMode) {
		float kcorrec=1.f;
		float r, g, b;
		if (autoch->get_active()) {
			r = nextredbw;
			g = nextgreenbw;
			b = nextbluebw;
		}
		else {
			r = mixerRed->getValue();
			g = mixerGreen->getValue();
			b = mixerBlue->getValue();
		}
		double mixR, mixG, mixB;
		Glib::ustring sSetting = getSettingString();
		Color::computeBWMixerConstants(sSetting, getFilterString(), r, g, b,
				mixerOrange->getValue(), mixerYellow->getValue(), mixerCyan->getValue(), mixerPurple->getValue(), mixerMagenta->getValue(),
				autoch->get_active(), enabledcc->get_active(), kcorrec, mixR, mixG, mixB);
		RGBLabels->set_text(
				Glib::ustring::compose(M("TP_BWMIX_RGBLABEL"),
				Glib::ustring::format(std::fixed, std::setprecision(1), r*100.),
				Glib::ustring::format(std::fixed, std::setprecision(1), g*100.),
				Glib::ustring::format(std::fixed, std::setprecision(1), b*100.),
				Glib::ustring::format(std::fixed, std::setprecision(0), ceil(kcorrec*100./*(r+g+b)*100.)*/)))
				);
		// We have to update the RGB sliders too if preset values has been chosen
		if (sSetting != "RGB-Abs" && sSetting != "RGB-Rel" && sSetting != "ROYGCBPM-Abs" && sSetting != "ROYGCBPM-Rel") {
			mixerRed->setValue(mixR);
			mixerGreen->setValue(mixG);
			mixerBlue->setValue(mixB);
		}
	}
}

void BlackWhite::setBatchMode (bool batchMode) {
	removeIfThere (autoHBox, autoch, false);
	autoch = Gtk::manage (new Gtk::CheckButton (M("TP_BWMIX_AUTOCH")));
	autoch->set_tooltip_markup (M("TP_BWMIX_AUTOCH_TIP"));
	autoconn = autoch->signal_toggled().connect( sigc::mem_fun(*this, &BlackWhite::autoch_toggled) );
	autoHBox->pack_start (*autoch);

	removeIfThere (mixerVBox, RGBLabels, false);
	delete RGBLabels;
	RGBLabels = NULL;

	ToolPanel::setBatchMode (batchMode);
	mixerRed->showEditedCB ();
	mixerOrange->showEditedCB ();
	mixerYellow->showEditedCB ();
	mixerGreen->showEditedCB ();
	mixerCyan->showEditedCB ();
	mixerBlue->showEditedCB ();
	mixerMagenta->showEditedCB ();
	mixerPurple->showEditedCB ();
	gammaRed->showEditedCB ();
	gammaGreen->showEditedCB ();
	gammaBlue->showEditedCB ();
	method->append_text (M("GENERAL_UNCHANGED"));
	filter->append_text (M("GENERAL_UNCHANGED"));
	setting->append_text (M("GENERAL_UNCHANGED"));
	luminanceCEG->setBatchMode (batchMode);
	beforeCurveCEG->setBatchMode (batchMode);
	beforeCurveMode->append_text (M("GENERAL_UNCHANGED"));
	afterCurveCEG->setBatchMode (batchMode);
//	afterCurveMode->append_text (M("GENERAL_UNCHANGED"));

	showLuminance();
	showFilter();
	showEnabledCC();
	showGamma();
	showMixer(7);
}

void BlackWhite::autoOpenCurve () {
	luminanceCurve->openIfNonlinear();
	beforeCurve->openIfNonlinear();
	afterCurve->openIfNonlinear();
}

void BlackWhite::setAdjusterBehavior (bool bwadd, bool bwgadd) {

	mixerRed->setAddMode(bwadd);
	mixerOrange->setAddMode(bwadd);
	mixerYellow->setAddMode(bwadd);
	mixerGreen->setAddMode(bwadd);
	mixerCyan->setAddMode(bwadd);
	mixerBlue->setAddMode(bwadd);
	mixerMagenta->setAddMode(bwadd);
	mixerPurple->setAddMode(bwadd);

	gammaRed->setAddMode(bwgadd);
	gammaGreen->setAddMode(bwgadd);
	gammaBlue->setAddMode(bwgadd);
}

void BlackWhite::trimValues (rtengine::procparams::ProcParams* pp) {

	mixerRed->trimValue (pp->blackwhite.mixerRed);
	mixerGreen->trimValue (pp->blackwhite.mixerGreen);
	mixerBlue->trimValue (pp->blackwhite.mixerBlue);
	gammaRed->trimValue (pp->blackwhite.gammaRed);
	gammaGreen->trimValue (pp->blackwhite.gammaGreen);
	gammaBlue->trimValue (pp->blackwhite.gammaBlue);
	mixerOrange->trimValue (pp->blackwhite.mixerOrange);
	mixerYellow->trimValue (pp->blackwhite.mixerYellow);
	mixerCyan->trimValue (pp->blackwhite.mixerCyan);
	mixerMagenta->trimValue (pp->blackwhite.mixerMagenta);
	mixerPurple->trimValue (pp->blackwhite.mixerPurple);
}

void BlackWhite::showLuminance() {
	luminanceCEG->show();
	luminanceSep->show();
}

void BlackWhite::hideLuminance() {
	if (!batchMode) {
		luminanceCEG->hide();
		luminanceSep->hide();
	}
}

void BlackWhite::showFilter() {
	filterHBox->show();
	filterSep->show();
}

void BlackWhite::hideFilter() {
	if (!batchMode) {
		filterHBox->hide();
		filterSep->hide();
	}
}

void BlackWhite::showEnabledCC() {
	enabledcc->show();
	enabledccSep->show();
}

void BlackWhite::hideEnabledCC() {
	if (!batchMode) {
		enabledcc->hide();
		enabledccSep->hide();
	}
}

void BlackWhite::showMixer(int nChannels, bool RGBIsSensitive) {
	if (!batchMode)
		RGBLabels->show();

	if (!batchMode && nChannels == 3) {
		mixerRed->show();   mixerRed->set_sensitive (RGBIsSensitive);
		mixerGreen->show(); mixerGreen->set_sensitive (RGBIsSensitive);
		mixerBlue->show();  mixerBlue->set_sensitive (RGBIsSensitive);
		filterSep2->hide();
		mixerOrange->hide();
		mixerYellow->hide();
		mixerCyan->hide();
		mixerMagenta->hide();
		mixerPurple->hide();
	}
	else {
		mixerRed->show();   mixerRed->set_sensitive (true);
		mixerGreen->show(); mixerGreen->set_sensitive (true);
		mixerBlue->show();  mixerBlue->set_sensitive (true);
		filterSep2->show();
		mixerOrange->show();
		mixerYellow->show();
		mixerCyan->show();
		mixerMagenta->show();
		mixerPurple->show();
	}
	mixerFrame->show();
}

void BlackWhite::hideMixer() {
	if (!batchMode)
		mixerFrame->hide();
}

void BlackWhite::showGamma() {
	gammaFrame->show();
}

void BlackWhite::hideGamma() {
	if (!batchMode)
		gammaFrame->hide();
}

Glib::ustring BlackWhite::getSettingString() {
	Glib::ustring retVal;
	if (setting->get_active_row_number()==0)
		retVal = "NormalContrast";
	else if (setting->get_active_row_number()==1)
		retVal = "HighContrast";
	else if (setting->get_active_row_number()==2)
		retVal = "Luminance";
	else if (setting->get_active_row_number()==3)
		retVal = "Landscape";
	else if (setting->get_active_row_number()==4)
		retVal = "Portrait";
	else if (setting->get_active_row_number()==5)
		retVal = "LowSensitivity";
	else if (setting->get_active_row_number()==6)
		retVal = "HighSensitivity";
	else if (setting->get_active_row_number()==7)
		retVal = "Panchromatic";
	else if (setting->get_active_row_number()==8)
		retVal = "HyperPanchromatic";
	else if (setting->get_active_row_number()==9)
		retVal = "Orthochromatic";
	else if (setting->get_active_row_number()==10)
		retVal = "RGB-Abs";
	else if (setting->get_active_row_number()==11)
		retVal = "RGB-Rel";
	else if (setting->get_active_row_number()==12)
		retVal = "ROYGCBPM-Abs";
	else if (setting->get_active_row_number()==13)
		retVal = "ROYGCBPM-Rel";
	else if (setting->get_active_row_number()==14)
		retVal = "InfraRed";
	return retVal;
}

Glib::ustring BlackWhite::getFilterString() {
	Glib::ustring retVal;
	if (filter->get_active_row_number()==0)
		retVal = "None";
	else if (filter->get_active_row_number()==1)
		retVal = "Red";
	else if (filter->get_active_row_number()==2)
		retVal = "Orange";
	else if (filter->get_active_row_number()==3)
		retVal = "Yellow";
	else if (filter->get_active_row_number()==4)
		retVal = "YellowGreen";
	else if (filter->get_active_row_number()==5)
		retVal = "Green";
	else if (filter->get_active_row_number()==6)
		retVal = "Cyan";
	else if (filter->get_active_row_number()==7)
		retVal = "Blue";
	else if (filter->get_active_row_number()==8)
		retVal = "Purple";
	return retVal;
}
