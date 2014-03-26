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
 */
#include "colorappearance.h"
#include <cmath>
#include <iomanip>
#include "guiutils.h"
#include "../rtengine/color.h"

using namespace rtengine;
using namespace rtengine::procparams;

ColorAppearance::ColorAppearance () : FoldableToolPanel(this) {
	CurveListener::setMulti(true);
	std::vector<GradientMilestone> milestones;
	milestones.push_back( GradientMilestone(0., 0., 0., 0.) );
	milestones.push_back( GradientMilestone(1., 1., 1., 1.) );

	enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
	enabled->set_active (false);
	pack_start (*enabled);


	// ------------------------ Process #1: Converting to CIECAM


	// Process 1 frame
	Gtk::Frame *p1Frame;
	// Vertical box container for the content of the Process 1 frame
	Gtk::VBox *p1VBox;

	p1Frame = Gtk::manage (new Gtk::Frame(M("TP_COLORAPP_LABEL_SCENE")) );
	p1Frame->set_border_width(0);
	p1Frame->set_label_align(0.025, 0.5);

	p1VBox = Gtk::manage ( new Gtk::VBox());
	p1VBox->set_border_width(4);
	p1VBox->set_spacing(2);

	degree  = Gtk::manage (new Adjuster (M("TP_COLORAPP_CIECAT_DEGREE"),    0.,  100.,  1.,   100.));
	if (degree->delay < 1000) degree->delay = 1000;
	degree->throwOnButtonRelease();
	degree->addAutoButton(M("TP_COLORAPP_DEGREE_AUTO_TOOLTIP"));
	degree->set_tooltip_markup (M("TP_COLORAPP_DEGREE_TOOLTIP"));
	p1VBox->pack_start (*degree);

	surrsource = Gtk::manage (new Gtk::CheckButton (M("TP_COLORAPP_SURSOURCE")));
	surrsource->set_tooltip_markup (M("TP_COLORAPP_SURSOURCE_TOOLTIP"));
	p1VBox->pack_start (*surrsource, Gtk::PACK_SHRINK);

	Gtk::HBox* wbmHBox = Gtk::manage (new Gtk::HBox ());
	wbmHBox->set_border_width (0);
	wbmHBox->set_spacing (2);
	wbmHBox->set_tooltip_markup (M("TP_COLORAPP_MODEL_TOOLTIP"));
	Gtk::Label* wbmLab = Gtk::manage (new Gtk::Label (M("TP_COLORAPP_MODEL")+":"));
	wbmHBox->pack_start (*wbmLab, Gtk::PACK_SHRINK);
	wbmodel = Gtk::manage (new MyComboBoxText ());
	wbmodel->append_text (M("TP_COLORAPP_WBRT"));
	wbmodel->append_text (M("TP_COLORAPP_WBCAM"));
	wbmodel->set_active (0);
	wbmHBox->pack_start (*wbmodel);
	p1VBox->pack_start (*wbmHBox);

	adapscen = Gtk::manage (new Adjuster (M("TP_COLORAPP_ADAPTSCENE"), 0.001, 16384., 0.001, 2000.));// EV -7  ==> EV 17
	if (adapscen->delay < 1000) adapscen->delay = 1000;
	adapscen->throwOnButtonRelease();
	adapscen->addAutoButton(M("TP_COLORAPP_ADAP_AUTO_TOOLTIP"));
	adapscen->set_tooltip_markup (M("TP_COLORAPP_ADAPTSCENE_TOOLTIP"));
	p1VBox->pack_start (*adapscen);

	p1Frame->add(*p1VBox);
	pack_start (*p1Frame, Gtk::PACK_EXPAND_WIDGET, 4);


	// ------------------------ Process #2: Modifying image inside CIECAM


	// Process 1 frame
	Gtk::Frame *p2Frame;
	// Vertical box container for the content of the Process 1 frame
	Gtk::VBox *p2VBox;

	p2Frame = Gtk::manage (new Gtk::Frame(M("TP_COLORAPP_LABEL_CAM02")) );
	p2Frame->set_border_width(0);
	p2Frame->set_label_align(0.025, 0.5);

	p2VBox = Gtk::manage ( new Gtk::VBox());
	p2VBox->set_border_width(4);
	p2VBox->set_spacing(2);

	Gtk::HBox* alHBox = Gtk::manage (new Gtk::HBox ());
	alHBox->set_border_width (0);
	alHBox->set_spacing (2);
	alHBox->set_tooltip_markup (M("TP_COLORAPP_ALGO_TOOLTIP"));
	Gtk::Label* alLabel = Gtk::manage (new Gtk::Label (M("TP_COLORAPP_ALGO")+":"));
	alHBox->pack_start (*alLabel, Gtk::PACK_SHRINK);
	algo = Gtk::manage (new MyComboBoxText ());
	algo->append_text (M("TP_COLORAPP_ALGO_JC"));
	algo->append_text (M("TP_COLORAPP_ALGO_JS"));
	algo->append_text (M("TP_COLORAPP_ALGO_QM"));
	algo->append_text (M("TP_COLORAPP_ALGO_ALL"));
	algo->set_active (0);
	alHBox->pack_start (*algo);
	p2VBox->pack_start (*alHBox);

	p2VBox->pack_start (*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_EXPAND_WIDGET, 4);

	jlight = Gtk::manage (new Adjuster (M("TP_COLORAPP_LIGHT"), -100.0, 100.0, 0.1, 0.));
	if (jlight->delay < 1000) jlight->delay = 1000;
	jlight->throwOnButtonRelease();
	jlight->set_tooltip_markup (M("TP_COLORAPP_LIGHT_TOOLTIP"));
	p2VBox->pack_start (*jlight);

	qbright = Gtk::manage (new Adjuster (M("TP_COLORAPP_BRIGHT"), -100.0, 100.0, 0.1, 0.));
	if (qbright->delay < 1000) qbright->delay = 1000;
	qbright->throwOnButtonRelease();
	qbright->set_tooltip_markup (M("TP_COLORAPP_BRIGHT_TOOLTIP"));
	p2VBox->pack_start (*qbright);

	chroma = Gtk::manage (new Adjuster (M("TP_COLORAPP_CHROMA"), -100.0, 100.0, 0.1, 0.));
	if (chroma->delay < 1000) chroma->delay = 1000;
	chroma->throwOnButtonRelease();
	chroma->set_tooltip_markup (M("TP_COLORAPP_CHROMA_TOOLTIP"));
	p2VBox->pack_start (*chroma);


	schroma = Gtk::manage (new Adjuster (M("TP_COLORAPP_CHROMA_S"), -100.0, 100.0, 0.1, 0.));
	if (schroma->delay < 1000) schroma->delay = 1000;
	schroma->throwOnButtonRelease();
	schroma->set_tooltip_markup (M("TP_COLORAPP_CHROMA_S_TOOLTIP"));
	p2VBox->pack_start (*schroma);

	mchroma = Gtk::manage (new Adjuster (M("TP_COLORAPP_CHROMA_M"), -100.0, 100.0, 0.1, 0.));
	if (mchroma->delay < 1000) mchroma->delay = 1000;
	mchroma->throwOnButtonRelease();
	mchroma->set_tooltip_markup (M("TP_COLORAPP_CHROMA_M_TOOLTIP"));
	p2VBox->pack_start (*mchroma);
	
	rstprotection = Gtk::manage ( new Adjuster (M("TP_COLORAPP_RSTPRO"), 0., 100., 0.1, 0.) );
	if (rstprotection->delay < 1000) rstprotection->delay = 1000;
	rstprotection->throwOnButtonRelease();
	rstprotection->set_tooltip_markup (M("TP_COLORAPP_RSTPRO_TOOLTIP"));
	p2VBox->pack_start (*rstprotection);

	contrast = Gtk::manage (new Adjuster (M("TP_COLORAPP_CONTRAST"), -100.0, 100.0, 0.1, 0.));
	if (contrast->delay < 1000) contrast->delay = 1000;
	contrast->throwOnButtonRelease();
	contrast->set_tooltip_markup (M("TP_COLORAPP_CONTRAST_TOOLTIP"));
	p2VBox->pack_start (*contrast);

	qcontrast = Gtk::manage (new Adjuster (M("TP_COLORAPP_CONTRAST_Q"), -100.0, 100.0, 0.1, 0.));
	if (qcontrast->delay < 1000) qcontrast->delay = 1000;
	qcontrast->throwOnButtonRelease();
	qcontrast->set_tooltip_markup (M("TP_COLORAPP_CONTRAST_Q_TOOLTIP"));
	p2VBox->pack_start (*qcontrast);

	
	colorh = Gtk::manage (new Adjuster (M("TP_COLORAPP_HUE"), -100.0, 100.0, 0.1, 0.));
	if (colorh->delay < 1000) colorh->delay = 1000;
	colorh->throwOnButtonRelease();
	colorh->set_tooltip_markup (M("TP_COLORAPP_HUE_TOOLTIP"));
	p2VBox->pack_start (*colorh);

	tonecie = Gtk::manage (new Gtk::CheckButton (M("TP_COLORAPP_TONECIE")));
	tonecie->set_tooltip_markup (M("TP_COLORAPP_TONECIE_TOOLTIP"));
	tonecieconn = tonecie->signal_toggled().connect( sigc::mem_fun(*this, &ColorAppearance::tonecie_toggled) );
	p2VBox->pack_start (*tonecie);
/*
	sharpcie = Gtk::manage (new Gtk::CheckButton (M("TP_COLORAPP_SHARPCIE")));
	sharpcie->set_tooltip_markup (M("TP_COLORAPP_SHARPCIE_TOOLTIP"));
	sharpcieconn = sharpcie->signal_toggled().connect( sigc::mem_fun(*this, &ColorAppearance::sharpcie_toggled) );
	p2VBox->pack_start (*sharpcie);
*/	
	p2VBox->pack_start (*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_EXPAND_WIDGET, 4);

	toneCurveMode = Gtk::manage (new MyComboBoxText ());
	toneCurveMode->append_text (M("TP_COLORAPP_TCMODE_LIGHTNESS"));
	toneCurveMode->append_text (M("TP_COLORAPP_TCMODE_BRIGHTNESS"));
	toneCurveMode->set_active (0);
	toneCurveMode->set_tooltip_text(M("TP_COLORAPP_TCMODE_LABEL1"));

	curveEditorG = new CurveEditorGroup (options.lastToneCurvesDir, M("TP_COLORAPP_CURVEEDITOR1"));
	curveEditorG->setCurveListener (this);
	curveEditorG->setTooltip(M("TP_COLORAPP_CURVEEDITOR1_TOOLTIP"));

	shape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "", toneCurveMode));



	tcmodeconn = toneCurveMode->signal_changed().connect( sigc::mem_fun(*this, &ColorAppearance::curveMode1Changed), true );

	toneCurveMode2 = Gtk::manage (new MyComboBoxText ());
	toneCurveMode2->append_text (M("TP_COLORAPP_TCMODE_LIGHTNESS"));
	toneCurveMode2->append_text (M("TP_COLORAPP_TCMODE_BRIGHTNESS"));
	toneCurveMode2->set_active (0);
	toneCurveMode2->set_tooltip_text(M("TP_COLORAPP_TCMODE_LABEL2"));

	curveEditorG2 = new CurveEditorGroup (options.lastToneCurvesDir, M("TP_COLORAPP_CURVEEDITOR2"));
	curveEditorG2->setCurveListener (this);

	shape2 = static_cast<DiagonalCurveEditor*>(curveEditorG2->addCurve(CT_Diagonal, "", toneCurveMode2));

	tcmode2conn = toneCurveMode2->signal_changed().connect( sigc::mem_fun(*this, &ColorAppearance::curveMode2Changed), true );

	toneCurveMode3 = Gtk::manage (new MyComboBoxText ());
	toneCurveMode3->append_text (M("TP_COLORAPP_TCMODE_CHROMA"));
	toneCurveMode3->append_text (M("TP_COLORAPP_TCMODE_SATUR"));
	toneCurveMode3->append_text (M("TP_COLORAPP_TCMODE_COLORF"));
	toneCurveMode3->set_active (0);
	toneCurveMode3->set_tooltip_text(M("TP_COLORAPP_TCMODE_LABEL3"));

	curveEditorG3 = new CurveEditorGroup (options.lastToneCurvesDir, M("TP_COLORAPP_CURVEEDITOR3"));
	curveEditorG3->setCurveListener (this);

	shape3 = static_cast<DiagonalCurveEditor*>(curveEditorG3->addCurve(CT_Diagonal, "", toneCurveMode3));
	shape3->setRangeLabels(
			M("TP_LABCURVE_CURVEEDITOR_CC_RANGE1"), M("TP_LABCURVE_CURVEEDITOR_CC_RANGE2"),
			M("TP_LABCURVE_CURVEEDITOR_CC_RANGE3"), M("TP_LABCURVE_CURVEEDITOR_CC_RANGE4")
	);
	shape3->setBottomBarColorProvider(this, 1);
	shape3->setLeftBarColorProvider(this, 1);
	shape3->setRangeDefaultMilestones(0.05, 0.2, 0.58);

	
//	shape3->setBottomBarColorProvider(this, 2);
//	shape3->setLeftBarColorProvider(this, 2);
//	shape3->setRangeDefaultMilestones(0.05, 0.2, 0.58);
	
	// The milestones are still the same than those define above
	//milestones.push_back( GradientMilestone(0., 0., 0., 0.) );
	//milestones.push_back( GradientMilestone(1., 1., 1., 1.) );
	shape->setBottomBarBgGradient(milestones);
	shape->setLeftBarBgGradient(milestones);
	shape2->setBottomBarBgGradient(milestones);
	shape2->setLeftBarBgGradient(milestones);
	
	std::vector<GradientMilestone> shape3Milestones;
	float R, G, B;
	for (int i=0; i<7; i++) {
		float x = float(i)*(1.0f/6.0);
		Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
		shape3Milestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
	}
	shape3->setBottomBarBgGradient(shape3Milestones);
	shape3->setLeftBarBgGradient(shape3Milestones);

	shape3->setRangeDefaultMilestones(0.05, 0.2, 0.58);	

	curveEditorG->curveListComplete();

	curveEditorG2->curveListComplete();
	curveEditorG2->setTooltip(M("TP_COLORAPP_CURVEEDITOR2_TOOLTIP"));

	curveEditorG3->curveListComplete();
	curveEditorG3->setTooltip(M("TP_COLORAPP_CURVEEDITOR3_TOOLTIP"));
	tcmode3conn = toneCurveMode3->signal_changed().connect( sigc::mem_fun(*this, &ColorAppearance::curveMode3Changed), true );

	p2VBox->pack_start( *curveEditorG, Gtk::PACK_SHRINK, 2);
	p2VBox->pack_start( *curveEditorG2, Gtk::PACK_SHRINK, 2);
	p2VBox->pack_start( *curveEditorG3, Gtk::PACK_SHRINK, 2);
	
	// ------------------------ Choice CIECAM data


	datacie = Gtk::manage (new Gtk::CheckButton (M("TP_COLORAPP_DATACIE")));
	datacie->set_tooltip_markup (M("TP_COLORAPP_DATACIE_TOOLTIP"));
	datacieconn = datacie->signal_toggled().connect( sigc::mem_fun(*this, &ColorAppearance::datacie_toggled) );
	p2VBox->pack_start (*datacie);

    //-------------------------
	
	
	
	p2Frame->add(*p2VBox);
	
	
	pack_start (*p2Frame, Gtk::PACK_EXPAND_WIDGET, 4);

	

	// ------------------------ Process #3: Converting back to Lab/RGB


	// Process 3 frame
	Gtk::Frame *p3Frame;
	// Vertical box container for the content of the Process 3 frame
	Gtk::VBox *p3VBox;

	p3Frame = Gtk::manage (new Gtk::Frame(M("TP_COLORAPP_LABEL_VIEWING")) ); // "Editing viewing conditions" ???
	p3Frame->set_border_width(0);
	p3Frame->set_label_align(0.025, 0.5);

	p3VBox = Gtk::manage ( new Gtk::VBox());
	p3VBox->set_border_width(4);
	p3VBox->set_spacing(2);

	adaplum = Gtk::manage (new Adjuster (M("TP_COLORAPP_ADAPTVIEWING"), 0.1,  1000., 0.1,   16.));
	if (adaplum->delay < 1000) adaplum->delay = 1000;
	adaplum->throwOnButtonRelease();
	adaplum->set_tooltip_markup (M("TP_COLORAPP_ADAPTVIEWING_TOOLTIP"));
	p3VBox->pack_start (*adaplum);

	Gtk::HBox* surrHBox = Gtk::manage (new Gtk::HBox ());
	surrHBox->set_border_width (0);
	surrHBox->set_spacing (2);
	surrHBox->set_tooltip_markup(M("TP_COLORAPP_SURROUND_TOOLTIP"));
	Gtk::Label* surrLabel = Gtk::manage (new Gtk::Label (M("TP_COLORAPP_SURROUND")+":"));
	surrHBox->pack_start (*surrLabel, Gtk::PACK_SHRINK);
	surround = Gtk::manage (new MyComboBoxText ());
	surround->append_text (M("TP_COLORAPP_SURROUND_AVER"));
	surround->append_text (M("TP_COLORAPP_SURROUND_DIM"));
	surround->append_text (M("TP_COLORAPP_SURROUND_DARK"));
	surround->append_text (M("TP_COLORAPP_SURROUND_EXDARK"));
	surround->set_active (1);
	surrHBox->pack_start (*surround);
	p3VBox->pack_start (*surrHBox);

	p3Frame->add(*p3VBox);
	pack_start (*p3Frame, Gtk::PACK_EXPAND_WIDGET, 4);


	// ------------------------ Lab Gamut control


	gamut = Gtk::manage (new Gtk::CheckButton (M("TP_COLORAPP_GAMUT")));
	gamut->set_tooltip_markup (M("TP_COLORAPP_GAMUT_TOOLTIP"));
	gamutconn = gamut->signal_toggled().connect( sigc::mem_fun(*this, &ColorAppearance::gamut_toggled) );
	pack_start (*gamut, Gtk::PACK_SHRINK);
	
	// ------------------------ Bad pixel control

/*
	badpix = Gtk::manage (new Gtk::CheckButton (M("TP_COLORAPP_BADPIX")));
	badpix->set_tooltip_markup (M("TP_COLORAPP_BADPIX_TOOLTIP"));
	badpixconn = badpix->signal_toggled().connect( sigc::mem_fun(*this, &ColorAppearance::badpix_toggled) );
	pack_start (*badpix, Gtk::PACK_SHRINK);
*/	
	badpixsl = Gtk::manage (new Adjuster (M("TP_COLORAPP_BADPIXSL"), 0,  2, 1,  0));
	if (badpixsl->delay < 1000) badpixsl->delay = 1000;
	badpixsl->throwOnButtonRelease();
	badpixsl->set_tooltip_markup (M("TP_COLORAPP_BADPIXSL_TOOLTIP"));
	pack_start (*badpixsl, Gtk::PACK_SHRINK);
	
	// ------------------------ Listening events


	surrconn = surrsource->signal_toggled().connect( sigc::mem_fun(*this, &ColorAppearance::surrsource_toggled) );
	enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &ColorAppearance::enabledChanged) );
	wbmodelconn = wbmodel->signal_changed().connect ( sigc::mem_fun(*this, &ColorAppearance::wbmodelChanged) );
	algoconn = algo->signal_changed().connect ( sigc::mem_fun(*this, &ColorAppearance::algoChanged) );
	surroundconn = surround->signal_changed().connect ( sigc::mem_fun(*this, &ColorAppearance::surroundChanged) );

	degree->setAdjusterListener  (this);
	adapscen->setAdjusterListener (this);
	adaplum->setAdjusterListener (this);
	badpixsl->setAdjusterListener (this);
	jlight->setAdjusterListener  (this);
	qbright->setAdjusterListener  (this);
	colorh->setAdjusterListener  (this);
	chroma->setAdjusterListener  (this);
	schroma->setAdjusterListener  (this);
	mchroma->setAdjusterListener  (this);
	contrast->setAdjusterListener  (this);
	qcontrast->setAdjusterListener  (this);
	rstprotection->setAdjusterListener  (this);

	show_all();
}

ColorAppearance::~ColorAppearance () {
	delete curveEditorG;
	delete curveEditorG2;
	delete curveEditorG3;	
}



bool ColorAppearance::bgTTipQuery(int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip) {
	return true;
}

bool ColorAppearance::srTTipQuery(int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip) {
	return true;
}

void ColorAppearance::read (const ProcParams* pp, const ParamsEdited* pedited) {

	disableListener ();
	tcmodeconn.block(true);
	tcmode2conn.block(true);
	tcmode3conn.block(true);
	shape->setCurve (pp->colorappearance.curve);
	shape2->setCurve (pp->colorappearance.curve2);
	shape3->setCurve (pp->colorappearance.curve3);
	toneCurveMode->set_active(pp->colorappearance.curveMode);
	toneCurveMode2->set_active(pp->colorappearance.curveMode2);
	toneCurveMode3->set_active(pp->colorappearance.curveMode3);
	curveMode3Changed(); // This will set the correct sensitive state of depending Adjusters
	if (pedited) {
		degree->setEditedState        (pedited->colorappearance.degree ? Edited : UnEdited);
		adapscen->setEditedState      (pedited->colorappearance.adapscen ? Edited : UnEdited);
		adaplum->setEditedState       (pedited->colorappearance.adaplum ? Edited : UnEdited);
		badpixsl->setEditedState      (pedited->colorappearance.badpixsl ? Edited : UnEdited);
		jlight->setEditedState        (pedited->colorappearance.jlight ? Edited : UnEdited);
		qbright->setEditedState       (pedited->colorappearance.qbright ? Edited : UnEdited);
		chroma->setEditedState        (pedited->colorappearance.chroma ? Edited : UnEdited);
		schroma->setEditedState       (pedited->colorappearance.schroma ? Edited : UnEdited);
		mchroma->setEditedState       (pedited->colorappearance.mchroma ? Edited : UnEdited);
		rstprotection->setEditedState (pedited->colorappearance.rstprotection ? Edited : UnEdited);
		contrast->setEditedState      (pedited->colorappearance.contrast ? Edited : UnEdited);
		qcontrast->setEditedState     (pedited->colorappearance.qcontrast ? Edited : UnEdited);
		colorh->setEditedState        (pedited->colorappearance.colorh ? Edited : UnEdited);
		surrsource->set_inconsistent  (!pedited->colorappearance.surrsource);
		gamut->set_inconsistent       (!pedited->colorappearance.gamut);
	//	badpix->set_inconsistent      (!pedited->colorappearance.badpix);
		datacie->set_inconsistent     (!pedited->colorappearance.datacie);
		tonecie->set_inconsistent     (!pedited->colorappearance.tonecie);
	//	sharpcie->set_inconsistent    (!pedited->colorappearance.sharpcie);

		degree->setAutoInconsistent   (multiImage && !pedited->colorappearance.autodegree);
		adapscen->setAutoInconsistent   (multiImage && !pedited->colorappearance.autoadapscen);
		
		enabled->set_inconsistent     (multiImage && !pedited->colorappearance.enabled);
		shape->setUnChanged (!pedited->colorappearance.curve);
		shape2->setUnChanged (!pedited->colorappearance.curve2);
		shape3->setUnChanged (!pedited->colorappearance.curve3);
		if (!pedited->colorappearance.curveMode) {
			toneCurveMode->set_active(2);
		}
		if (!pedited->colorappearance.curveMode2) {
			toneCurveMode2->set_active(2);
		}
		if (!pedited->colorappearance.curveMode3) {
			toneCurveMode3->set_active(3);
		}


		}

	enaConn.block (true);
	enabled->set_active (pp->colorappearance.enabled);
	enaConn.block (false);

	surroundconn.block(true);
	if (pedited && !pedited->colorappearance.surround)
		surround->set_active (4);
	else if (pp->colorappearance.surround=="Average")
		surround->set_active (0);
	else if (pp->colorappearance.surround=="Dim")
		surround->set_active (1);
	else if (pp->colorappearance.surround=="Dark")
		surround->set_active (2);
	else if (pp->colorappearance.surround=="ExtremelyDark")
		surround->set_active (3);
	surroundconn.block(false);
	// Have to be manually called to handle initial state update
	surroundChanged();

	wbmodelconn.block(true);
	if (pedited && !pedited->colorappearance.wbmodel)
		wbmodel->set_active (2);
	else if (pp->colorappearance.wbmodel=="RawT")
		wbmodel->set_active (0);
	else if (pp->colorappearance.wbmodel=="RawTCAT02")
		wbmodel->set_active (1);
	wbmodelconn.block(false);
	// Have to be manually called to handle initial state update
	wbmodelChanged();

	algoconn.block(true);
	if (pedited && !pedited->colorappearance.algo)
		algo->set_active (4);
	else if (pp->colorappearance.algo=="JC")
		algo->set_active (0);
	else if (pp->colorappearance.algo=="JS")
		algo->set_active (1);
	else if (pp->colorappearance.algo=="QM")
		algo->set_active (2);
	else if (pp->colorappearance.algo=="ALL")
		algo->set_active (3);
	algoconn.block(false);
	// Have to be manually called to handle initial state update
	algoChanged();

	surrconn.block (true);
	surrsource->set_active (pp->colorappearance.surrsource);
	surrconn.block (false);
	gamutconn.block (true);
	gamut->set_active (pp->colorappearance.gamut);
	gamutconn.block (false);
//	badpixconn.block (true);
//	badpix->set_active (pp->colorappearance.badpix);
//	badpixconn.block (false);
	datacieconn.block (true);
	datacie->set_active (pp->colorappearance.datacie);
	datacieconn.block (false);
	tonecieconn.block (true);
	tonecie->set_active (pp->colorappearance.tonecie);
	tonecieconn.block (false);
//	sharpcieconn.block (true);
//	sharpcie->set_active (pp->colorappearance.sharpcie);
//	sharpcieconn.block (false);

	lastsurr=pp->colorappearance.surrsource;
	lastgamut=pp->colorappearance.gamut;
//	lastbadpix=pp->colorappearance.badpix;
	lastdatacie=pp->colorappearance.datacie;
	lasttonecie=pp->colorappearance.tonecie;
//	lastsharpcie=pp->colorappearance.sharpcie;

	lastEnabled = pp->colorappearance.enabled;
	lastAutoDegree = pp->colorappearance.autodegree;
	lastAutoAdapscen = pp->colorappearance.autoadapscen;

	degree->setValue (pp->colorappearance.degree);
	degree->setAutoValue(pp->colorappearance.autodegree);
	adapscen->setValue (pp->colorappearance.adapscen);
	adapscen->setAutoValue (pp->colorappearance.autoadapscen);
	
	adaplum->setValue (pp->colorappearance.adaplum);
	badpixsl->setValue (pp->colorappearance.badpixsl);
	jlight->setValue (pp->colorappearance.jlight);
	qbright->setValue (pp->colorappearance.qbright);
	chroma->setValue (pp->colorappearance.chroma);
	schroma->setValue (pp->colorappearance.schroma);
	mchroma->setValue (pp->colorappearance.mchroma);
	rstprotection->setValue (pp->colorappearance.rstprotection);
	contrast->setValue (pp->colorappearance.contrast);
	qcontrast->setValue (pp->colorappearance.qcontrast);
	colorh->setValue (pp->colorappearance.colorh);

	tcmode3conn.block(false);
	tcmode2conn.block(false);
	tcmodeconn.block(false);
	enableListener ();
}
void ColorAppearance::autoOpenCurve  () {
	shape->openIfNonlinear();
	shape2->openIfNonlinear();
	shape3->openIfNonlinear();

}


void ColorAppearance::write (ProcParams* pp, ParamsEdited* pedited) {

	pp->colorappearance.degree        = degree->getValue ();
	pp->colorappearance.autodegree    = degree->getAutoValue ();
	pp->colorappearance.enabled       = enabled->get_active();
	pp->colorappearance.adapscen      = adapscen->getValue ();
	pp->colorappearance.autoadapscen  = adapscen->getAutoValue ();
	pp->colorappearance.adaplum       = adaplum->getValue ();
	pp->colorappearance.badpixsl      = badpixsl->getValue ();
	pp->colorappearance.jlight        = jlight->getValue ();
	pp->colorappearance.qbright       = qbright->getValue ();
	pp->colorappearance.chroma        = chroma->getValue ();
	pp->colorappearance.schroma       = schroma->getValue ();
	pp->colorappearance.mchroma       = mchroma->getValue ();
	pp->colorappearance.contrast      = contrast->getValue ();
	pp->colorappearance.qcontrast     = qcontrast->getValue ();
	pp->colorappearance.colorh        = colorh->getValue ();
	pp->colorappearance.rstprotection = rstprotection->getValue ();
	pp->colorappearance.surrsource    = surrsource->get_active();
	pp->colorappearance.gamut         = gamut->get_active();
//	pp->colorappearance.badpix        = badpix->get_active();
	pp->colorappearance.datacie       = datacie->get_active();
	pp->colorappearance.tonecie       = tonecie->get_active();
//	pp->colorappearance.sharpcie      = sharpcie->get_active();
	pp->colorappearance.curve         = shape->getCurve ();
	pp->colorappearance.curve2        = shape2->getCurve ();
	pp->colorappearance.curve3        = shape3->getCurve ();
	
	int tcMode = toneCurveMode->get_active_row_number();
	if      (tcMode == 0) pp->colorappearance.curveMode = ColorAppearanceParams::TC_MODE_LIGHT;
	else if (tcMode == 1) pp->colorappearance.curveMode = ColorAppearanceParams::TC_MODE_BRIGHT;
	
	tcMode = toneCurveMode2->get_active_row_number();
	if      (tcMode == 0) pp->colorappearance.curveMode2 = ColorAppearanceParams::TC_MODE_LIGHT;
	else if (tcMode == 1) pp->colorappearance.curveMode2 = ColorAppearanceParams::TC_MODE_BRIGHT;
	
	int tcMode3 = toneCurveMode3->get_active_row_number();
	if      (tcMode3 == 0) pp->colorappearance.curveMode3 = ColorAppearanceParams::TC_MODE_CHROMA;
	else if (tcMode3 == 1) pp->colorappearance.curveMode3 = ColorAppearanceParams::TC_MODE_SATUR;
	else if (tcMode3 == 2) pp->colorappearance.curveMode3 = ColorAppearanceParams::TC_MODE_COLORF;

	if (pedited) {
		pedited->colorappearance.degree        = degree->getEditedState ();
		pedited->colorappearance.adapscen      = adapscen->getEditedState ();
		pedited->colorappearance.adaplum       = adaplum->getEditedState ();
		pedited->colorappearance.badpixsl      = badpixsl->getEditedState ();
		pedited->colorappearance.jlight        = jlight->getEditedState ();
		pedited->colorappearance.qbright       = qbright->getEditedState ();
		pedited->colorappearance.chroma        = chroma->getEditedState ();
		pedited->colorappearance.schroma       = schroma->getEditedState ();
		pedited->colorappearance.mchroma       = mchroma->getEditedState ();
		pedited->colorappearance.contrast      = contrast->getEditedState ();
		pedited->colorappearance.qcontrast     = qcontrast->getEditedState ();
		pedited->colorappearance.colorh        = colorh->getEditedState ();
		pedited->colorappearance.rstprotection = rstprotection->getEditedState ();
		pedited->colorappearance.autodegree    = !degree->getAutoInconsistent();
		pedited->colorappearance.autoadapscen  = !adapscen->getAutoInconsistent();
		pedited->colorappearance.enabled       = !enabled->get_inconsistent();
		pedited->colorappearance.surround      = surround->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->colorappearance.wbmodel       = wbmodel->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->colorappearance.algo          = algo->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->colorappearance.surrsource    = !surrsource->get_inconsistent();
		pedited->colorappearance.gamut         = !gamut->get_inconsistent();
	//	pedited->colorappearance.badpix        = !badpix->get_inconsistent();
		pedited->colorappearance.datacie       = !datacie->get_inconsistent();
		pedited->colorappearance.tonecie       = !tonecie->get_inconsistent();
	//	pedited->colorappearance.sharpcie      = !sharpcie->get_inconsistent();
		pedited->colorappearance.curve         = !shape->isUnChanged ();
		pedited->colorappearance.curve2        = !shape2->isUnChanged ();
		pedited->colorappearance.curve3        = !shape3->isUnChanged ();
		pedited->colorappearance.curveMode     = toneCurveMode->get_active_row_number() != 2;
		pedited->colorappearance.curveMode2    = toneCurveMode2->get_active_row_number() != 2;
		pedited->colorappearance.curveMode3    = toneCurveMode3->get_active_row_number() != 3;
	}
	if (surround->get_active_row_number()==0)
		pp->colorappearance.surround = "Average";
	else if (surround->get_active_row_number()==1)
		pp->colorappearance.surround = "Dim";
	else if (surround->get_active_row_number()==2)
		pp->colorappearance.surround = "Dark";
	else if (surround->get_active_row_number()==3)
		pp->colorappearance.surround = "ExtremelyDark";

	if (wbmodel->get_active_row_number()==0)
		pp->colorappearance.wbmodel = "RawT";
	else if (wbmodel->get_active_row_number()==1)
		pp->colorappearance.wbmodel = "RawTCAT02";

	if (algo->get_active_row_number()==0)
		pp->colorappearance.algo = "JC";
	else if (algo->get_active_row_number()==1)
		pp->colorappearance.algo = "JS";
	else if (algo->get_active_row_number()==2)
		pp->colorappearance.algo = "QM";
	else if (algo->get_active_row_number()==3)
		pp->colorappearance.algo = "ALL";

}
void ColorAppearance::curveChanged (CurveEditor* ce) {

	if (listener) {
		if (ce == shape)
			listener->panelChanged (EvCATCurve1, M("HISTORY_CUSTOMCURVE"));
		else if (ce == shape2)
			listener->panelChanged (EvCATCurve2, M("HISTORY_CUSTOMCURVE"));
		else if (ce == shape3)
			listener->panelChanged (EvCATCurve3, M("HISTORY_CUSTOMCURVE"));
	}
}

void ColorAppearance::curveMode1Changed () {
	if (listener)  Glib::signal_idle().connect (sigc::mem_fun(*this, &ColorAppearance::curveMode1Changed_));
}

bool ColorAppearance::curveMode1Changed_ () {
	if (listener) listener->panelChanged (EvCATCurveMode1, toneCurveMode->get_active_text());
	return false;
}

void ColorAppearance::curveMode2Changed () {
	if (listener)  Glib::signal_idle().connect (sigc::mem_fun(*this, &ColorAppearance::curveMode2Changed_));
}

bool ColorAppearance::curveMode2Changed_ () {
	if (listener) listener->panelChanged (EvCATCurveMode2, toneCurveMode2->get_active_text());
	return false;
}

void ColorAppearance::curveMode3Changed () {
	int tcMode3 = toneCurveMode3->get_active_row_number();
	if      (tcMode3 == 0) {chroma->set_sensitive(true);  schroma->set_sensitive(true); }
	else if (tcMode3 == 2) {chroma->set_sensitive(false); schroma->set_sensitive(false);}
	else if (tcMode3 == 1) {chroma->set_sensitive(false); schroma->set_sensitive(true); }

	if (listener)  Glib::signal_idle().connect (sigc::mem_fun(*this, &ColorAppearance::curveMode3Changed_));
}

bool ColorAppearance::curveMode3Changed_ () {
	if (listener) {
		listener->panelChanged (EvCATCurveMode3, toneCurveMode3->get_active_text());
	}
	return false;
}

void ColorAppearance::surrsource_toggled () {

	if (batchMode) {
		if (surrsource->get_inconsistent()) {
			surrsource->set_inconsistent (false);
			surrconn.block (true);
			surrsource->set_active (false);
			surrconn.block (false);
		}
		else if (lastsurr)
			surrsource->set_inconsistent (true);

		lastsurr = surrsource->get_active ();
	}

	if (listener) {
		if (surrsource->get_active ())
			listener->panelChanged (EvCATsurr, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvCATsurr, M("GENERAL_DISABLED"));
	}
}

void ColorAppearance::gamut_toggled () {

	if (batchMode) {
		if (gamut->get_inconsistent()) {
			gamut->set_inconsistent (false);
			gamutconn.block (true);
			gamut->set_active (false);
			gamutconn.block (false);
		}
		else if (lastgamut)
			gamut->set_inconsistent (true);

		lastgamut = gamut->get_active ();
	}
	if (listener) {
		if (gamut->get_active ())
			listener->panelChanged (EvCATgamut, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvCATgamut, M("GENERAL_DISABLED"));
	}
	

}
/*
void ColorAppearance::badpix_toggled () {

	if (batchMode) {
		if (badpix->get_inconsistent()) {
			badpix->set_inconsistent (false);
			badpixconn.block (true);
			badpix->set_active (false);
			badpixconn.block (false);
		}
		else if (lastbadpix)
			badpix->set_inconsistent (true);

		lastbadpix = badpix->get_active ();
	}
	if (listener) {
		if (badpix->get_active ())
			listener->panelChanged (EvCATbadpix, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvCATbadpix, M("GENERAL_DISABLED"));
	}
	

}
*/
void ColorAppearance::datacie_toggled () {

	if (batchMode) {
		if (datacie->get_inconsistent()) {
			datacie->set_inconsistent (false);
			datacieconn.block (true);
			datacie->set_active (false);
			datacieconn.block (false);
		}
		else if (lastdatacie)
			datacie->set_inconsistent (true);

		lastdatacie = datacie->get_active ();
	}

	if (listener) {
		if (datacie->get_active ())
			listener->panelChanged (EvCATdatacie, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvCATdatacie, M("GENERAL_DISABLED"));
	}
}
void ColorAppearance::tonecie_toggled () {

	if (batchMode) {
		if (tonecie->get_inconsistent()) {
			tonecie->set_inconsistent (false);
			tonecieconn.block (true);
			tonecie->set_active (false);
			tonecieconn.block (false);
		}
		else if (lasttonecie)
			tonecie->set_inconsistent (true);

		lasttonecie = tonecie->get_active ();
	}
	if (listener) {
		if (tonecie->get_active ())
			listener->panelChanged (EvCATtonecie, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvCATtonecie, M("GENERAL_DISABLED"));
	}

}
/*
void ColorAppearance::sharpcie_toggled () {

	if (batchMode) {
		if (sharpcie->get_inconsistent()) {
			sharpcie->set_inconsistent (false);
			sharpcieconn.block (true);
			sharpcie->set_active (false);
			sharpcieconn.block (false);
		}
		else if (lastsharpcie)
			sharpcie->set_inconsistent (true);

		lastsharpcie = sharpcie->get_active ();
	}
	if (listener) {
		if (sharpcie->get_active ())
			listener->panelChanged (EvCATsharpcie, M("GENERAL_ENABLED"));
		else
			listener->panelChanged (EvCATsharpcie, M("GENERAL_DISABLED"));
	}

}
*/

void ColorAppearance::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

	degree->setDefault (defParams->colorappearance.degree);
	adapscen->setDefault (defParams->colorappearance.adapscen);
	adaplum->setDefault (defParams->colorappearance.adaplum);
	badpixsl->setDefault (defParams->colorappearance.badpixsl);
	jlight->setDefault (defParams->colorappearance.jlight);
	qbright->setDefault (defParams->colorappearance.qbright);
	chroma->setDefault (defParams->colorappearance.chroma);
	schroma->setDefault (defParams->colorappearance.schroma);
	mchroma->setDefault (defParams->colorappearance.mchroma);
	rstprotection->setDefault (defParams->colorappearance.rstprotection);
	contrast->setDefault (defParams->colorappearance.contrast);
	qcontrast->setDefault (defParams->colorappearance.qcontrast);
	colorh->setDefault (defParams->colorappearance.colorh);

	if (pedited) {
		degree->setDefaultEditedState (pedited->colorappearance.degree ? Edited : UnEdited);
		adapscen->setDefaultEditedState (pedited->colorappearance.adapscen ? Edited : UnEdited);
		adaplum->setDefaultEditedState (pedited->colorappearance.adaplum ? Edited : UnEdited);
		badpixsl->setDefaultEditedState (pedited->colorappearance.badpixsl ? Edited : UnEdited);
		jlight->setDefaultEditedState (pedited->colorappearance.jlight ? Edited : UnEdited);
		qbright->setDefaultEditedState (pedited->colorappearance.qbright ? Edited : UnEdited);
		chroma->setDefaultEditedState (pedited->colorappearance.chroma ? Edited : UnEdited);
		schroma->setDefaultEditedState (pedited->colorappearance.schroma ? Edited : UnEdited);
		mchroma->setDefaultEditedState (pedited->colorappearance.mchroma ? Edited : UnEdited);
		rstprotection->setDefaultEditedState (pedited->colorappearance.rstprotection ? Edited : UnEdited);
		contrast->setDefaultEditedState (pedited->colorappearance.contrast ? Edited : UnEdited);
		qcontrast->setDefaultEditedState (pedited->colorappearance.qcontrast ? Edited : UnEdited);
		colorh->setDefaultEditedState (pedited->colorappearance.colorh ? Edited : UnEdited);
		
	}
	else {
		degree->setDefaultEditedState (Irrelevant);
		adapscen->setDefaultEditedState (Irrelevant);
		adaplum->setDefaultEditedState (Irrelevant);
		badpixsl->setDefaultEditedState (Irrelevant);
		jlight->setDefaultEditedState (Irrelevant);
		qbright->setDefaultEditedState (Irrelevant);
		chroma->setDefaultEditedState (Irrelevant);
		schroma->setDefaultEditedState (Irrelevant);
		mchroma->setDefaultEditedState (Irrelevant);
		contrast->setDefaultEditedState (Irrelevant);
		qcontrast->setDefaultEditedState (Irrelevant);
		rstprotection->setDefaultEditedState (Irrelevant);
		colorh->setDefaultEditedState (Irrelevant);
		
	}
}
int autoCamChangedUI (void* data) {
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
    (static_cast<ColorAppearance*>(data))->autoCamComputed_ ();
    return 0;
}
void ColorAppearance::autoCamChanged (double ccam) 
{  
    nextCcam = ccam;
    g_idle_add (autoCamChangedUI, this);
}

bool ColorAppearance::autoCamComputed_ () {

    disableListener ();
//	degree->setEnabled (true);
	degree->setValue (nextCcam);
    enableListener ();

    return false;
}
int adapCamChangedUI (void* data) {
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
    (static_cast<ColorAppearance*>(data))->adapCamComputed_ ();
    return 0;
}
void ColorAppearance::adapCamChanged (double cadap) 
{
    nextCadap = cadap;
    g_idle_add (adapCamChangedUI, this);
}

bool ColorAppearance::adapCamComputed_ () {

    disableListener ();
//	degree->setEnabled (true);
	adapscen->setValue (nextCadap);
    enableListener ();

    return false;
}


void ColorAppearance::colorForValue (double valX, double valY, int callerId, ColorCaller *caller) {

	float R, G, B;
	if (callerId == 1) {    // cc - bottom bar

		float value = (1.f - 0.7f) * float(valX) + 0.7f;
		// whole hue range
		// Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
		Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
	}
	caller->ccRed = double(R);
	caller->ccGreen = double(G);
	caller->ccBlue = double(B);
}

void ColorAppearance::adjusterChanged (Adjuster* a, double newval) {

	if (listener && (multiImage||enabled->get_active()) ) {
		if(a==degree)
			listener->panelChanged (EvCATDegree, a->getTextValue());
		else if(a==adapscen)
			listener->panelChanged (EvCATAdapscen, a->getTextValue());
		else if(a==adaplum)
			listener->panelChanged (EvCATAdapLum, a->getTextValue());
		else if(a==badpixsl)
			listener->panelChanged (EvCATbadpix, a->getTextValue());			
		else if(a==jlight)
			listener->panelChanged (EvCATJLight, a->getTextValue());
		else if(a==qbright)
			listener->panelChanged (EvCATQbright, a->getTextValue());
		else if(a==chroma)
			listener->panelChanged (EvCATChroma, a->getTextValue());
		else if(a==schroma)
			listener->panelChanged (EvCATSChroma, a->getTextValue());
		else if(a==mchroma)
			listener->panelChanged (EvCATMChroma, a->getTextValue());
		else if(a==rstprotection)
			listener->panelChanged (EvCATRstpro, a->getTextValue());
		else if(a==contrast)
			listener->panelChanged (EvCATContrast, a->getTextValue());
		else if(a==colorh)
			listener->panelChanged (EvCAThue, a->getTextValue());
		else if(a==qcontrast)
			listener->panelChanged (EvCATQContrast, a->getTextValue());
			
	}
}

void ColorAppearance::adjusterAutoToggled (Adjuster* a, bool newval) {

	if (multiImage) {
		if (degree->getAutoInconsistent()) {
			degree->setAutoInconsistent(false);
			degree->setAutoValue(false);
		}
		else if (lastAutoDegree)
			degree->setAutoInconsistent(true);

		lastAutoDegree = degree->getAutoValue();
		
		if (adapscen->getAutoInconsistent()) {
			adapscen->setAutoInconsistent(false);
			adapscen->setAutoValue(false);
		}
		else if (lastAutoAdapscen)
			adapscen->setAutoInconsistent(true);

		lastAutoAdapscen = adapscen->getAutoValue();
		
	}

	if (listener && (multiImage||enabled->get_active()) ) {
	
		if(a==degree) {
			if (degree->getAutoInconsistent())
				listener->panelChanged (EvCATAutoDegree, M("GENERAL_UNCHANGED"));
			else if (degree->getAutoValue())
				listener->panelChanged (EvCATAutoDegree, M("GENERAL_ENABLED"));
			else
				listener->panelChanged (EvCATAutoDegree, M("GENERAL_DISABLED"));
		}
		if(a==adapscen) {
			if (adapscen->getAutoInconsistent())
				listener->panelChanged (EvCATAutoAdap, M("GENERAL_UNCHANGED"));
			else if (adapscen->getAutoValue())
				listener->panelChanged (EvCATAutoAdap, M("GENERAL_ENABLED"));
			else
				listener->panelChanged (EvCATAutoAdap, M("GENERAL_DISABLED"));
		}
		
		
	}
}
void ColorAppearance::enabledChanged () {

	if (multiImage) {
		if (enabled->get_inconsistent()) {
			enabled->set_inconsistent (false);
			enaConn.block (true);
			enabled->set_active (false);
			enaConn.block (false);
		}
		else if (lastEnabled)
			enabled->set_inconsistent (true);

		lastEnabled = enabled->get_active ();
	}

	if (listener) {
		if (enabled->get_inconsistent())
			listener->panelChanged (EvCATEnabled, M("GENERAL_UNCHANGED"));
		else if (enabled->get_active ()) {
			listener->panelChanged (EvCATEnabled, M("GENERAL_ENABLED"));
			    curveEditorG->set_sensitive (true);
				toneCurveMode->set_sensitive (true);
			}
		else
			{listener->panelChanged (EvCATEnabled, M("GENERAL_DISABLED"));
			}
	}
}

void ColorAppearance::surroundChanged () {

	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvCATMethodsur, surround->get_active_text ());
	}
}

void ColorAppearance::wbmodelChanged () {

	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvCATMethodWB, wbmodel->get_active_text ());
	}
}


void ColorAppearance::algoChanged () {

	if ( algo->get_active_row_number()==0 ) {
		contrast->show();
		rstprotection->show();
		qcontrast->hide();
		jlight->show();
		mchroma->hide();
		chroma->show();
		schroma->hide();
		qbright->hide();
		colorh->hide();
		tonecie->hide();
	//	sharpcie->hide();
		curveEditorG->show();
		curveEditorG2->show();
		curveEditorG3->show();
	}
	else if ( algo->get_active_row_number()==1 ) {
		rstprotection->show();
		contrast->show();
		qcontrast->hide();
		jlight->show();
		mchroma->hide();
		chroma->hide();
		schroma->show();
		qbright->hide();
		colorh->hide();
		tonecie->hide();
//		sharpcie->hide();
		curveEditorG->show();
		curveEditorG2->show();
		curveEditorG3->show();
	}
	else if ( algo->get_active_row_number()==2 ) {
		contrast->hide();
		rstprotection->show();
		qcontrast->show();
		jlight->hide();
		mchroma->show();
		chroma->hide();
		schroma->hide();
		qbright->show();
		colorh->hide();
		tonecie->show();
	//	sharpcie->show();
	//	sharpcie->hide();
		curveEditorG->show();
		curveEditorG2->show();
		curveEditorG3->show();
	}
	else if ( algo->get_active_row_number()>=3 ) {  // ">=3" because everything has to be visible with the "(unchanged)" option too
		contrast->show();
		rstprotection->show();
		qcontrast->show();
		jlight->show();
		mchroma->show();
		chroma->show();
		schroma->show();
		qbright->show();
		colorh->show();
		tonecie->show();
//		sharpcie->show();
//		sharpcie->hide();
		curveEditorG->show();
		curveEditorG2->show();
		curveEditorG3->show();
	}

	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvCATMethodalg, algo->get_active_text ());
	}
}

void ColorAppearance::setBatchMode (bool batchMode) {

	ToolPanel::setBatchMode (batchMode);

	degree->showEditedCB ();
	adapscen->showEditedCB ();
	adaplum->showEditedCB ();
	badpixsl->showEditedCB ();
	jlight->showEditedCB ();
	qbright->showEditedCB ();
	chroma->showEditedCB ();
	schroma->showEditedCB ();
	mchroma->showEditedCB ();
	rstprotection->showEditedCB ();
	contrast->showEditedCB ();
	qcontrast->showEditedCB ();
	colorh->showEditedCB ();

	surround->append_text (M("GENERAL_UNCHANGED"));
	wbmodel->append_text (M("GENERAL_UNCHANGED"));
	algo->append_text (M("GENERAL_UNCHANGED"));
	toneCurveMode->append_text (M("GENERAL_UNCHANGED"));
	toneCurveMode2->append_text (M("GENERAL_UNCHANGED"));
	toneCurveMode3->append_text (M("GENERAL_UNCHANGED"));

	curveEditorG->setBatchMode (batchMode);
	curveEditorG2->setBatchMode (batchMode);
	curveEditorG3->setBatchMode (batchMode);
}

void ColorAppearance::updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve, LUTu & histCLurve, LUTu & histLLCurve, LUTu & histLCAM,  LUTu & histCCAM, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma){

	shape->updateBackgroundHistogram (histLCAM);
	shape3->updateBackgroundHistogram (histCCAM);
}



void ColorAppearance::setAdjusterBehavior (bool degreeadd, bool adapscenadd, bool adaplumadd, bool badpixsladd, bool jlightadd, bool chromaadd, bool contrastadd, bool rstprotectionadd, bool qbrightadd, bool qcontrastadd, bool schromaadd, bool mchromaadd, bool colorhadd) {

	degree->setAddMode(degreeadd);
	adapscen->setAddMode(adapscenadd);
	adaplum->setAddMode(adaplumadd);
	badpixsl->setAddMode(badpixsladd);
	jlight->setAddMode(jlightadd);
	qbright->setAddMode(qbrightadd);
	chroma->setAddMode(chromaadd);
	schroma->setAddMode(schromaadd);
	mchroma->setAddMode(mchromaadd);
	rstprotection->setAddMode(rstprotectionadd);
	contrast->setAddMode(contrastadd);
	qcontrast->setAddMode(qcontrastadd);
	colorh->setAddMode(colorhadd);
}

void ColorAppearance::trimValues (rtengine::procparams::ProcParams* pp) {

	degree->trimValue(pp->colorappearance.degree);
	adapscen->trimValue(pp->colorappearance.adapscen);
	adaplum->trimValue(pp->colorappearance.adaplum);
	badpixsl->trimValue(pp->colorappearance.badpixsl);
	jlight->trimValue(pp->colorappearance.jlight);
	qbright->trimValue(pp->colorappearance.qbright);
	chroma->trimValue(pp->colorappearance.chroma);
	schroma->trimValue(pp->colorappearance.schroma);
	mchroma->trimValue(pp->colorappearance.mchroma);
	rstprotection->trimValue(pp->colorappearance.rstprotection);
	contrast->trimValue(pp->colorappearance.contrast);
	qcontrast->trimValue(pp->colorappearance.qcontrast);
	colorh->trimValue(pp->colorappearance.colorh);
}
