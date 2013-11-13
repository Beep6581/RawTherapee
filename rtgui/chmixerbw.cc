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
#include "chmixerbw.h"
#include "rtimage.h"
#include "../rtengine/color.h"
#include <iomanip>
#include <cmath>
#include "guiutils.h"

using namespace rtengine;
using namespace rtengine::procparams;


ChMixerbw::ChMixerbw (): Gtk::VBox(), FoldableToolPanel(this) {
	CurveListener::setMulti(true);
    set_border_width(4);
	enabledLm = Gtk::manage (new Gtk::CheckButton (M("TP_BWMIX_ENABLED_LM")));	
	enabledLm->set_active (false);
	
	pack_start(*enabledLm, Gtk::PACK_SHRINK, 0);
	enabledLm->show ();

	Gtk::HBox* metHBox = Gtk::manage (new Gtk::HBox ());
	metHBox->set_border_width (0);
	metHBox->set_spacing (2);
	metHBox->set_tooltip_markup (M("TP_BWMIX_MET_TOOLTIP"));
	Gtk::Label* metLabel = Gtk::manage (new Gtk::Label (M("TP_BWMIX_MET")+":"));
	metHBox->pack_start (*metLabel, Gtk::PACK_SHRINK);
	met = Gtk::manage (new MyComboBoxText ());
	met->append_text (M("TP_BWMIX_MET0"));
	met->append_text (M("TP_BWMIX_MET1"));
	met->append_text (M("TP_BWMIX_MET2"));
	met->append_text (M("TP_BWMIX_MET3"));

	met->set_active (0);
	metHBox->pack_start (*met);
	pack_start (*metHBox);
	
    imgIcon[0] = Gtk::manage (new RTImage ("Chanmixer-R.png"));
    imgIcon[1] = Gtk::manage (new RTImage ("Chanmixer-G.png"));
    imgIcon[2] = Gtk::manage (new RTImage ("Chanmixer-B.png"));
    imgIcon[3] = Gtk::manage (new RTImage ("Chanmixer-O.png"));
    imgIcon[4] = Gtk::manage (new RTImage ("Chanmixer-Y.png"));
    imgIcon[5] = Gtk::manage (new RTImage ("Chanmixer-C.png"));
    imgIcon[6] = Gtk::manage (new RTImage ("Chanmixer-M.png"));
    imgIcon[7] = Gtk::manage (new RTImage ("Chanmixer-P.png"));
	
    imgIcon[8] = Gtk::manage (new RTImage ("Chanmixer-Rgamma.png"));
    imgIcon[9] = Gtk::manage (new RTImage ("Chanmixer-Ggamma.png"));
    imgIcon[10] = Gtk::manage (new RTImage ("Chanmixer-Bgamma.png"));

  std::vector<GradientMilestone> bottomMilestonesbw;
  bottomMilestonesbw.push_back( GradientMilestone(0., 0., 0., 0.) );
  bottomMilestonesbw.push_back( GradientMilestone(1., 1., 1., 1.) );
	
	
	std::vector<GradientMilestone> bottomMilestones;
	float R, G, B;
	// -0.1 rad < Hue < 1.6 rad
	for (int i=0; i<7; i++) {
		float x = float(i)*(1.0f/6.0);
		Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
		bottomMilestones.push_back( GradientMilestone(double(x), double(R), double(G), double(B)) );
	}
	curveEditorG = new CurveEditorGroup (options.lastBWCurvesDir, M("TP_BWMIX_CHANNEL"));
	curveEditorG->setCurveListener (this);
	vshape = static_cast<FlatCurveEditor*>(curveEditorG->addCurve(CT_Flat, M("TP_BWMIX_VAL")));
	vshape->setBottomBarBgGradient(bottomMilestones);
	vshape->setCurveColorProvider(this, 3);
	vshape->setTooltip(M("TP_BWMIX_CURVEEDITOR_LH_TOOLTIP"));
	
	curveEditorG->curveListComplete();

	pack_start (*curveEditorG, Gtk::PACK_SHRINK, 4);

		
	Gtk::HSeparator *hsep1 = Gtk::manage (new  Gtk::HSeparator());
	hsep1->show ();
	pack_start (*hsep1);

	enabled = Gtk::manage (new Gtk::CheckButton (M("TP_BWMIX_ENABLED")));
	
	enabled->set_active (true);
	enabled->set_tooltip_markup (M("TP_BWMIX_TOOLTIP"));
	
	pack_start(*enabled, Gtk::PACK_SHRINK, 0);
	enabled->show ();

	abox = Gtk::manage (new Gtk::HBox ());
	abox->set_border_width (2);

	autoch = Gtk::manage (new Gtk::ToggleButton (M("TP_BWMIX_AUTOCH")));
	autoch->set_tooltip_markup (M("TP_BWMIX_AUTOCH_TIP"));
	autoconn = autoch->signal_toggled().connect( sigc::mem_fun(*this, &ChMixerbw::autoch_toggled) );
	
	neutral = Gtk::manage (new Gtk::Button (M("TP_BWMIX_NEUTRAL")));
	neutral->set_tooltip_text (M("TP_BWMIX_NEUTRAL_TIP"));
	neutralconn = neutral->signal_pressed().connect( sigc::mem_fun(*this, &ChMixerbw::neutral_pressed) );
	neutral->show();
	
	abox->pack_start (*autoch);
	abox->pack_end (*neutral);
	abox->pack_end (*Gtk::manage (new Gtk::Label (" "))); //spacer
	pack_start (*abox);

	pack_start (*Gtk::manage (new  Gtk::HSeparator()));
	
	
	Gtk::HBox* setHBox = Gtk::manage (new Gtk::HBox ());
	setHBox->set_border_width (0);
	setHBox->set_spacing (2);
	setHBox->set_tooltip_markup (M("TP_BWMIX_SETTING_TOOLTIP"));
	setLabel = Gtk::manage (new Gtk::Label (M("TP_BWMIX_SETTING")+":"));
	
	setHBox->pack_start (*setLabel, Gtk::PACK_SHRINK);
	set = Gtk::manage (new MyComboBoxText ());
	set->append_text (M("TP_BWMIX_SET0"));
	set->append_text (M("TP_BWMIX_SET1"));
	set->append_text (M("TP_BWMIX_SET2"));
	set->append_text (M("TP_BWMIX_SET3"));
	set->append_text (M("TP_BWMIX_SET4"));
	set->append_text (M("TP_BWMIX_SET5"));
	set->append_text (M("TP_BWMIX_SET6"));
	set->append_text (M("TP_BWMIX_SET7"));
	set->append_text (M("TP_BWMIX_SET8"));
	set->append_text (M("TP_BWMIX_SET9"));
	set->append_text (M("TP_BWMIX_SET10"));
	set->append_text (M("TP_BWMIX_SET11"));
	set->append_text (M("TP_BWMIX_SET12"));
	set->append_text (M("TP_BWMIX_SET13"));
	set->append_text (M("TP_BWMIX_SET14"));

	set->set_active (0);
	setHBox->pack_start (*set);
	pack_start (*setHBox);
	
	Gtk::HBox* filHBox = Gtk::manage (new Gtk::HBox ());
	filHBox->set_border_width (0);
	filHBox->set_spacing (2);
	filHBox->set_tooltip_markup (M("TP_BWMIX_FILTER_TOOLTIP"));
	filLabel = Gtk::manage (new Gtk::Label (M("TP_BWMIX_FILTER")+":"));
	filHBox->pack_start (*filLabel, Gtk::PACK_SHRINK);
	fil = Gtk::manage (new MyComboBoxText ());
	fil->append_text (M("TP_BWMIX_FILTER0"));
	fil->append_text (M("TP_BWMIX_FILTER1"));
	fil->append_text (M("TP_BWMIX_FILTER2"));
	fil->append_text (M("TP_BWMIX_FILTER3"));
	fil->append_text (M("TP_BWMIX_FILTER4"));
	fil->append_text (M("TP_BWMIX_FILTER5"));
	fil->append_text (M("TP_BWMIX_FILTER6"));
	fil->append_text (M("TP_BWMIX_FILTER7"));
	fil->append_text (M("TP_BWMIX_FILTER8"));
	
	fil->set_active (0);
	filHBox->pack_start (*fil);
	pack_start (*filHBox);
	
	rlabel = Gtk::manage (new Gtk::Label ());
    rlabel->set_markup (Glib::ustring( M("TP_BWMIX_RED")));
    rlabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*rlabel, Gtk::PACK_SHRINK, 0);
	
	bwred= Gtk::manage(new Adjuster (imgIcon[0],-100,200,1,33));
	if (bwred->delay < 50) bwred->delay = 50;
	bwred->setAdjusterListener (this);
	bwred->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	
	bwred->show();
	pack_start( *bwred, Gtk::PACK_SHRINK, 0);

	orlabel = Gtk::manage (new Gtk::Label ());
    orlabel->set_markup (Glib::ustring( M("TP_BWMIX_ORA")));
    orlabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*orlabel, Gtk::PACK_SHRINK, 0);
	
	bworan= Gtk::manage(new Adjuster (imgIcon[3],-100,200,1,33));
	if (bworan->delay < 50) bworan->delay = 50;
	bworan->setAdjusterListener (this);
	bworan->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	bworan->show();
	pack_start( *bworan, Gtk::PACK_SHRINK, 0);

	ylabel = Gtk::manage (new Gtk::Label ());
    ylabel->set_markup (Glib::ustring( M("TP_BWMIX_YEL")));
    ylabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*ylabel, Gtk::PACK_SHRINK, 0);
	
	bwyell= Gtk::manage(new Adjuster (imgIcon[4],-100,200,1,33));
	if (bwyell->delay < 50) bwyell->delay = 50;
	bwyell->setAdjusterListener (this);
	bwyell->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	bwyell->show();
	pack_start( *bwyell, Gtk::PACK_SHRINK, 0);

	
	glabel = Gtk::manage (new Gtk::Label ());
    glabel->set_markup (Glib::ustring( M("TP_BWMIX_GREEN")));
    glabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*glabel, Gtk::PACK_SHRINK, 0);
	
	bwgreen= Gtk::manage(new Adjuster (imgIcon[1],-100,200,1,33));
	if (bwgreen->delay < 50) bwgreen->delay = 50;
	
	bwgreen->setAdjusterListener (this);
	bwgreen->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	
	bwgreen->show();
	pack_start( *bwgreen, Gtk::PACK_SHRINK, 0);

	clabel = Gtk::manage (new Gtk::Label ());
    clabel->set_markup (Glib::ustring( M("TP_BWMIX_CYAN")));
    clabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*clabel, Gtk::PACK_SHRINK, 0);
	
	bwcyan= Gtk::manage(new Adjuster (imgIcon[5],-100,200,1,33));
	if (bwcyan->delay < 50) bwcyan->delay = 50;
	bwcyan->setAdjusterListener (this);
	bwcyan->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	bwcyan->show();
	pack_start( *bwcyan, Gtk::PACK_SHRINK, 0);
	
	blabel = Gtk::manage (new Gtk::Label ());
    blabel->set_markup (Glib::ustring( M("TP_BWMIX_BLUE")));
    blabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*blabel, Gtk::PACK_SHRINK, 0);

	bwblue= Gtk::manage(new Adjuster (imgIcon[2],-100,200,1,33));
	if (bwblue->delay < 50) bwblue->delay = 50;
	
	bwblue->setAdjusterListener (this);
	bwblue->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	
	bwblue->show();
	pack_start( *bwblue, Gtk::PACK_SHRINK, 0);

	mlabel = Gtk::manage (new Gtk::Label ());
    mlabel->set_markup (Glib::ustring( M("TP_BWMIX_MAG")));
    mlabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*mlabel, Gtk::PACK_SHRINK, 0);
	
	bwmag= Gtk::manage(new Adjuster (imgIcon[6],-100,200,1,33));
	if (bwmag->delay < 50) bwmag->delay = 50;
	bwmag->setAdjusterListener (this);
	bwmag->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	bwmag->show();
	pack_start( *bwmag, Gtk::PACK_SHRINK, 0);

	plabel = Gtk::manage (new Gtk::Label ());
    plabel->set_markup (Glib::ustring( M("TP_BWMIX_PUR")));
    plabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*plabel, Gtk::PACK_SHRINK, 0);
	
	bwpur= Gtk::manage(new Adjuster (imgIcon[7],-100,200,1,33));
	if (bwpur->delay < 50) bwpur->delay = 50;
	bwpur->setAdjusterListener (this);
	bwpur->set_tooltip_markup (M("TP_BWMIX_RGB_TOOLTIP"));
	bwpur->show();
	pack_start( *bwpur, Gtk::PACK_SHRINK, 0);

    Gamlabel = Gtk::manage (new Gtk::Label ());
    Gamlabel->set_markup (Glib::ustring("<b>") + M("TP_BWMIX_GAMMA") + Glib::ustring("</b>"));
    Gamlabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*Gamlabel, Gtk::PACK_SHRINK, 0);

	
	rglabel = Gtk::manage (new Gtk::Label ());
    rglabel->set_markup (Glib::ustring( M("TP_BWMIX_GAM_RED")));
    rglabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*rglabel, Gtk::PACK_SHRINK, 0);
	
	bwredgam= Gtk::manage(new Adjuster (imgIcon[8],-100,100,1,0));
	if (bwredgam->delay < 50) bwredgam->delay = 50;
	
	bwredgam->setAdjusterListener (this);
	bwredgam->set_tooltip_markup (M("TP_BWMIX_GAM_TOOLTIP"));
	bwredgam->show();
	pack_start( *bwredgam, Gtk::PACK_SHRINK, 0);
	
	gglabel = Gtk::manage (new Gtk::Label ());
    gglabel->set_markup (Glib::ustring( M("TP_BWMIX_GAM_GREEN")));
    gglabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*gglabel, Gtk::PACK_SHRINK, 0);

	bwgreengam= Gtk::manage(new Adjuster (imgIcon[9],-100,100,1,0));
	if (bwgreengam->delay < 50) bwgreengam->delay = 50;
	bwgreengam->setAdjusterListener (this);
	bwgreengam->set_tooltip_markup (M("TP_BWMIX_GAM_TOOLTIP"));	
	bwgreengam->show();
	pack_start( *bwgreengam, Gtk::PACK_SHRINK, 0);

	bglabel = Gtk::manage (new Gtk::Label ());
    bglabel->set_markup (Glib::ustring( M("TP_BWMIX_GAM_BLUE")));
    bglabel->set_alignment(Gtk::ALIGN_LEFT);
    pack_start (*bglabel, Gtk::PACK_SHRINK, 0);
	
	bwbluegam= Gtk::manage(new Adjuster (imgIcon[10],-100,100,1,0));
	if (bwbluegam->delay < 50) bwbluegam->delay = 50;
	bwbluegam->setAdjusterListener (this);
	bwbluegam->set_tooltip_markup (M("TP_BWMIX_GAM_TOOLTIP"));
	bwbluegam->show();
	pack_start( *bwbluegam, Gtk::PACK_SHRINK, 0);

	enaLmconn    = enabledLm->signal_toggled().connect( sigc::mem_fun(*this, &ChMixerbw::enabledLm_toggled) );
	
	enaconn    = enabled->signal_toggled().connect( sigc::mem_fun(*this, &ChMixerbw::enabled_toggled) );
	filconn = fil->signal_changed().connect ( sigc::mem_fun(*this, &ChMixerbw::filChanged) );
	setconn = set->signal_changed().connect ( sigc::mem_fun(*this, &ChMixerbw::setChanged) );
	metconn = met->signal_changed().connect ( sigc::mem_fun(*this, &ChMixerbw::metChanged) );

//----------- Curve 1 ------------------------------
  pack_start (*Gtk::manage (new  Gtk::HSeparator()));

  toneCurveBW = Gtk::manage (new MyComboBoxText ());
  toneCurveBW->append_text (M("TP_BWMIX_TCMODE_STANDARD"));
  toneCurveBW->append_text (M("TP_BWMIX_TCMODE_WEIGHTEDSTD"));
  toneCurveBW->append_text (M("TP_BWMIX_TCMODE_FILMLIKE"));
  toneCurveBW->append_text (M("TP_BWMIX_TCMODE_SATANDVALBLENDING"));
  toneCurveBW->set_active (0);
  toneCurveBW->set_tooltip_text(M("TP_BWMIX_TCMODE_LABEL1"));

  curveEditorGBW = new CurveEditorGroup (options.lastBWCurvesDir, M("TP_BWMIX_CURVEEDITOR1"));
  curveEditorGBW->setCurveListener (this);

  shape = static_cast<DiagonalCurveEditor*>(curveEditorGBW->addCurve(CT_Diagonal, "", toneCurveBW));
  shape->setBottomBarBgGradient(bottomMilestonesbw);
  shape->setLeftBarBgGradient(bottomMilestonesbw);
  shape->setTooltip(M("TP_BWMIX_CURVEEDITOR_BEFORE_TOOLTIP"));

  // This will add the reset button at the end of the curveType buttons
  curveEditorGBW->curveListComplete();

  pack_start( *curveEditorGBW, Gtk::PACK_SHRINK, 2);

  tcmodeconn = toneCurveBW->signal_changed().connect( sigc::mem_fun(*this, &ChMixerbw::curveMode1Changed), true );

//----------- Curve 2 ------------------------------
  toneCurveBW2 = Gtk::manage (new MyComboBoxText ());
  toneCurveBW2->append_text (M("TP_BWMIX_TCMODE_STANDARD"));
//  toneCurveBW2->append_text (M("TP_BWMIX_TCMODE_WEIGHTEDSTD"));
  toneCurveBW2->set_active (0);
  toneCurveBW2->set_tooltip_text(M("TP_BWMIX_TCMODE_LABEL2"));

  curveEditorGBW2 = new CurveEditorGroup (options.lastBWCurvesDir, M("TP_BWMIX_CURVEEDITOR2"));
  curveEditorGBW2->setCurveListener (this);

  shape2 = static_cast<DiagonalCurveEditor*>(curveEditorGBW2->addCurve(CT_Diagonal, "", toneCurveBW2));
  shape2->setBottomBarBgGradient(bottomMilestonesbw);
  shape2->setLeftBarBgGradient(bottomMilestonesbw);
  shape2->setTooltip(M("TP_BWMIX_CURVEEDITOR_AFTER_TOOLTIP"));

  curveEditorGBW2->curveListComplete();

  pack_start( *curveEditorGBW2, Gtk::PACK_SHRINK, 2);

  tcmodeconn2 = toneCurveBW2->signal_changed().connect( sigc::mem_fun(*this, &ChMixerbw::curveMode1Changed2), true );
		
  show_all();
}
ChMixerbw::~ChMixerbw () {
	delete curveEditorG;
	delete curveEditorGBW;
	delete curveEditorGBW2;
}

int BWChangedUI (void* data) {
    GThreadLock lock; 
    (static_cast<ChMixerbw*>(data))->BWComputed_ ();
    return 0;
}

void ChMixerbw::BWChanged  (double redbw, double greenbw, double bluebw){
    nextredbw = redbw;
    nextgreenbw = greenbw;
    nextbluebw = bluebw;
    g_idle_add (BWChangedUI, this);
}

bool ChMixerbw::BWComputed_ () {

    disableListener ();
    bwred->setValue (nextredbw);
    bwgreen->setValue (nextgreenbw);
    bwblue->setValue (nextbluebw);
    enableListener ();

    return false;
}

void ChMixerbw::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();
	metconn.block(true);
    autoconn.block (true);
    autoch->set_active (pp->chmixerbw.autoc);
    lastAuto = pp->chmixerbw.autoc;
	
	if (pedited && !pedited->chmixerbw.met)
		met->set_active (4);
	else if (pp->chmixerbw.met=="No")
		met->set_active (0);
	else if (pp->chmixerbw.met=="De")
		met->set_active (1);
	else if (pp->chmixerbw.met=="Le")
		met->set_active (2);
	else if (pp->chmixerbw.met=="Ch")
		met->set_active (3);
	metconn.block(false);
	metChanged();
	
	
	filconn.block(true);
	if (pedited && !pedited->chmixerbw.fil)
		fil->set_active (9);
	else if (pp->chmixerbw.fil=="No")
		fil->set_active (0);
	else if (pp->chmixerbw.fil=="Re")
		fil->set_active (1);
	else if (pp->chmixerbw.fil=="Or")
		fil->set_active (2);
	else if (pp->chmixerbw.fil=="Ye")
		fil->set_active (3);
	else if (pp->chmixerbw.fil=="Yg")
		fil->set_active (4);
	else if (pp->chmixerbw.fil=="Gr")
		fil->set_active (5);
	else if (pp->chmixerbw.fil=="Cy")
		fil->set_active (6);
	else if (pp->chmixerbw.fil=="Bl")
		fil->set_active (7);
	else if (pp->chmixerbw.fil=="Pu")
		fil->set_active (8);	
	filconn.block(false);
	filChanged();

	setconn.block(true);
	if (pedited && !pedited->chmixerbw.set)
		set->set_active (15);
	else if (pp->chmixerbw.set=="Nc")
		set->set_active (0);
	else if (pp->chmixerbw.set=="Hc")
		set->set_active (1);
	else if (pp->chmixerbw.set=="Lu")
		set->set_active (2);
	else if (pp->chmixerbw.set=="La")
		set->set_active (3);
	else if (pp->chmixerbw.set=="Po")
		set->set_active (4);
	else if (pp->chmixerbw.set=="Ls")
		set->set_active (5);
	else if (pp->chmixerbw.set=="Hs")
		set->set_active (6);		
	else if (pp->chmixerbw.set=="Pa")
		set->set_active (7);		
	else if (pp->chmixerbw.set=="Hp")
		set->set_active (8);
	else if (pp->chmixerbw.set=="Or")
		set->set_active (9);		
	else if (pp->chmixerbw.set=="Ma")
		set->set_active (10);
	else if (pp->chmixerbw.set=="Mr")
		set->set_active (11);
	else if (pp->chmixerbw.set=="Fa")
		set->set_active (12);
	else if (pp->chmixerbw.set=="Fr")
		set->set_active (13);
	else if (pp->chmixerbw.set=="Ir")
		set->set_active (14);

		
	setconn.block(false);
	setChanged();
	

    if (pedited) {
        vshape->setUnChanged (!pedited->chmixerbw.vcurve);
        shape->setUnChanged (!pedited->chmixerbw.curve);
        shape2->setUnChanged (!pedited->chmixerbw.curve2);
        autoch->set_inconsistent (!pedited->chmixerbw.autoc);
		enabledLm->set_inconsistent  (!pedited->chmixerbw.enabledLm);
		enabled->set_inconsistent  (!pedited->chmixerbw.enabled);
		bwred->setEditedState (pedited->chmixerbw.bwred ? Edited : UnEdited);
		bwgreen->setEditedState (pedited->chmixerbw.bwgreen ? Edited : UnEdited);
		bwblue->setEditedState (pedited->chmixerbw.bwblue ? Edited : UnEdited);
		bwredgam->setEditedState (pedited->chmixerbw.bwredgam ? Edited : UnEdited);
		bwgreengam->setEditedState (pedited->chmixerbw.bwgreengam ? Edited : UnEdited);
		bwbluegam->setEditedState (pedited->chmixerbw.bwbluegam ? Edited : UnEdited);
		bworan->setEditedState (pedited->chmixerbw.bworan ? Edited : UnEdited);
		bwyell->setEditedState (pedited->chmixerbw.bwyell ? Edited : UnEdited);
		bwcyan->setEditedState (pedited->chmixerbw.bwcyan ? Edited : UnEdited);
		bwmag->setEditedState (pedited->chmixerbw.bwmag ? Edited : UnEdited);
		bwpur->setEditedState (pedited->chmixerbw.bwpur ? Edited : UnEdited);
        if (!pedited->chmixerbw.curveMode) {
            toneCurveBW->set_active(4);
        }
        if (!pedited->chmixerbw.curveMode2) {
            toneCurveBW2->set_active(1);
        }	
	}
    autoconn.block (false);	
	enaconn.block (true);
	enabled->set_active (pp->chmixerbw.enabled);
	enaconn.block (false);
	lastEnabled = pp->chmixerbw.enabled;
	enaLmconn.block (true);
	enabledLm->set_active (pp->chmixerbw.enabledLm);
	enaLmconn.block (false);
	lastEnabledLm = pp->chmixerbw.enabledLm;
	bwred->setValue (pp->chmixerbw.bwred);
	bwgreen->setValue (pp->chmixerbw.bwgreen);
	bwblue->setValue (pp->chmixerbw.bwblue);
	bwredgam->setValue (pp->chmixerbw.bwredgam);
	bwgreengam->setValue (pp->chmixerbw.bwgreengam);
	bwbluegam->setValue (pp->chmixerbw.bwbluegam);
	bworan->setValue (pp->chmixerbw.bworan);
	bwyell->setValue (pp->chmixerbw.bwyell);
	bwcyan->setValue (pp->chmixerbw.bwcyan);
	bwmag->setValue (pp->chmixerbw.bwmag);
	bwpur->setValue (pp->chmixerbw.bwpur);
    vshape->setCurve         (pp->chmixerbw.vcurve);
    shape->setCurve         (pp->chmixerbw.curve);
    toneCurveBW->set_active(pp->chmixerbw.curveMode);
    shape2->setCurve         (pp->chmixerbw.curve2);
    toneCurveBW2->set_active(pp->chmixerbw.curveMode2);
    enableListener ();
}

void ChMixerbw::write (ProcParams* pp, ParamsEdited* pedited) {
	pp->chmixerbw.enabled    = enabled->get_active ();
	pp->chmixerbw.enabledLm    = enabledLm->get_active ();
    pp->chmixerbw.autoc = autoch->get_active();
	pp->chmixerbw.bwred = bwred->getValue ();
	pp->chmixerbw.bwgreen = bwgreen->getValue ();
	pp->chmixerbw.bwblue = bwblue->getValue ();
	pp->chmixerbw.bwredgam = bwredgam->getValue ();
	pp->chmixerbw.bwgreengam = bwgreengam->getValue ();
	pp->chmixerbw.bwbluegam = bwbluegam->getValue ();
	pp->chmixerbw.bworan = bworan->getValue ();
	pp->chmixerbw.bwyell = bwyell->getValue ();
	pp->chmixerbw.bwcyan = bwcyan->getValue ();
	pp->chmixerbw.bwmag = bwmag->getValue ();
	pp->chmixerbw.bwpur = bwpur->getValue ();
    pp->chmixerbw.vcurve = vshape->getCurve ();
    pp->chmixerbw.curve = shape->getCurve ();
    pp->chmixerbw.curve2 = shape2->getCurve ();
 
    int tcMode = toneCurveBW->get_active_row_number();
    if      (tcMode == 0) pp->chmixerbw.curveMode = ChannelMixerbwParams::TC_MODE_STD_BW;
    else if (tcMode == 1) pp->chmixerbw.curveMode = ChannelMixerbwParams::TC_MODE_WEIGHTEDSTD_BW;
    else if (tcMode == 2) pp->chmixerbw.curveMode = ChannelMixerbwParams::TC_MODE_FILMLIKE_BW;
    else if (tcMode == 3) pp->chmixerbw.curveMode = ChannelMixerbwParams::TC_MODE_SATANDVALBLENDING_BW;

    tcMode = toneCurveBW2->get_active_row_number();
    if      (tcMode == 0) pp->chmixerbw.curveMode2 = ChannelMixerbwParams::TC_MODE_STD_BW;
  //  else if (tcMode == 1) pp->chmixerbw.curveMode2 = ChannelMixerbwParams::TC_MODE_WEIGHTEDSTD;
	
    if (pedited) { 
        pedited->chmixerbw.vcurve = !vshape->isUnChanged ();
		pedited->chmixerbw.enabledLm    = !enabledLm->get_inconsistent();
        pedited->chmixerbw.curve = !shape->isUnChanged ();
        pedited->chmixerbw.autoc    = !autoch->get_inconsistent();
		pedited->chmixerbw.enabled    = !enabled->get_inconsistent();
		pedited->chmixerbw.bwred = bwred->getEditedState ();
		pedited->chmixerbw.bwgreen = bwgreen->getEditedState ();
		pedited->chmixerbw.bwblue = bwblue->getEditedState ();
		pedited->chmixerbw.bwredgam = bwredgam->getEditedState ();
		pedited->chmixerbw.bwgreengam = bwgreengam->getEditedState ();
		pedited->chmixerbw.bwbluegam = bwbluegam->getEditedState ();
		pedited->chmixerbw.fil   = fil->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->chmixerbw.set   = set->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->chmixerbw.met   = met->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->chmixerbw.bworan = bworan->getEditedState ();
		pedited->chmixerbw.bwyell = bwyell->getEditedState ();
		pedited->chmixerbw.bwcyan = bwcyan->getEditedState ();
		pedited->chmixerbw.bwmag = bwmag->getEditedState ();
		pedited->chmixerbw.bwpur = bwpur->getEditedState ();
        pedited->chmixerbw.curveMode  = toneCurveBW->get_active_row_number() != 4;
        pedited->chmixerbw.curveMode2  = toneCurveBW2->get_active_row_number() != 1;
	}
	if (met->get_active_row_number()==0)
		pp->chmixerbw.met = "No";
	else if (met->get_active_row_number()==1)
		pp->chmixerbw.met = "De";
	else if (met->get_active_row_number()==2)
		pp->chmixerbw.met = "Le";
	else if (met->get_active_row_number()==3)
		pp->chmixerbw.met = "Ch";
	
	
	if (set->get_active_row_number()==0)
		pp->chmixerbw.set = "Nc";
	else if (set->get_active_row_number()==1)
		pp->chmixerbw.set = "Hc";
	else if (set->get_active_row_number()==2)
		pp->chmixerbw.set = "Lu";
	else if (set->get_active_row_number()==3)
		pp->chmixerbw.set = "La";
	else if (set->get_active_row_number()==4)
		pp->chmixerbw.set = "Po";
	else if (set->get_active_row_number()==5)
		pp->chmixerbw.set = "Ls";
	else if (set->get_active_row_number()==6)
		pp->chmixerbw.set = "Hs";		
	else if (set->get_active_row_number()==7)
		pp->chmixerbw.set = "Pa";		
	else if (set->get_active_row_number()==8)
		pp->chmixerbw.set = "Hp";	
	else if (set->get_active_row_number()==9)
		pp->chmixerbw.set = "Or";		
	else if (set->get_active_row_number()==10)
		pp->chmixerbw.set = "Ma";
	else if (set->get_active_row_number()==11)
		pp->chmixerbw.set = "Mr";
	else if (set->get_active_row_number()==12)
		pp->chmixerbw.set = "Fa";
	else if (set->get_active_row_number()==13)
		pp->chmixerbw.set = "Fr";
	else if (set->get_active_row_number()==14)
		pp->chmixerbw.set = "Ir";
		
	
	if (fil->get_active_row_number()==0)
		pp->chmixerbw.fil = "No";
	else if (fil->get_active_row_number()==1)
		pp->chmixerbw.fil = "Re";
	else if (fil->get_active_row_number()==2)
		pp->chmixerbw.fil = "Or";
	else if (fil->get_active_row_number()==3)
		pp->chmixerbw.fil = "Ye";
	else if (fil->get_active_row_number()==4)
		pp->chmixerbw.fil = "Yg";
	else if (fil->get_active_row_number()==5)
		pp->chmixerbw.fil = "Gr";
	else if (fil->get_active_row_number()==6)
		pp->chmixerbw.fil = "Cy";	
	else if (fil->get_active_row_number()==7)
		pp->chmixerbw.fil = "Bl";
	else if (fil->get_active_row_number()==8)
		pp->chmixerbw.fil = "Pu";	
}
void ChMixerbw::curveChanged (CurveEditor* ce) {
	if (listener) {
        if (ce == shape)
            listener->panelChanged (EvToneCurvebw1, M("HISTORY_CUSTOMCURVE"));
        if (ce == shape2)
            listener->panelChanged (EvToneCurvebw2, M("HISTORY_CUSTOMCURVE"));
		if (ce == vshape)
			listener->panelChanged (EvBWequalV, M("HISTORY_CUSTOMCURVE"));	
	}
}

void ChMixerbw::curveMode1Changed () {
    if (listener)  Glib::signal_idle().connect (sigc::mem_fun(*this, &ChMixerbw::curveMode1Changed_));
}
bool ChMixerbw::curveMode1Changed_ () {
    if (listener) listener->panelChanged (EvToneCurveBWMode1, toneCurveBW->get_active_text());
    return false;
}
void ChMixerbw::curveMode1Changed2 () {
    if (listener)  Glib::signal_idle().connect (sigc::mem_fun(*this, &ChMixerbw::curveMode1Changed2_));
}
bool ChMixerbw::curveMode1Changed2_ () {
    if (listener) listener->panelChanged (EvToneCurveBWMode2, toneCurveBW2->get_active_text());
    return false;
}

void ChMixerbw::colorForValue (double valX, double valY, int callerId, ColorCaller* caller) {

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


void ChMixerbw::setChanged () {

	if ( set->get_active_row_number()==10 || set->get_active_row_number()==11 ) {
	bwred->show();
	bwred->set_sensitive (true);
	bwgreen->show();
	bwgreen->set_sensitive (true);
	bwblue->show();
	bwblue->set_sensitive (true);
	bwredgam->show();
	bwgreengam->show();
	bwbluegam->show();
	bworan->hide();
	bwyell->hide();
	bwcyan->hide();
	bwmag->hide();
	bwpur->hide();
	orlabel->hide();
	ylabel->hide();
	clabel->hide();
	mlabel->hide();
	plabel->hide();	
	Gamlabel->show();
	rlabel->show();
	glabel->show();
	blabel->show();
	rglabel->show();
	gglabel->show();
	bglabel->show();
	enabled->hide();
	fil->set_sensitive (true);
	
	}
	else if ( set->get_active_row_number()==12 || set->get_active_row_number()==13 ) {
	bwred->show();
	bwred->set_sensitive (true);
	bwgreen->show();
	bwgreen->set_sensitive (true);
	bwblue->show();
	bwblue->set_sensitive (true);
	bwredgam->show();
	bwgreengam->show();
	bwbluegam->show();
	bworan->show();
	bwyell->show();
	bwcyan->show();
	bwmag->show();
	bwpur->show();
	orlabel->show();
	ylabel->show();
	clabel->show();
	mlabel->show();
	plabel->show();
	Gamlabel->show();
	rlabel->show();
	glabel->show();
	blabel->show();
	rglabel->show();
	gglabel->show();
	bglabel->show();
	enabled->show();
	fil->set_sensitive (true);
	
	}
	else if ( set->get_active_row_number()==14 ) {
 	fil->set_active (0);	
	fil->set_sensitive (false);
	}
	else  {
	bwred->show();
	bwred->set_sensitive (false);
	bwgreen->show();
	bwgreen->set_sensitive (false);
	bwblue->show();
	bwblue->set_sensitive (false);
	rlabel->show();
	glabel->show();
	blabel->show();
	rglabel->show();
	gglabel->show();
	bglabel->show();
	Gamlabel->show();
	bwredgam->show();
	bwgreengam->show();
	bwbluegam->show();	
	bworan->hide();
	bwyell->hide();
	bwcyan->hide();
	bwmag->hide();
	bwpur->hide();
	orlabel->hide();
	ylabel->hide();
	clabel->hide();
	mlabel->hide();
	plabel->hide();	
	enabled->hide();
	fil->set_sensitive (true);
	
	}

	if (listener && (multiImage||enabledLm->get_active())) {
		listener->panelChanged (EvBWset, set->get_active_text ());
	}
}


void ChMixerbw::filChanged () {
	if (listener && (multiImage||enabledLm->get_active())) {	
		listener->panelChanged (EvBWfil, fil->get_active_text ());
	}
}

void ChMixerbw::metChanged () {
	 if(met->get_active_row_number()==3) {
			set->show();
			setLabel->show();
			//enabled->show();
			enabledLm->show();
			curveEditorG->hide();
			curveEditorGBW->show();
			curveEditorGBW2->show();
			autoch->show();
			neutral->show();
			fil->show();
			filLabel->show();

		if(set->get_active_row_number()==10 || set->get_active_row_number()==11 || set->get_active_row_number()==12 || set->get_active_row_number()==13){
			bwred->show();
			bwgreen->show();
			bwblue->show();
			rlabel->show();
			glabel->show();
			blabel->show();
			rglabel->show();
			gglabel->show();
			bglabel->show();
			bwredgam->show();
			bwgreengam->show();
			bwbluegam->show();
			Gamlabel->show();	
			enabled->hide();		
		}
		if(set->get_active_row_number()==12 || set->get_active_row_number()==13) {
			bworan->show();
			bwyell->show();
			bwcyan->show();
			bwmag->show();
			bwpur->show();
			orlabel->show();
			ylabel->show();
			clabel->show();
			mlabel->show();
			plabel->show();
			enabled->show();
		}
	}
		
	 else if(met->get_active_row_number()==2) {
		curveEditorG->show();
		curveEditorGBW->show();
		curveEditorGBW2->show();
		autoch->hide();
		neutral->hide();
		set->hide();
		setLabel->hide();
		fil->hide();
		filLabel->hide();
		enabled->hide();
		enabledLm->show();
		bwred->hide();
		bwgreen->hide();
		bwblue->hide();
		rlabel->hide();
		glabel->hide();
		blabel->hide();
		bworan->hide();
		bwyell->hide();
		bwcyan->hide();
		bwmag->hide();
		bwpur->hide();
		orlabel->hide();
		ylabel->hide();
		clabel->hide();
		mlabel->hide();
		plabel->hide();	
		rglabel->show();
		gglabel->show();
		bglabel->show();
		bwredgam->show();
		bwgreengam->show();
		bwbluegam->show();
		Gamlabel->show();	
		}
	 else if(met->get_active_row_number()==1) {
		curveEditorGBW->show();
		curveEditorGBW2->show();
		autoch->hide();
		neutral->hide();
		rglabel->show();
		gglabel->show();
		bglabel->show();
		bwredgam->show();
		bwgreengam->show();
		bwbluegam->show();
		Gamlabel->show();
		set->hide();
		setLabel->hide();
		fil->hide();
		filLabel->hide();
		enabled->hide();
		enabledLm->show();
		curveEditorG->hide();	
		bwred->hide();
		bwgreen->hide();
		bwblue->hide();
		rlabel->hide();
		glabel->hide();
		blabel->hide();
		bworan->hide();
		bwyell->hide();
		bwcyan->hide();
		bwmag->hide();
		bwpur->hide();
		orlabel->hide();
		ylabel->hide();
		clabel->hide();
		mlabel->hide();
		plabel->hide();	
		}
	else {
		autoch->hide();
		neutral->hide();	
		set->hide();
		setLabel->hide();
		fil->hide();
		filLabel->hide();
		enabled->hide();
		enabledLm->show();
		curveEditorG->hide();	
		bwred->hide();
		bwgreen->hide();
		bwblue->hide();
		rlabel->hide();
		glabel->hide();
		blabel->hide();
		rglabel->hide();
		gglabel->hide();
		bglabel->hide();
		bwredgam->hide();
		bwgreengam->hide();
		bwbluegam->hide();		
		Gamlabel->hide();	
		bworan->hide();
		bwyell->hide();
		bwcyan->hide();
		bwmag->hide();
		bwpur->hide();
		orlabel->hide();
		ylabel->hide();
		clabel->hide();
		mlabel->hide();
		plabel->hide();	
		curveEditorGBW->hide();
		curveEditorGBW2->hide();
	 }
	if (listener && (multiImage||enabledLm->get_active())) {
		listener->panelChanged (EvBWmet, met->get_active_text ());
	}
}

void ChMixerbw::enabledLm_toggled () {
	
	if (batchMode) {
		if (enabledLm->get_inconsistent()) {
			enabledLm->set_inconsistent (false);
			enaLmconn.block (true);
			enabledLm->set_active (false);
			enaLmconn.block (false);
		}
		else if (lastEnabledLm)
			enabledLm->set_inconsistent (true);

		lastEnabledLm = enabledLm->get_active ();
	}

	if (listener) {
			listener->panelChanged (EvBWChmixEnabledLm, M("GENERAL_DISABLED"));
			}
			
	}
void ChMixerbw::neutral_pressed () {
// This method deselects auto chmixer 
// and sets "neutral" values to params

    if (batchMode) {
        autoch->set_inconsistent (false);
        autoconn.block (true);
        autoch->set_active (false);
        autoconn.block (false);

        lastAuto = autoch->get_active ();
    }
    else { //!batchMode
        autoch->set_active (false);
        autoch->set_inconsistent (false);
    }
    bwred->setValue(33);
    bwgreen->setValue(33);
    bwblue->setValue(33);
    bworan->setValue(33);
    bwyell->setValue(33);
    bwmag->setValue(33);
    bwpur->setValue(33);
    bwcyan->setValue(33);
 	set->set_active (11);
 	fil->set_active (0);

    listener->panelChanged (EvNeutralBW, M("GENERAL_ENABLED"));
}


void ChMixerbw::enabled_toggled () {
	
	if (batchMode) {
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
		if (enabled->get_active ()){
			listener->panelChanged (EvBWChmixEnabled, M("GENERAL_ENABLED"));
			}
		else {		
			listener->panelChanged (EvBWChmixEnabled, M("GENERAL_DISABLED"));
			}
			
	}
}


void ChMixerbw::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

	bwred->setDefault (defParams->chmixerbw.bwred);
 	bwgreen->setDefault (defParams->chmixerbw.bwgreen);
	bwblue->setDefault (defParams->chmixerbw.bwblue);
	bwredgam->setDefault (defParams->chmixerbw.bwredgam);
 	bwgreengam->setDefault (defParams->chmixerbw.bwgreengam);
	bwbluegam->setDefault (defParams->chmixerbw.bwbluegam);
	bworan->setDefault (defParams->chmixerbw.bworan);
 	bwyell->setDefault (defParams->chmixerbw.bwyell);
	bwcyan->setDefault (defParams->chmixerbw.bwcyan);
	bwmag->setDefault (defParams->chmixerbw.bwmag);
 	bwpur->setDefault (defParams->chmixerbw.bwpur);
   
    if (pedited) {
		bwred->setDefaultEditedState (pedited->chmixerbw.bwred ? Edited : UnEdited);
		bwgreen->setDefaultEditedState (pedited->chmixerbw.bwgreen ? Edited : UnEdited);
		bwblue->setDefaultEditedState (pedited->chmixerbw.bwblue ? Edited : UnEdited);
		bwredgam->setDefaultEditedState (pedited->chmixerbw.bwredgam ? Edited : UnEdited);
		bwgreengam->setDefaultEditedState (pedited->chmixerbw.bwgreengam ? Edited : UnEdited);
		bwbluegam->setDefaultEditedState (pedited->chmixerbw.bwbluegam ? Edited : UnEdited);
		bworan->setDefaultEditedState (pedited->chmixerbw.bworan ? Edited : UnEdited);
		bwyell->setDefaultEditedState (pedited->chmixerbw.bwyell ? Edited : UnEdited);
		bwcyan->setDefaultEditedState (pedited->chmixerbw.bwcyan ? Edited : UnEdited);
		bwmag->setDefaultEditedState (pedited->chmixerbw.bwmag ? Edited : UnEdited);
		bwpur->setDefaultEditedState (pedited->chmixerbw.bwpur ? Edited : UnEdited);
		}
    else {
		bwred->setDefaultEditedState (Irrelevant);
		bwgreen->setDefaultEditedState (Irrelevant);
		bwblue->setDefaultEditedState (Irrelevant);		
		bwredgam->setDefaultEditedState (Irrelevant);
		bwgreengam->setDefaultEditedState (Irrelevant);
		bwbluegam->setDefaultEditedState (Irrelevant);	
		bworan->setDefaultEditedState (Irrelevant);
		bwyell->setDefaultEditedState (Irrelevant);
		bwcyan->setDefaultEditedState (Irrelevant);	
		bwmag->setDefaultEditedState (Irrelevant);
		bwpur->setDefaultEditedState (Irrelevant);
		}
}

void ChMixerbw::autoch_toggled () {

    if (batchMode) {
        if (autoch->get_inconsistent()) {
            autoch->set_inconsistent (false);
            autoconn.block (true);
            autoch->set_active (false);
            autoconn.block (false);
        }
        else if (lastAuto)
            autoch->set_inconsistent (true);

        lastAuto = autoch->get_active ();

        bwred->setEditedState (UnEdited);
        bwgreen->setEditedState (UnEdited);
        bwblue->setEditedState (UnEdited);
        bworan->setEditedState (UnEdited);
        bwyell->setEditedState (UnEdited);
        bwpur->setEditedState (UnEdited);
        bwmag->setEditedState (UnEdited);
        bwcyan->setEditedState (UnEdited);
		
        if (bwred->getAddMode())
            bwred->setValue (33);
        if (bwgreen->getAddMode())
            bwgreen->setValue (33);
        if (bwblue->getAddMode())
            bwblue->setValue (33);
        if (bworan->getAddMode())
            bworan->setValue (33);
        if (bwyell->getAddMode())
            bwyell->setValue (33);
        if (bwmag->getAddMode())
            bwmag->setValue (33);
        if (bwpur->getAddMode())
            bwpur->setValue (33);
        if (bwmag->getAddMode())
            bwpur->setValue (33);
			set->set_active (11);
			fil->set_active (0);

        if (listener) {
            if (!autoch->get_inconsistent()) {
                if (autoch->get_active ()) {
				
                    listener->panelChanged (EvAutoch, M("GENERAL_ENABLED"));}
                else
                    listener->panelChanged (EvFixedch, M("GENERAL_DISABLED"));
            }
        }
    }
    else if (/* !batchMode && */ listener) {
        if (autoch->get_active()) {
			bwred->setValue(33);
			bwgreen->setValue(33);
			bwblue->setValue(33);
			bworan->setValue(33);
			bwyell->setValue(33);
			bwmag->setValue(33);
			bwpur->setValue(33);
			bwcyan->setValue(33);
			set->set_active (11);
			fil->set_active (0);
						
            listener->panelChanged (EvAutoch, M("GENERAL_ENABLED"));
        }
        else {
            listener->panelChanged (EvFixedch, M("GENERAL_DISABLED"));
        }
    }

}


void ChMixerbw::adjusterChanged (Adjuster* a, double newval) {

    if (autoch->get_active() && (a==bwred || a==bwgreen || a==bwblue || a==bworan || a==bwyell || a==bwmag || a==bwpur || a==bwcyan )) {
        autoconn.block(true);
        autoch->set_active (false);
        autoconn.block(false);
        autoch->set_inconsistent (false);
    }

    if (listener  && (multiImage||enabledLm->get_active())) {
		Glib::ustring value = a->getTextValue();		
		if (a == bwred)
			listener->panelChanged (EvBWred,     value );
		else if (a == bwgreen)
			listener->panelChanged (EvBWgreen, value );
		else if (a == bwblue)
			listener->panelChanged (EvBWblue, value );
		else if (a == bwgreengam)
			listener->panelChanged (EvBWgreengam, value );
		else if (a == bwbluegam)
			listener->panelChanged (EvBWbluegam, value );
		else if (a == bwredgam)
			listener->panelChanged (EvBWredgam, value );
		else if (a == bworan)
			listener->panelChanged (EvBWoran, value );
		else if (a == bwyell)
			listener->panelChanged (EvBWyell, value );
		else if (a == bwcyan)
			listener->panelChanged (EvBWcyan, value );
		else if (a == bwmag)
			listener->panelChanged (EvBWmag, value );
		else if (a == bwpur)
			listener->panelChanged (EvBWpur, value );
    }
}

void ChMixerbw::setBatchMode (bool batchMode) {
    removeIfThere (abox, autoch, false);
    autoch = Gtk::manage (new Gtk::CheckButton (M("TP_BWMIX_AUTOCH")));
    autoch->set_tooltip_markup (M("TP_BWMIX_AUTOCH_TIP"));
    autoconn = autoch->signal_toggled().connect( sigc::mem_fun(*this, &ChMixerbw::autoch_toggled) );
    abox->pack_start (*autoch);

    ToolPanel::setBatchMode (batchMode);
	bwred->showEditedCB ();
	bwgreen->showEditedCB ();
	bwblue->showEditedCB ();
	bwredgam->showEditedCB ();
	bwgreengam->showEditedCB ();
	bwbluegam->showEditedCB ();
	bworan->showEditedCB ();
	bwyell->showEditedCB ();
	bwcyan->showEditedCB ();
	bwmag->showEditedCB ();
	bwpur->showEditedCB ();
	met->append_text (M("GENERAL_UNCHANGED"));
	fil->append_text (M("GENERAL_UNCHANGED"));
	set->append_text (M("GENERAL_UNCHANGED"));
    curveEditorG->setBatchMode (batchMode);
    curveEditorGBW->setBatchMode (batchMode);
	toneCurveBW->append_text (M("GENERAL_UNCHANGED"));
    curveEditorGBW2->setBatchMode (batchMode);
	toneCurveBW2->append_text (M("GENERAL_UNCHANGED"));
}

void ChMixerbw::autoOpenCurve () {
    vshape->openIfNonlinear();
    shape->openIfNonlinear();
    shape2->openIfNonlinear();
}


void ChMixerbw::setAdjusterBehavior (bool bwadd, bool bwgadd, bool bwfadd) {

	
	bwred->setAddMode(bwadd);
	bwgreen->setAddMode(bwadd);
	bwblue->setAddMode(bwadd);
	
	bworan->setAddMode(bwfadd);
	bwyell->setAddMode(bwfadd);
	bwcyan->setAddMode(bwfadd);
	bwmag->setAddMode(bwfadd);
	bwpur->setAddMode(bwfadd);

	bwredgam->setAddMode(bwgadd);
	bwgreengam->setAddMode(bwgadd);
	bwbluegam->setAddMode(bwgadd);
	
}

void ChMixerbw::trimValues (rtengine::procparams::ProcParams* pp) {

	bwred->trimValue (pp->chmixerbw.bwred);
	bwgreen->trimValue (pp->chmixerbw.bwgreen);
	bwblue->trimValue (pp->chmixerbw.bwblue);
	bwredgam->trimValue (pp->chmixerbw.bwredgam);
	bwgreengam->trimValue (pp->chmixerbw.bwgreengam);
	bwbluegam->trimValue (pp->chmixerbw.bwbluegam);
	bworan->trimValue (pp->chmixerbw.bworan);
	bwyell->trimValue (pp->chmixerbw.bwyell);
	bwcyan->trimValue (pp->chmixerbw.bwcyan);
	bwmag->trimValue (pp->chmixerbw.bwmag);
	bwpur->trimValue (pp->chmixerbw.bwpur);
}
